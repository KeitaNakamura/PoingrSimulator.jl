module PenetrateIntoGround

using PoingrSimulator
using PoingrSimulator: Input, Input_Phase
using Poingr
using GeometricObjects

using Serialization

function preprocess_input!(input::Input)
    input.Material = input.SoilLayer
    input.BoundaryCondition.left   = ContactMohrCoulomb(:slip)
    input.BoundaryCondition.right  = ContactMohrCoulomb(:slip)
    input.BoundaryCondition.bottom = ContactMohrCoulomb(:sticky)
    input.BoundaryCondition.top    = ContactMohrCoulomb(:slip)
    for rigidbody in input.RigidBody
        @assert rigidbody.control == true
        # for calculation of effective mass in collision
        rigidbody.density = Inf
        rigidbody.model.m = Inf
    end
    @assert isempty(input.BoundaryCondition.Dirichlet)
end

function initialize(input::Input)
    NodeState = @NamedTuple begin
        m::Float64
        m′::Float64
        v::Vec{2, Float64}
        v_n::Vec{2, Float64}
        m_contacted::Float64
        vᵣ::Vec{2, Float64}
        fc::Vec{2, Float64}
        d::Vec{2, Float64}
        μ::Vec{2, Float64} # [μ, c]
        poly_coef::Vec{3, Float64}
        poly_mat::Mat{3, 3, Float64, 9}
    end
    L = isa(input.General.transfer, LinearWLS) ? 3 : 2
    PointState = @NamedTuple begin
        m::Float64
        V::Float64
        x::Vec{2, Float64}
        x0::Vec{2, Float64}
        v::Vec{2, Float64}
        b::Vec{2, Float64}
        σ::SymmetricSecondOrderTensor{3, Float64, 6}
        ϵ::SymmetricSecondOrderTensor{3, Float64, 6}
        F::SecondOrderTensor{3, Float64, 9}
        J::Float64
        ∇v::SecondOrderTensor{3, Float64, 9}
        P::Float64
        C::Mat{2, L, Float64, 2*L}
        r::Vec{2, Float64}
        index::Int
        matindex::Int
    end

    # General
    coordinate_system = input.General.coordinate_system
    (xmin, xmax), (ymin, ymax) = input.General.domain
    dx = input.General.grid_space
    g = input.General.gravity

    # SoilLayer
    soillayers = input.SoilLayer
    H = sum(layer -> layer.thickness, soillayers) # ground surface
    @assert ymin + H ≤ ymax

    # Advanced
    α = input.Advanced.contact_threshold_scale
    nptsincell = input.Advanced.npoints_in_cell

    grid = Grid(NodeState, input.General.interpolation, xmin:dx:xmax, ymin:dx:ymax; coordinate_system)
    pointstate = generate_pointstate((x,y) -> y < ymin + H, PointState, grid; n = nptsincell)
    rigidbody = only(input.RigidBody).model

    bottom = ymin
    for i in length(soillayers):-1:1 # from low to high
        layer = soillayers[i]
        Threads.@threads for p in 1:length(pointstate)
            y = pointstate.x[p][2]
            if bottom ≤ y ≤ bottom + layer.thickness
                pointstate.matindex[p] = i
            end
        end
        bottom += layer.thickness
    end

    # initialize variables of points
    Threads.@threads for p in 1:length(pointstate)
        layerindex = pointstate.matindex[p]
        σ_y = 0.0
        for layer in soillayers[begin:layerindex-1]
            ρ₀ = layer.density
            σ_y += -ρ₀ * g * layer.thickness
        end
        h = sum(layer -> layer.thickness, soillayers[layerindex:end])
        layer = soillayers[layerindex]
        y = pointstate.x[p][2]
        ρ0 = layer.density
        K0 = layer.K0
        σ_y += -ρ0 * g * (h - y)
        σ_x = K0 * σ_y
        pointstate.σ[p] = (@Mat [σ_x 0.0 0.0
                                 0.0 σ_y 0.0
                                 0.0 0.0 σ_x]) |> symmetric
        pointstate.m[p] = ρ0 * pointstate.V[p]
    end
    @. pointstate.x0 = pointstate.x
    @. pointstate.b = Vec(0.0, -g)

    if only(input.RigidBody).reset_position
        y0 = minimum(x -> x[2], coordinates(rigidbody))
        δ = sqrt(eps(Float64))
        translate!(rigidbody, Vec(0.0, (ymin - y0) + H + (α-1)*(dx/nptsincell)/2 + δ))
    else
        PoingrSimulator.remove_invalid_pointstate!(pointstate, input)
    end

    t = 0.0
    t, grid, pointstate, rigidbody, deepcopy(rigidbody)
end

function main(input::Input, phase::Input_Phase, t, grid, pointstate, rigidbody, rigidbody0)

    # General/Output
    dx = input.General.grid_space
    t_start = t
    t_stop = phase.time_stop
    t_step = input.Output.time_interval

    # SoilLayer
    matmodels = map(x -> x.model, input.SoilLayer)

    ################
    # Output files #
    ################

    outdir = input.Output.directory
    outputs = Dict{String, Any}()
    if input.Output.paraview
        mkpath(joinpath(outdir, "paraview"))
        outputs["paraview_file"] = joinpath(outdir, "paraview", "output")
        paraview_collection(vtk_save, outputs["paraview_file"])
    end
    if input.Output.history
        outputs["history_file"] = joinpath(outdir, "history.csv")
        open(outputs["history_file"], "w") do io
            write(io, "depth,force\n")
        end
    end
    if input.Output.snapshots || input.Output.snapshot_last
        mkpath(joinpath(outdir, "snapshots"))
    end
    if isdefined(input.Injection, :main_output)
        input.Injection.main_output_initialize((;
            input,
            t,
            grid,
            pointstate,
            rigidbody,
            rigidbody0,
        ))
    end

    ##################
    # Run simulation #
    ##################

    cache = MPCache(grid, pointstate.x)
    logger = Logger(t_start, t_stop, t_step; input.General.showprogress)
    update!(logger, t)
    writeoutput(outputs, input, grid, pointstate, rigidbody, rigidbody0, t, logindex(logger))

    try
        while !isfinised(logger, t)
            dt = phase.CFL * PoingrSimulator.safe_minimum(pointstate) do pt
                PoingrSimulator.timestep(matmodels[pt.matindex], pt, dx)
            end
            PoingrSimulator.advancestep!(grid, pointstate, [rigidbody], cache, dt, input, phase)

            if input.Output.quickview
                update!(logger, t += dt; print = PoingrSimulator.quickview_sparsity_pattern(cache.spat))
            else
                update!(logger, t += dt)
            end

            if islogpoint(logger)
                if input.Advanced.reorder_pointstate
                    Poingr.reorder_pointstate!(pointstate, cache)
                end
                writeoutput(outputs, input, grid, pointstate, rigidbody, rigidbody0, t, logindex(logger))
            end
        end
    catch e
        writeoutput(outputs, input, grid, pointstate, rigidbody, rigidbody0, t, "error")
        rethrow()
    end

    if input.Output.snapshot_last
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot_last"),
            (; t, grid, pointstate, rigidbody, rigidbody0)
        )
    end

    t
end

function writeoutput(
        outputs::Dict{String, Any},
        input::Input,
        grid::Grid,
        pointstate::AbstractVector,
        rigidbody::GeometricObject,
        rigidbody0::GeometricObject,
        t::Real,
        output_index,
    )
    if input.Output.paraview
        compress = true
        paraview_file = outputs["paraview_file"]
        paraview_collection(paraview_file, append = true) do pvd
            vtk_multiblock(string(paraview_file, output_index)) do vtm
                vtk_points(vtm, pointstate.x; compress) do vtk
                    PoingrSimulator.write_vtk_points(vtk, pointstate)
                end
                vtk_grid(vtm, rigidbody)
                if input.Output.paraview_grid
                    vtk_grid(vtm, grid; compress) do vtk
                        vtk["nodal contact force"] = vec(grid.state.fc)
                        vtk["nodal contact distance"] = vec(grid.state.d)
                        vtk["nodal friction"] = vec(grid.state.μ[1])
                    end
                end
                pvd[t] = vtm
            end
        end
    end

    if input.Output.history
        history_file = outputs["history_file"]
        open(history_file, "a") do io
            (_, _), (ymin, _) = input.General.domain
            H = ymin + sum(layer -> layer.thickness, input.SoilLayer) # ground surface
            tip = minimum(x -> x[2], coordinates(rigidbody))
            depth = H - tip
            force = -sum(grid.state.fc)[2]
            if input.General.coordinate_system isa Axisymmetric
                force *= 2π
            end
            write(io, "$depth,$force\n")
        end
    end

    if input.Output.snapshots
        serialize(
            joinpath(input.Output.directory, "snapshots", "snapshot$output_index"),
            (; t, grid, pointstate, rigidbody, rigidbody0)
        )
    end

    if isdefined(input.Injection, :main_output)
        input.Injection.main_output((;
            input,
            grid,
            pointstate,
            rigidbody,
            rigidbody0,
            t,
            output_index,
        ))
    end
end

end
