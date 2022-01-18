using PoingrSimulator
using Test

using LinearAlgebra # norm
using CSV           # history file
using JLD2          # snapshots

# vtk
using ReadVTK
using NaturalSort

const fix_results = false

function check_results(inputtoml::String)
    @assert endswith(inputtoml, ".toml")
    testcase = joinpath(basename(dirname(inputtoml)), basename(inputtoml))
    println("\n>> ", testcase)
    @testset "$testcase" begin
        @time PoingrSimulator.main(inputtoml)

        INPUT = PoingrSimulator.build_INPUT(inputtoml)
        proj_dir = INPUT.General.project_directory
        output_dir = INPUT.Output.directory

        testname = first(splitext(basename(inputtoml)))
        if haskey(INPUT.General, :restart)
            testname = replace(testname, "_restart" => "")
        end

        # vtk files
        vtk_file = joinpath(
            "paraview",
            sort(
                filter(
                    file -> endswith(file, "_1.vtu"),
                    only(walkdir(joinpath(output_dir, "paraview")))[3]
                ),
                lt = natural
            )[end],
        )

        if fix_results && !haskey(INPUT.General, :restart)
            cp(joinpath(output_dir, vtk_file),
               joinpath(proj_dir, "output", "$testname.vtu"); force = true)
        else
            # check results
            expected = VTKFile(joinpath(proj_dir, "output", "$testname.vtu")) # expected output
            result = VTKFile(joinpath(output_dir, vtk_file))
            expected_points = get_points(expected)
            result_points = get_points(result)
            @assert size(expected_points) == size(result_points)
            @test all(eachindex(expected_points)) do i
                norm(expected_points[i] - result_points[i]) < 0.05*INPUT.General.grid_space
            end
        end

        # snapshots file
        if !fix_results && !haskey(INPUT.General, :restart)
            nsteps = floor(Int, INPUT.General.total_time / INPUT.Output.interval)
            jldopen(joinpath(output_dir, "snapshots.jld2"), "r") do file
                @test keys(file) == string.(0:nsteps)
                for i in keys(file)
                    @test file[i] isa NamedTuple
                end
            end
        end

        # test history.csv if exists
        if isfile(joinpath(output_dir, "history.csv")) && !haskey(INPUT.General, :restart)
            if fix_results
                cp(joinpath(output_dir, "history.csv"),
                   joinpath(proj_dir, "output", "$testname.csv"); force = true)
            else
                # check results
                output = CSV.File(joinpath(proj_dir, "output", "$testname.csv")) # expected output
                history = CSV.File(joinpath(output_dir, "history.csv"))
                for name in propertynames(output)
                    output_col = output[name]
                    history_col = history[name]
                    @test output_col ≈ history_col atol=1e-8 rtol=1e-3
                end
            end
        end
    end
end

@testset "$module_name" for module_name in ("PenetrateIntoGround", "FreeRun")
    for (root, dirs, files) in walkdir(module_name)
        for file in files
            path = joinpath(root, file)
            splitext(dirname(path))[2] != ".tmp" && endswith(path, ".toml") && check_results(path)
        end
    end
end
