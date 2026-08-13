#!/usr/bin/env julia

include(joinpath(@__DIR__, "src", "GraphSplit.jl"))
using .GraphSplit

function usage()
    println("Usage: julia --project=. run_graphsplit.jl [graphsplit.toml] [--build-tt-only] [--force-rebuild-tt]")
end

function main(args=ARGS)
    config_file = "graphsplit.toml"
    build_only = false
    force_rebuild = false
    for argument in args
        if argument == "--build-tt-only"
            build_only = true
        elseif argument == "--force-rebuild-tt"
            build_only = true
            force_rebuild = true
        elseif argument in ("-h", "--help")
            usage()
            return 0
        elseif endswith(lowercase(argument), ".toml")
            config_file = argument
        else
            usage()
            error("Unknown argument: $argument")
        end
    end

    cfg = GraphSplit.load_config(config_file)
    if build_only
        table = GraphSplit.build_travel_times(cfg; force=force_rebuild)
        table isa GraphSplit.TravelTimeTable && close(table.io)
    else
        GraphSplit.run(cfg)
    end
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main())
