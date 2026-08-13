#!/usr/bin/env julia

include(joinpath(@__DIR__, "src", "GraphSplit.jl"))
using .GraphSplit

function usage()
    println("Usage: julia --project=. build_travel_times.jl [graphsplit.toml] [--force]")
end

function main(args=ARGS)
    config_file = "graphsplit.toml"
    force = false
    for argument in args
        if argument == "--force"
            force = true
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
    table = GraphSplit.build_travel_times(cfg; force=force)
    table isa GraphSplit.TravelTimeTable && close(table.io)
    return 0
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main())
