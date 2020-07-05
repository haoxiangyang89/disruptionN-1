using SDDP

import DelimitedFiles
import Gurobi
import PowerModels
import Statistics

const GRB = Gurobi.Env()

include("data.jl")

reset_limits(pm, ::Nothing, ::Nothing) = nothing

function reset_limits(pm, reset_key, reset_data)
    return set_bounds(pm, reset_key, reset_data...)
end

function set_bounds(pm, ω::Int, pgl, pgu, qgl, qgu)
    pg = PowerModels.var(pm, :pg)[ω]
    gq = PowerModels.var(pm, :qg)[ω]
    previous = (
        lower_bound(pg),
        upper_bound(pg),
        lower_bound(qg),
        upper_bound(qg),
    )
    set_lower_bound(pg, pgl)
    set_upper_bound(pg, pgu)
    set_lower_bound(qg, qgl)
    set_upper_bound(qg, qgu)
    return previous
end

function set_bounds(pm, ω::Int, pgl, pgu, qgl, qgu)
    for key in keys(PowerModels.var(pm, :p))
        l, i, j = key[1]
        if (ω != i => j && ω != j => i)
            continue
        end
        p = PowerModels.var(pm, :p)[key]
        q = PowerModels.var(pm, :q)[key]


        fix(, 0.0; force = true)
        fix(PowerModels.var(pm, :q)[key], 0.0; force = true)
    end

    pg = PowerModels.var(pm, :pg)[ω]
    gq = PowerModels.var(pm, :qg)[ω]
    previous = (
        lower_bound(pg),
        upper_bound(pg),
        lower_bound(qg),
        upper_bound(qg),
    )
    set_lower_bound(pg, pgl)
    set_upper_bound(pg, pgu)
    set_lower_bound(qg, qgl)
    set_upper_bound(qg, qgu)
    return previous
end

function main()
    data = build_data(
        "../data/case13_ieee.m",
        "../data/testDataB_13.csv",
        "../data/testDataP_13.csv",
        "../data/testDataQ_13.csv",
        "../data/testProbRead_13.csv",
    )
    disruption_length = 2
    model = SDDP.PolicyGraph(
        build_graph(data.T, data.ρ, disruption_length),
        sense = :Min,
        lower_bound = 0.0,
        optimizer = () -> Gurobi.Optimizer(GRB),
    ) do sp, node
        stage, disruption = node.stage, node.disruption_length
        set_silent(sp)

        # Modify the pd and qd for the buses
        for load in values(data.network_data["load"])
            load["pq"] = data.P[load["load_bus"], stage]
            load["gq"] = data.Q[load["load_bus"], stage]
        end

        # TODO: account for the fact that the disruption periods correspond to
        # `disruption` time periods.
        pm = PowerModels.instantiate_model(
            data.network_data,
            PowerModels.SOCWRPowerModel,
            PowerModels.build_opf;
            jump_model = sp
        )
        # TODO: state variables
        generators = collect(keys(data.network_data["gen"]))
        @variable(
            sp, s[generators, [:P, :Q]] >= 0, SDDP.State, initial_value = 0.0
        )

        @stageobjective(sp, objective_function(sp))

        if disruption > 0
            reset_key, reset_data = nothing, nothing
            SDDP.parameterize(sp, data.support, data.nominal_probability) do ω
                # TODO: Reset the line and generator limits because the previous
                # scenario will still be set!
                reset_limits(pm, reset_key, reset_data)
                if ω isa Pair
                    # Enforce the line disruption

                else
                    @assert ω isa Int
                    # Enforce the generator disruption
                    reset_key = ω
                    reset_data = set_bounds(pm, ω, 0.0, 0.0, 0.0, 0.0)
                end
            end
        end
    end
    return model
end
