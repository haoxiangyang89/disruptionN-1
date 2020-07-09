using SDDP

import DelimitedFiles
import Gurobi
import PowerModels
import Statistics

const GRB = Gurobi.Env()

include("data.jl")

function main()
    data = build_data(
        "../data/case13_ieee.m",
        "../data/testDataB_13.csv",
        "../data/testDataP_13.csv",
        "../data/testDataQ_13.csv",
        "../data/testProbRead_13.csv",
    )
    gen_data = data.network_data["gen"]
    gens = collect(keys(gen_data))

    line_data = data.network_data["branch"]
    lines = collect(keys(line_data))
    from(l) = string(line_data[l]["f_bus"])
    to(l) = string(line_data[l]["t_bus"])

    buses = 1:size(data.P, 1)

    batteries = 1:size(data.B, 1)
    battery_bus = Dict(string(Int, data.B[b, 2]) for b in batteries)

    disruption_length = 2
    model = SDDP.PolicyGraph(
        build_graph(data.T, data.ρ - 1, disruption_length),
        sense = :Min,
        lower_bound = 0.0,
        optimizer = () -> Gurobi.Optimizer(GRB),
    ) do subproblem, node
        t, disruption = node.stage, node.disruption_length
        set_silent(sp)

        # State variables
        @variable(
            subproblem,
            gen_data[g]["pmin"] <= Sp[g = gens] <= gen_data[g]["pmax"],
            SDDP.State,
            initial_value = 0.0
        )
        @variable(
            subproblem,
            gen_data[g]["qmin"] <= Sq[g = gens] <= gen_data[g]["qmax"],
            SDDP.State,
            initial_value = 0.0
        )
        @variable(subproblem, w[batteries], SDDP.State, initial_value = 0.0)
        @variable(subproblem, u[batteries], SDDP.State, initial_value = 0.0)

        if t > 1
            # u is a first-stage decision only
            @constraint(subproblem, [b = batteries], u[b].out == u[b].in)
        end

        # Control variables
        @variables(subproblem, begin
            V[buses]
            P[lines]
            Q[lines]
            W[lines]
            zp[batteries]
            zq[batteries]
            L[buses, [:p, :q], [:+, :-]] >= 0
        end)

        # Kirchoff's Current Law
        @constraint(
            subproblem,
            [i in buses],
            sum(zp[b] for b in batteries if battery_bus[b] == i) +
                sum(Sp[g].out for g in gens if gen_data[g]["gen_bus"] == i) -
                data.P[parse(Int, i), t] +
                L[i, :p, :+] - L[i, :p, :-] ==
                sum(P[l, t] for l in lines if from(l) == i)
        )
        @constraint(
            subproblem,
            [i in buses],
            sum(zq[b] for b in batteries if battery_bus[b] == i) +
                sum(Sq[g].out for g in gens if gen_data[g]["gen_bus"] == i) -
                data.Q[parse(Int, i), t] +
                L[i, :q, :+] - L[i, :q, :-] ==
                sum(Q[l, t] for l in lines if from(l) == i)
        )

        @constraints(subproblem, begin
            # LinDist Flow power flow equations
            [l = lines], V[to(l)] == V[from(l)] - 2 * (
                line_data[l]["br_r"] * P[l] + line_data[l]["br_x"] * Q[l]
            )
            # Thermal limits
            [l = lines], [W[l], P[l], Q[l]] in SecondOrderCone()
            [b = batteries], [u.out[b], zp[b], zq[b]] in SecondOrderCone()
            # Generator ramp up limit
            [g = gens], s[g, :p].out - s[g, :p].in <= gen_data[g]["ramp_10"]
            # Generator ramp down limit
            [g = gens], s[g, :p].out - s[g, :p].in >= -gen_data[g]["ramp_10"]
        end)

        # Battery efficiency

        @stageobjective(sp, 0.0)
        SDDP.parameterize(sp, data.support, data.nominal_probability) do ω
            if ω isa Int
                # Change generator bounds
            else
                # Change line bounds
            end
        end
    end
    return model
end
