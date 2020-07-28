using SDDP

import DelimitedFiles
import Gurobi
import PowerModels
import Statistics

const GRB = Gurobi.Env()

include("data.jl");

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
    rev_line = Dict();
    for l in lines
        rev_line[(line_data[l]["f_bus"], line_data[l]["t_bus"])] = l;
    end

    bus_data = data.network_data["bus"]
    buses = collect(keys(bus_data))

    battery_data = readBattery(data.B);
    batteries = battery_data.IDList;
    battery_bus = battery_data.LocDict;

    disruption_length = 2
    # time period = 15 minutes
    Δt = 0.25

    model = SDDP.PolicyGraph(
        build_graph(data.T, data.ρ - 1, disruption_length),
        sense = :Min,
        lower_bound = 0.0,
        optimizer = () -> Gurobi.Optimizer(GRB),
    ) do subproblem, node
        t, disruption = node.stage, node.disruption_length
        #set_silent(sp)

        # State variables
        @variable(
            subproblem,
            Sp[g in gens],
            SDDP.State,
            initial_value = 0.0
        )
        @constraint(subproblem, spUb[g in gens], Sp[g].out <= gen_data[g]["pmax"]);
        @constraint(subproblem, spLb[g in gens], Sp[g].out >= gen_data[g]["pmin"]);
        @variable(
            subproblem,
            Sq[g in gens],
            SDDP.State,
            initial_value = 0.0
        )
        @constraint(subproblem, sqUb[g in gens], Sq[g].out <= gen_data[g]["qmax"]);
        @constraint(subproblem, sqLb[g in gens], Sq[g].out >= gen_data[g]["qmin"]);
        @variable(
            subproblem,
            0 <= w[i in batteries] <= battery_data.capacity[i],
            SDDP.State,
            initial_value = 0.0
        )
        @variable(
            subproblem,
            0 <= u[i in batteries] <= battery_data.uCap[i],
            SDDP.State,
            initial_value = 0.0
        )

        if t > 1
            # u is a first-stage decision only
            @constraint(subproblem, [b = batteries], u[b].out == u[b].in)
        end

        # Control variables
        @variables(subproblem, begin
            bus_data[i]["vmin"]^2 <= V[i = buses] <= bus_data[i]["vmax"]^2
            P[lines]
            Q[lines]
            W[lines]
            zp[batteries]
            zq[batteries]
            y[batteries]
            L[buses, [:p, :q], [:+, :-]] >= 0
        end)

        # battery inventory transition
        @constraint(
            subproblem,
            [b in batteries],
            w[b].out == w[b].in - y[b]*Δt
        )

        # Kirchoff's Current Law
        @constraint(
            subproblem,
            [i in buses],
            sum(zp[b] for b in batteries if battery_bus[b] == i) +
                sum(Sp[g].out for g in gens if gen_data[g]["gen_bus"] == i) -
                data.P[parse(Int, i), t] +
                L[i, :p, :+] - L[i, :p, :-] ==
                sum(P[l] for l in lines if from(l) == i) -
                sum(P[l] for l in lines if to(l) == i)
        )
        @constraint(
            subproblem,
            [i in buses],
            sum(zq[b] for b in batteries if battery_bus[b] == i) +
                sum(Sq[g].out for g in gens if gen_data[g]["gen_bus"] == i) -
                data.Q[parse(Int, i), t] +
                L[i, :q, :+] - L[i, :q, :-] ==
                sum(Q[l] for l in lines if from(l) == i) -
                sum(Q[l] for l in lines if to(l) == i)
        )

        # LinDist Flow power flow equations
        @constraint(subproblem, linDistFlow_ub[l = lines], V[from(l)] - V[to(l)] - 2 * (
            line_data[l]["br_r"] * P[l] + line_data[l]["br_x"] * Q[l]) >= 0)
        @constraint(subproblem, linDistFlow_lb[l = lines], V[from(l)] - V[to(l)] - 2 * (
            line_data[l]["br_r"] * P[l] + line_data[l]["br_x"] * Q[l]) <= 0)

        @constraints(subproblem, begin
            # Thermal limits
            [l = lines], [W[l], P[l], Q[l]] in SecondOrderCone()
            [b = batteries], [u[b].out, zp[b], zq[b]] in SecondOrderCone()
            # battery efficiency curve
            [b = batteries, l in 1:length(battery_data.ηα[b])], zp[b] <= battery_data.ηα[b][l]*y[b] + battery_data.ηβ[b][l]
        end)
        # thermal limits under the normal condition
        @constraint(subproblem, thermal[l = lines], W[l] == line_data[l]["rate_a"]);

        # Battery efficiency
        # Generator ramp up limit
        @constraints(subproblem, begin
            genRamp_up[g = gens], Sp[g].out - Sp[g].in <= gen_data[g]["ramp_agc"]
            genRamp_down[g = gens], Sp[g].out - Sp[g].in >= -gen_data[g]["ramp_agc"]
        end)

        # Generator ramp down limit


        @stageobjective(subproblem, 0.0)
        SDDP.parameterize(subproblem, data.support, data.nominal_probability) do ω
            if ω isa Int
                # Change generator bounds
                for g in gens
                    if g == ω
                        JuMP.set_normalized_rhs(spUb[ω], 0)
                        JuMP.set_normalized_rhs(sqUb[ω], 0)
                        JuMP.set_normalized_rhs(spLb[ω], 0)
                        JuMP.set_normalized_rhs(sqLb[ω], 0)
                        JuMP.set_normalized_rhs(genRamp_up[ω], 1000*gen_data[g]["ramp_agc"])
                        JuMP.set_normalized_rhs(genRamp_down[ω], -1000*gen_data[g]["ramp_agc"])
                    else
                        JuMP.set_normalized_rhs(spUb[ω], gen_data[g]["pmax"])
                        JuMP.set_normalized_rhs(sqUb[ω], gen_data[g]["qmax"])
                        JuMP.set_normalized_rhs(spLb[ω], gen_data[g]["pmin"])
                        JuMP.set_normalized_rhs(sqLb[ω], gen_data[g]["qmin"])
                        JuMP.set_normalized_rhs(genRamp_up[ω], gen_data[g]["ramp_agc"])
                        JuMP.set_normalized_rhs(genRamp_down[ω], -gen_data[g]["ramp_agc"])
                    end
                end
            else
                # Change line bounds
                if (ω.first,ω.second) in keys(rev_line)
                    line_down = rev_line((ω.first,ω.second));
                else
                    line_down = rev_line((ω.second,ω.first));
                end
                for l in lines
                    if l == line_down
                        JuMP.set_normalized_rhs(linDistFlow_ub[line_down], -1000*line_data[l]["rate_a"]);
                        JuMP.set_normalized_rhs(linDistFlow_lb[line_down], 1000*line_data[l]["rate_a"]);
                        JuMP.set_normalized_rhs(thermal[l], 0.0);
                    else
                        JuMP.set_normalized_rhs(linDistFlow_ub[line_down], 0);
                        JuMP.set_normalized_rhs(linDistFlow_lb[line_down], 0);
                        JuMP.set_normalized_rhs(thermal[l], line_data[l]["rate_a"]);
                    end
                end
            end
        end
    end
    return model
end
