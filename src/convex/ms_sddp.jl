# sddp with every scenario possible to occur at every period
using SDDP;

function create_sddp(T, Δt, support, nominal_probability)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end

    model = SDDP.LinearPolicyGraph(
        stages = T,
        sense = :Min,
        lower_bound = 0.0,
        optimizer = optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "QCPDual" => 1, "NumericFocus" => 3),
    ) do subproblem, node
        t = node

        # State variables
        @variable(
            subproblem,
            sp[g in fData.genIDList],
            SDDP.State,
            initial_value = 0.0
        )
        @constraint(subproblem, spUb[g in fData.genIDList], sp[g].out <= fData.Pmax[g]);
        @constraint(subproblem, spLb[g in fData.genIDList], sp[g].out >= fData.Pmin[g]);
        @variable(
            subproblem,
            0 <= w[i in bData.IDList] <= bData.cap[i],
            SDDP.State,
            initial_value = 0.0
        )
        @variable(
            subproblem,
            0 <= u[i in bData.IDList] <= bData.uCap[i],
            SDDP.State,
            initial_value = 0.0
        )

        if t > 1
            # u is a first-stage decision only
            @constraint(subproblem, [b in bData.IDList], u[b].out == u[b].in)
        end

        # Control variables
        @variables(subproblem, begin
            fData.Vmin[i]^2 <= V[i in fData.IDList] <= fData.Vmax[i]^2
            sq[fData.genIDList]
            P[fData.brList]
            Q[fData.brList]
            W[fData.brList]
            zp[bData.IDList]
            zq[bData.IDList]
            y[bData.IDList]
            L[fData.IDList, [:p, :q], [:+, :-]] >= 0
        end)

        @constraint(subproblem, sqUb[g in fData.genIDList], sq[g] <= fData.Qmax[g]);
        @constraint(subproblem, sqLb[g in fData.genIDList], sq[g] >= fData.Qmin[g]);

        # battery inventory transition
        @constraint(
            subproblem,
            [b in bData.IDList],
            w[b].out == w[b].in - y[b]*Δt
        )

        sphatsum = Dict();
        for i in fData.IDList
            sphatsum[i] = @expression(subproblem, 0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i] += sp[j].out;
                end
            end
        end
        sqhatsum = Dict();
        for i in fData.IDList
            sqhatsum[i] = @expression(subproblem, 0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i] += sq[j];
                end
            end
        end
        # Kirchoff's Current Law
        @constraint(
            subproblem,
            [i in fData.IDList],
            sum(zp[b] for b in bData.IDList if bData.Loc[b] == i) +
                sphatsum[i] -
                dData.pd[i][t] +
                L[i, :p, :+] - L[i, :p, :-] ==
                sum(P[l] for l in fData.brList if l[1] == i)
        )
        @constraint(
            subproblem,
            [i in fData.IDList],
            sum(zq[b] for b in bData.IDList if bData.Loc[b] == i) +
                sqhatsum[i] -
                dData.qd[i][t] +
                L[i, :q, :+] - L[i, :q, :-] ==
                sum(Q[l] for l in fData.brList if l[1] == i)
        )

        # LinDist Flow power flow equations
        @constraint(subproblem, linDistFlow_ub[l in fData.brList], V[l[1]] - V[l[2]] - 2 * (
            Rdict[l] * P[l] + Xdict[l] * Q[l]) >= 0)
        @constraint(subproblem, linDistFlow_lb[l in fData.brList], V[l[1]] - V[l[2]] - 2 * (
            Rdict[l] * P[l] + Xdict[l] * Q[l]) <= 0)

        # reverse flow consistency
        @constraint(subproblem, pequal[k in fData.brList], P[k] == -P[(k[2],k[1],k[3])]);
        @constraint(subproblem, qequal[k in fData.brList], Q[k] == -Q[(k[2],k[1],k[3])]);

        @constraints(subproblem, begin
            # Thermal limits
            [l in fData.brList], [W[l]; [P[l], Q[l]]] in SecondOrderCone()
            [b in bData.IDList], [u[b].out; [zp[b], zq[b]]] in SecondOrderCone()
            # battery efficiency curve
            [b in bData.IDList, l in 1:length(bData.ηα[b])], zp[b] <= bData.ηα[b][l]*y[b] + bData.ηβ[b][l]
        end)
        # thermal limits under the normal condition
        @constraint(subproblem, thermal[l in fData.brList], W[l] == fData.rateA[l]);

        @variable(subproblem, Pp[l in fData.brList] >= 0);
        @variable(subproblem, Qp[l in fData.brList] >= 0);
        @constraint(subproblem, Pp1[l in fData.brList], Pp[l] >= P[l]);
        @constraint(subproblem, Pp2[l in fData.brList], Pp[l] >= -P[l]);
        @constraint(subproblem, Qp1[l in fData.brList], Qp[l] >= Q[l]);
        @constraint(subproblem, Qp2[l in fData.brList], Qp[l] >= -Q[l]);
        @constraint(subproblem, addCons[l in fData.brList], Pp[l] + Qp[l] >= 0.00001);

        @variable(subproblem, zpp[b in bData.IDList] >= 0);
        @variable(subproblem, zqp[b in bData.IDList] >= 0);
        @constraint(subproblem, zpp1[b in bData.IDList], zpp[b] >= zp[b]);
        @constraint(subproblem, zpp2[b in bData.IDList], zpp[b] >= -zp[b]);
        @constraint(subproblem, zqp1[b in bData.IDList], zqp[b] >= zq[b]);
        @constraint(subproblem, zqp2[b in bData.IDList], zqp[b] >= -zq[b]);
        @constraint(subproblem, addConsz[b in bData.IDList], zpp[b] + zqp[b] >= 0.00001);

        # Battery efficiency
        # Generator ramp up limit
        @constraints(subproblem, begin
            genRamp_up[g in fData.genIDList], sp[g].out - sp[g].in <= fData.RU[g]
            genRamp_down[g in fData.genIDList], sp[g].out - sp[g].in >= fData.RD[g]
        end)

        # calculate the cost
        # @stageobjective(subproblem, 0.0);
        if t == 1
            @stageobjective(subproblem, sum(bData.cost[i]*u[i].out for i in bData.IDList) +
                sum(fData.cz*(L[i, :p, :+] + L[i, :p, :-] + L[i, :q, :+] + L[i, :q, :-]) for i in fData.IDList) +
                sum(fData.cp[g].params[1]*sp[g].out^2 + fData.cp[g].params[2]*sp[g].out for g in fData.genIDList if fData.cp[g].n == 3) +
                sum(fData.cp[g].params[1]*sp[g].out for g in fData.genIDList if fData.cp[g].n == 2)
                );
        else
            @stageobjective(subproblem, sum(fData.cz*(L[i, :p, :+] + L[i, :p, :-] + L[i, :q, :+] + L[i, :q, :-]) for i in fData.IDList) +
                sum(fData.cp[g].params[1]*sp[g].out^2 + fData.cp[g].params[2]*sp[g].out for g in fData.genIDList if fData.cp[g].n == 3) +
                sum(fData.cp[g].params[1]*sp[g].out for g in fData.genIDList if fData.cp[g].n == 2)
                );
        end
        SDDP.parameterize(subproblem, support, nominal_probability) do ω
            if ω isa Int
                # Change generator bounds
                for g in fData.genIDList
                    if g == ω
                        JuMP.set_normalized_rhs(spUb[ω], 0.0)
                        JuMP.set_normalized_rhs(sqUb[ω], 0.0)
                        JuMP.set_normalized_rhs(spLb[ω], 0.0)
                        JuMP.set_normalized_rhs(sqLb[ω], 0.0)
                        JuMP.set_normalized_rhs(genRamp_up[ω], 1000*fData.RU[g])
                        JuMP.set_normalized_rhs(genRamp_down[ω], 1000*fData.RD[g])
                    else
                        JuMP.set_normalized_rhs(spUb[g], fData.Pmax[g])
                        JuMP.set_normalized_rhs(sqUb[g], fData.Qmax[g])
                        JuMP.set_normalized_rhs(spLb[g], fData.Pmin[g])
                        JuMP.set_normalized_rhs(sqLb[g], fData.Qmin[g])
                        JuMP.set_normalized_rhs(genRamp_up[g], fData.RU[g])
                        JuMP.set_normalized_rhs(genRamp_down[g], fData.RD[g])
                    end
                end
                for l in fData.brList
                    JuMP.set_normalized_rhs(linDistFlow_ub[l], 0.0);
                    JuMP.set_normalized_rhs(linDistFlow_lb[l], 0.0);
                    JuMP.set_normalized_rhs(thermal[l], fData.rateA[l]);
                end
            else
                # Change line bounds
                for l in fData.brList
                    if ((l[1] == ω[1])&&(l[2] == ω[2]))||((l[1] == ω[2])&&(l[2] == ω[1]))
                        JuMP.set_normalized_rhs(linDistFlow_ub[l], -1000*fData.rateA[l]);
                        JuMP.set_normalized_rhs(linDistFlow_lb[l], 1000*fData.rateA[l]);
                        JuMP.set_normalized_rhs(thermal[l], 0.0);
                    else
                        JuMP.set_normalized_rhs(linDistFlow_ub[l], 0.0);
                        JuMP.set_normalized_rhs(linDistFlow_lb[l], 0.0);
                        JuMP.set_normalized_rhs(thermal[l], fData.rateA[l]);
                    end
                end
                for g in fData.genIDList
                    JuMP.set_normalized_rhs(spUb[g], fData.Pmax[g])
                    JuMP.set_normalized_rhs(sqUb[g], fData.Qmax[g])
                    JuMP.set_normalized_rhs(spLb[g], fData.Pmin[g])
                    JuMP.set_normalized_rhs(sqLb[g], fData.Qmin[g])
                    JuMP.set_normalized_rhs(genRamp_up[g], fData.RU[g])
                    JuMP.set_normalized_rhs(genRamp_down[g], fData.RD[g])
                end
            end
        end
    end

    return model;
end

include("loadMod.jl");
const GUROBI_ENV = Gurobi.Env();
Δt = 0.25;
N = 5;
T = 8;
τ = 4;
caseList = [13,33,123];
readInData(ci,caseList,T,τ);
model = create_sddp(T,Δt,[1,3,(2,7)],[1/3,1/3,1/3]);
SDDP.train(model; iteration_limit = 100, duality_handler = SDDP.ContinuousConicDuality())
