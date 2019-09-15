using JuMP, Gurobi, Ipopt;

# forward pass of the SDDP algorithm
function noDisruptionBuild(Δt, T, fData, bData, dData, pDistr, cutDict, solveOpt = true)
    # precalculate data
    T = dData.T;
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    # construct the first stage without disruption occurring
    mp = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "ma27"));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in 1:T]);
    @variable(mp, q[k in fData.brList, t in 1:T]);
    @variable(mp, fData.Vmin[i] <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]);
    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lq[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, u[i in bData.IDList] >= 0);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lp[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lq[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, lineThermal[k in fData.brList, t in 1:T], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    @constraint(mp, powerflow[k in fData.brList, t in 1:T], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bThermal[i in bData.IDList, t in 1:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvmax[i in bData.IDList, t in 1:T], w[i,t] <= bData.cap[i]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # set up the cuts, here tp is the disruption time
    for tp in 2:T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*(sp[i] - cutDict[tp,ω][l].sphat[i]) for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*(w[i] - cutDict[tp,ω][l].what[i]) +
                            cutDict[tp,ω][l].μ[i]*(u[i] - cutDict[tp,ω][l].uhat[i]) for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList));
    for tp in 1:maximum(keys(pDistr.tDistrn))
        dExpr = @expression(mp, 0);
        if tp < T
            for t in 1:tp
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lp[i,t] + lq[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*(dExpr + sum(pDistr.ωDistrn[ω]*θ[tp + 1,ω] for ω in Ω));
        else
            for t in 1:T
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lp[i,t] + lq[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*dExpr;
        end
    end
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        optimize!(mp);
        mpObj = objective_value(mp);
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        for i in fData.genIDList
            for t in 1:T
                solSp[i,t] = value(sp[i,t]);
                solSq[i,t] = value(sq[i,t]);
            end
        end
        for i in bData.IDList
            solu[i] = value(u[i]);
            for t in 1:T
                solw[i,t] = value(w[i,t]);
            end
        end
        for i in fData.IDList
            for t in 1:T
                solLp[i,t] = value(lp[i,t]);
                solLq[i,t] = value(lq[i,t]);
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq);
        return sol,mpObj;
    else
        return mp;
    end
end

function fBuild(td, ωd, currentSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict, solveOpt = true)
    # precalculate data
    T = dData.T;
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    mp = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "ma27"));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList,t in (td - 1):T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i] <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]);
    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lq[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, u[i in bData.IDList] >= 0);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω]);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lp[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lq[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, lineThermal[k in fData.brList, t in td:T; ((k[1],k[2]) != ωd)&((k[2],k[1]) != ωd)], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    @constraint(mp, lineDown[k in fData.brList, t in td:T; ((k[1],k[2]) == ωd)|((k[2],k[1]) == ωd)], p[k,t]^2 + q[k,t]^2 <= 0);
    @constraint(mp, powerflow[k in fData.brList, t in td:T; ((k[1],k[2]) != ωd)&((k[2],k[1]) != ωd)], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; i != ωd], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; i != ωd], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bThermal[i in bData.IDList, t in td:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvmax[i in bData.IDList, t in td:T], w[i,t] <= bData.cap[i]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the cuts, here tp is the disruption time
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*(sp[i] - cutDict[tp,ω][l].sphat[i]) for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*(w[i] - cutDict[tp,ω][l].what[i]) +
                            cutDict[tp,ω][l].μ[i]*(u[i] - cutDict[tp,ω][l].uhat[i]) for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    objExpr = @expression(mp, 0);
    for tp in 1:maximum(keys(pDistr.tDistrn))
        if tp <= T - (td + τ)
            dExpr = @expression(mp, 0);
            for t in td:(tp + td + τ - 1)
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lp[i,t] + lq[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*(dExpr + sum(pDistr.ωDistrn[ω]*θ[tp + td + τ,ω] for ω in Ω));
        else
            dExpr = @expression(mp, 0);
            for t in td:T
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lp[i,t] + lq[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*dExpr;
        end
    end
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        optimize!(mp);
        mpObj = objective_value(mp);
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        for i in fData.genIDList
            for t in td:T
                solSp[i,t] = value(sp[i,t]);
                solSq[i,t] = value(sq[i,t]);
            end
        end
        for i in bData.IDList
            solu[i] = value(u[i]);
            for t in td:T
                solw[i,t] = value(w[i,t]);
            end
        end
        for i in fData.IDList
            for t in td:T
                solLp[i,t] = value(lp[i,t]);
                solLq[i,t] = value(lq[i,t]);
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq);
        return sol,mpObj;
    else
        return mp;
    end
end

function constructForwardM(td, ωd, sol, τ, Δt, T, fData, bData, dData, pDistr, cutDict)
    # construct the math program given the state variables and current stage
    if td == 1
        # if it is the no-disruption problem
        sol,objV = noDisruptionBuild(Δt, T, fData, bData, dData, pDistr, cutDict);
    else
        # if it is f_{ht}^ω
        sol,objV = fBuild(td, ωd, sol, τ, Δt, T, fData, bData, dData, pDistr, cutDict);
    end
    return sol,objV;
end

function exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict)
    # execution of forward pass
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    objV = 0;
    currentLB = 0;
    for n in 1:N
        # for each trial path
        disT = 1;
        ωd = 0;
        costn = 0;
        solHist = [];
        currentSol = solData(Dict(),Dict(),Dict(),Dict(),Dict(),Dict());
        while disT <= T
            # solve the current stage problem, state variables are passed
            nowT = disT;
            currentSol,objV = constructForwardM(disT, ωd, currentSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict);
            push!(solHist,(currentSol,nowT,ωd));

            # generate disruption
            tp,ωd = genScenario(pDistr);
            if nowT == 1
                currentLB = objV;
                disT += tp;
                # calculate the cost of the solution until the next disruption time
                costn += sum(sum(fData.cz*(currentSol.lp[i,t] + currentSol.lq[i,t]) for i in fData.IDList) for t in nowT:(disT - 1)) +
                    sum(currentSol.u[i]*bData.cost[i] for i in bData.IDList);
                for t in nowT:(disT - 1)
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            costn += fData.cp[i].params[1]*(currentSol.sp[i,t]^2) + fData.cp[i].params[2]*currentSol.sp[i,t];
                        elseif fData.cp[i].n == 2
                            costn += fData.cp[i].params[1]*currentSol.sp[i,t];
                        end
                    end
                end
            else
                disT += tp + τ;
                disT = min(disT, T + 1);
                # calculate the cost of the solution until the next disruption time
                costn += sum(sum(fData.cz*(currentSol.lp[i,t] + currentSol.lq[i,t]) for i in fData.IDList) for t in nowT:(disT - 1));
                for t in nowT:(disT - 1)
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            costn += fData.cp[i].params[1]*(currentSol.sp[i,t]^2) + fData.cp[i].params[2]*currentSol.sp[i,t];
                        elseif fData.cp[i].n == 2
                            costn += fData.cp[i].params[1]*currentSol.sp[i,t];
                        end
                    end
                end
            end
        end
        solDict[n] = solHist;
        costDict[n] = costn;
    end
    return solDict, currentLB, costDict;
end
