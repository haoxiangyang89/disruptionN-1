# deterministic optimization models
function detBuild(Δt, T, fData, bData, dData, solveOpt = true)
    # deterministic formulation
    T = dData.T;
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end

    # construct the first stage without disruption occurring
    mp = Model();

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
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]^2);
    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lpp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, u[i in bData.IDList] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t]  +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, pequal[k in fData.brList, t in 1:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in 1:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    @constraint(mp, lineThermal[k in fData.brList, t in 1:T], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    @constraint(mp, powerflow[k in fData.brList, t in 1:T], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bThermal[i in bData.IDList, t in 1:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvmax[i in bData.IDList, t in 1:T], w[i,t] <= bData.cap[i]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # deterministic objective function
    objExpr = calObjDet(mp, T, fData, bData, sp, lpp, lqp, lpm, lqm, u);

    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma97", acceptable_tol = 1e-8, print_level = 0));

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
                if (value(lpp[i,t]) > 1e-8)|(value(lpm[i,t]) > 1e-8)
                    solLp[i,t] = value(lpp[i,t]) - value(lpm[i,t]);
                else
                    solLp[i,t] = 0;
                end
                if (value(lqp[i,t]) > 1e-8)|(value(lqm[i,t]) > 1e-8)
                    solLq[i,t] = value(lqp[i,t]) - value(lqm[i,t]);
                else
                    solLq[i,t] = 0;
                end
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq);
        return sol,mpObj;
    else
        return mp;
    end
end

function fDetBuild(td, ωd, currentSol, τ, Δt, T, fData, bData, dData, solveOpt = true)
    # precalculate data
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Bparams = Dict();
    for t in td:T
        # create B parameters
        for k in fData.brList
            # if the line is disrupted and it is within disruption time
            if (((k[1],k[2]) == ωd)|((k[2],k[1]) == ωd))&(t <= td + τ)
                Bparams[k,t] = 0;
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i == ωd)&(t <= td + τ)
                Bparams[i,t] = 0;
            else
                Bparams[i,t] = 1;
            end
        end
    end
    for i in fData.genIDList
        Bparams[i,td - 1] = 1;
    end

    mp = Model();

    # set up the variables
    @variable(mp, fData.Pmin[i]*Bparams[i,t] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]*Bparams[i,t]);
    @variable(mp, fData.Qmin[i]*Bparams[i,t] <= sq[i in fData.genIDList,t in td:T] <= fData.Qmax[i]*Bparams[i,t]);
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
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]^2);
    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lpp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, u[i in bData.IDList] >= 0);

    # set up the constraints
    bigM = 1000;
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, pequal[k in fData.brList, t in td:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in td:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    @constraint(mp, lineThermal2[k in fData.brList, t in td:T;Bparams[k,t] == 0], p[k,t] == 0);
    @constraint(mp, lineThermal3[k in fData.brList, t in td:T;Bparams[k,t] == 0], q[k,t] == 0);
    @constraint(mp, powerflow1[k in fData.brList, t in td:T;Bparams[k,t] == 1], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bThermal[i in bData.IDList, t in td:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvmax[i in bData.IDList, t in td:T], w[i,t] <= bData.cap[i]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the objective function
    objExpr = @expression(mp, fData.cz*sum(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList) for t in td:T));
    for t in td:T
        for i in fData.genIDList
            # add generator cost
            if fData.cp[i].n == 3
                objExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
            elseif fData.cp[i].n == 2
                objExpr += fData.cp[i].params[1]*sp[i,t];
            end
        end
    end
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma97", acceptable_tol = 1e-8, print_level = 0));

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
                if (value(lpp[i,t]) > 1e-8)|(value(lpm[i,t]) > 1e-8)
                    solLp[i,t] = value(lpp[i,t]) - value(lpm[i,t]);
                else
                    solLp[i,t] = 0;
                end
                if (value(lqp[i,t]) > 1e-8)|(value(lqm[i,t]) > 1e-8)
                    solLq[i,t] = value(lqp[i,t]) - value(lqm[i,t]);
                else
                    solLq[i,t] = 0;
                end
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq);
        return sol,mpObj;
    else
        return mp;
    end
end

function constructDetM(td, ωd, sol, τ, Δt, T, fData, bData, dData)
    # construct the math program given the state variables and current stage
    if td == 1
        # if it is the no-disruption problem
        sol,objV = detBuild(Δt, T, fData, bData, dData);
    else
        # if it is f_{ht}^ω
        sol,objV = fDetBuild(td, ωd, sol, τ, Δt, T, fData, bData, dData);
    end
    return sol,objV;
end

function buildPathDet(τ, T, Δt, fData, bData, dData, pDistr)
    disT = 1;
    ωd = 0;
    costn = 0;
    solHist = [];
    currentSol = solData(Dict(),Dict(),Dict(),Dict(),Dict(),Dict());
    while disT <= T
        # solve the current stage problem, state variables are passed
        nowT = disT;
        currentSol,objV = constructDetM(disT, ωd, currentSol, τ, Δt, T, fData, bData, dData);
        push!(solHist,(currentSol,nowT,ωd));

        # generate disruption
        tp,ωd = genScenario(pDistr);
        if nowT == 1
            disT += tp;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn += sum(sum(fData.cz*(abs(currentSol.lp[i,t]) + abs(currentSol.lq[i,t])) for i in fData.IDList) for t in nowT:(disT - 1)) +
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
            costn += sum(sum(fData.cz*(abs(currentSol.lp[i,t]) + abs(currentSol.lq[i,t])) for i in fData.IDList) for t in nowT:(disT - 1));
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
    return [solHist,costn];
end

function exeDet(τ, T, Δt, fData, bData, dData, pDistr, N)
    # execution of forward pass
    # input: N: the number of trial points;
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    returnData = pmap(i -> buildPathDet(τ, T, Δt, fData, bData, dData, pDistr), 1:N);
    for n in 1:N
        solDict[n] = returnData[n][1];
        costDict[n] = returnData[n][2];
    end
    return solDict, costDict;
end
