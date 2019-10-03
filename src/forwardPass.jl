using JuMP, Gurobi, Ipopt;

# forward pass of the SDDP algorithm
function noDisruptionBuild(Δt, T, fData, bData, dData, pDistr, cutDict, solveOpt = true)
    # precalculate data
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

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
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
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
    objExpr = calObj1(mp, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ, u);
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", acceptable_tol = 1e-8, print_level = 0));
        mpObj = objective_value(mp);
        println("First stage, solving status $(termination_status(mp))");
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        for i in fData.genIDList
            for t in 1:T
                if value(sp[i,t]) > 1e-6
                    solSp[i,t] = value(sp[i,t]);
                else
                    solSp[i,t] = 0;
                end
                if value(sq[i,t]) > 1e-6
                    solSq[i,t] = value(sq[i,t]);
                else
                    solSq[i,t] = 0;
                end
            end
        end
        for i in bData.IDList
            solu[i] = value(u[i]);
            for t in 1:T
                if value(w[i,t]) > 1e-6
                    solw[i,t] = value(w[i,t]);
                else
                    solw[i,t] = 0;
                end
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

function fBuild(td, ωd, currentSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict, solveOpt = true)
    # precalculate data
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
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
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω] >= 0);

    # set up the constraints
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
    objExpr = calObj(mp, td, τ, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ);
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", acceptable_tol = 1e-8, print_level = 0));
        mpObj = objective_value(mp);
        println("Disruption time $(td), scenario $(ωd), solving status $(termination_status(mp))");
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        for i in fData.genIDList
            for t in td:T
                if value(sp[i,t]) > 1e-6
                    solSp[i,t] = value(sp[i,t]);
                else
                    solSp[i,t] = 0;
                end
                if value(sq[i,t]) > 1e-6
                    solSq[i,t] = value(sq[i,t]);
                else
                    solSq[i,t] = 0;
                end
            end
        end
        for i in bData.IDList
            solu[i] = value(u[i]);
            for t in td:T
                if value(w[i,t]) > 1e-6
                    solw[i,t] = value(w[i,t]);
                else
                    solw[i,t] = 0;
                end
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

function buildPath(τ, T, Δt, fData, bData, dData, pDistr, cutDict, pathList = [])
    disT = 1;
    ωd = 0;
    costn = 0;
    solHist = [];
    currentLB = 0;
    currentSol = solData(Dict(),Dict(),Dict(),Dict(),Dict(),Dict());
    iter = 1;
    while disT <= T
        # solve the current stage problem, state variables are passed
        nowT = disT;
        currentSol,objV = constructForwardM(disT, ωd, currentSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict);
        push!(solHist,(currentSol,nowT,ωd));

        # generate disruption
        if pathList == []
            tp,ωd = genScenario(pDistr);
        else
            tp,ωd = pathList[iter];
        end
        iter += 1;
        if nowT == 1
            currentLB = objV;
            disT += tp;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn += sum(currentSol.u[i]*bData.cost[i] for i in bData.IDList);
            costn = calCostF(costn, currentSol, T, fData, nowT, disT);
        else
            disT += tp + τ;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn = calCostF(costn, currentSol, T, fData, nowT, disT);
        end
    end
    println("Path Built!");
    return [solHist,currentLB,costn];
end

function exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict, pathDict = Dict())
    # execution of forward pass
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    objV = 0;
    currentLB = 0;
    # for n in 1:N
    #     # for each trial path
    #     returnData = buildPath(τ, T, Δt, fData, bData, dData, pDistr, cutDict);
    #     solDict[n] = returnData[1];
    #     costDict[n] = returnData[2];
    # end
    if pathDict == Dict()
        for i in 1:N
            pathDict[i] = [];
        end
    end
    returnData = pmap(i -> buildPath(τ, T, Δt, fData, bData, dData, pDistr, cutDict, pathDict[i]), 1:N);
    for n in 1:N
        solDict[n] = returnData[n][1];
        costDict[n] = returnData[n][3];
    end
    currentLB = returnData[1][2];
    return solDict, currentLB, costDict;
end
