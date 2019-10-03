# backward pass of the SDDP algorithm
function fBuild_D(td, ωd, currentPath, τ, Δt, T, fData, bData, dData, pDistr, cutDict, solveOpt = true)
    prevtpInd = maximum([i for i in 1:length(currentPath) if currentPath[i][2] < td]);
    currentSol = currentPath[prevtpInd][1];
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

    #mp = Model(solver = GurobiSolver(GUROBI_ENV, OutputFlag = 0, QCPDual = 1));
    mp = Model(solver = IpoptSolver(print_level = 0, linear_solver = "ma27", constr_viol_tol = 1e-10));

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
    @constraint(mp, bThermal1[i in bData.IDList, t in td:T], zp[i,t] <= u[i]);
    @constraint(mp, bThermal2[i in bData.IDList, t in td:T], zq[i,t] <= u[i]);
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
                        sum(cutDict[tp,ω][l].λ[i]*(sp[i,tp - 1] - cutDict[tp,ω][l].sphat[i]) for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*(w[i,tp - 1] - cutDict[tp,ω][l].what[i]) +
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
        #optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", print_level = 0, acceptable_tol = 1e-8, max_iter = 10000));
        statusMp = solve(mp);
        println(statusMp, " ", td, " ", ωd);
        # obtain the primal solutions & obj value
        vhat = getobjectivevalue(mp);
        # obtain the solutions
        solSp = Dict();
        solw = Dict();
        solu = Dict();
        for i in fData.genIDList
            solSp[i] = getvalue(sp[i,td - 1]);
        end
        for i in bData.IDList
            solu[i] = getvalue(u[i]);
            solw[i] = getvalue(w[i,td - 1]);
        end
        # obtain the dual solutions
        dsolλ = Dict();
        dsolγ = Dict();
        dsolμ = Dict();
        for i in fData.genIDList
            dsolλ[i] = getdual(spIni[i]);
            if abs(dsolλ[i]) < 1e-4
                if dsolλ[i] > 0
                    vhat += dsolλ[i]*(-solSp[i]);
                elseif dsolλ[i] < 0
                    vhat += dsolλ[i]*(fData.Pmax[i] - solSp[i]);
                end
                dsolλ[i] = 0;
            end
        end
        for i in bData.IDList
            dsolγ[i] = getdual(bInvIni[i]);
            if abs(dsolγ[i]) < 1e-4
                if dsolγ[i] > 0
                    vhat += dsolγ[i]*(-solw[i]);
                elseif dsolγ[i] < 0
                    vhat += dsolγ[i]*(bData.cap[i] - solw[i]);
                end
            end
            dsolμ[i] = getdual(uIni[i]);
            if abs(dsolμ[i]) < 1e-4
                if dsolμ[i] > 0
                    vhat += dsolμ[i]*(-solu[i]);
                elseif dsolμ[i] < 0
                    vhat += dsolμ[i]*(bData.uCap[i] - solu[i]);
                end
            end
        end

        cutTemp = cutData(statusMp,dsolλ,dsolγ,dsolμ,vhat,solSp,solw,solu);
        return cutTemp;
    else
        return mp;
    end
end

function constructBackwardM(td, τ, T, Δt, fData, pDistr, bData, dData, trialPaths, matchedTrial, cutDict)
    # construct the math program given the state variables and current stage
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    paraSet = Iterators.product(Ω,matchedTrial);

    cutCurrentData = pmap(item -> fBuild_D(td, item[1], trialPaths[item[2]], τ, Δt, T, fData, bData, dData, pDistr, cutDict), paraSet);
    # cutCurrentData is a list
    for ω in Ω
        itemInd = 0;
        for item in paraSet
            itemInd += 1;
            if item[1] == ω
                if (cutCurrentData[itemInd].solStatus == :Optimal)
                    if (td,ω) in keys(cutDict)
                        push!(cutDict[td,ω],cutCurrentData[itemInd]);
                    else
                        cutDict[td,ω] = [cutCurrentData[itemInd]];
                    end
                end
            end
        end
    end

    # for ω in Ω
    #     # solve the later stage problem
    #     cutCurrent = fBuild_D(td, ω, prevSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict);
    #     if (td,ω) in keys(cutDict)
    #         push!(cutDict[td,ω],cutCurrent);
    #     else
    #         cutDict[td,ω] = [cutCurrent];
    #     end
    # end
    return cutDict;
end

function exeBackward(τ, T, Δt, fData, pDistr, bData, dData, trialPaths, cutDict)
    # execution of forward pass
    # input: trialPaths: the collection of trial points
    #        cutDict: previously generated cuts
    # output: update the cutDict
    tpDict = Dict();
    for n in keys(trialPaths)
        tpDict[n] = [trialPaths[n][i][2] for i in 1:length(trialPaths[n])];
    end
    for t in T:-1:2
        matchedTrial = [];
        for n in keys(trialPaths)
            if t in tpDict[n]
                push!(matchedTrial,n);
            end
        end
        if trialPaths != []
            cutDict = constructBackwardM(t, τ, T, Δt, fData, pDistr, bData, dData, trialPaths, matchedTrial, cutDict);
        end
        println("Time $(t) Passed");
    end
    return cutDict;
end
