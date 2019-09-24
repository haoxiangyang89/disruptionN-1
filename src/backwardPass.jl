# backward pass of the SDDP algorithm
function fBuild_D(td, ωd, currentSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict, solveOpt = true)
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

    mp = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "ma27"));

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

    bigM = 1000;
    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]^2);
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
    @constraint(mp, pequal[k in fData.brList, t in td:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in td:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    @constraint(mp, lineThermal[k in fData.brList, t in td:T], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2*Bparams[k,t]);
    @constraint(mp, powerflow1[k in fData.brList, t in td:T], v[k[2],t] <= v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]) + (1 - Bparams[k,t])*bigM);
    @constraint(mp, powerflow2[k in fData.brList, t in td:T], v[k[2],t] >= v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]) - (1 - Bparams[k,t])*bigM);
    @constraint(mp, rampUp[i in fData.genIDList, t in td:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i] + bigM*(1 - Bparams[i,t]));
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i] - bigM*(1 - Bparams[i,t]));
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
                        sum(cutDict[tp,ω][l].λ[i]*(sp[i,tp - 1] - cutDict[tp,ω][l].sphat[i]) for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*(w[i,tp - 1] - cutDict[tp,ω][l].what[i]) +
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
        # obtain the dual solutions
        dsolλ = Dict();
        dsolγ = Dict();
        dsolμ = Dict();
        for i in fData.genIDList
            dsolλ[i] = dual(spIni[i]);
        end
        for i in bData.IDList
            dsolγ[i] = dual(bInvIni[i]);
            dsolμ[i] = dual(uIni[i]);
        end

        # obtain the primal solutions & obj value
        vhat = objective_value(mp);
        # obtain the solutions
        solSp = Dict();
        solw = Dict();
        solu = Dict();
        for i in fData.genIDList
            solSp[i] = value(sp[i,td - 1]);
        end
        for i in bData.IDList
            solu[i] = value(u[i]);
            solw[i] = value(w[i,td - 1]);
        end

        cutTemp = cutData(dsolλ,dsolγ,dsolμ,vhat,solSp,solw,solu);
        return cutTemp;
    else
        return mp;
    end
end

function constructBackwardM(td, τ, T, Δt, fData, pDistr, bData, dData, prevSol, cutDict)
    # construct the math program given the state variables and current stage
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    cutCurrentData = pmap(ω -> fBuild_D(td, ω, prevSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict), Ω);
    for ωInd in 1:length(Ω)
        ω = Ω[ωInd];
        if (td,ω) in keys(cutDict)
            push!(cutDict[td,ω],cutCurrentData[ωInd]);
        else
            cutDict[td,ω] = [cutCurrentData[ωInd]];
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
        for n in keys(trialPaths)
            currentPath = trialPaths[n];
            if t in tpDict[n]
                # if the disruption time is in the trial path, generate cuts
                # using the solution from the previous disruption
                prevtpInd = maximum([i for i in 1:length(currentPath) if currentPath[i][2] < t]);
                prevSol = currentPath[prevtpInd][1];
                cutDict = constructBackwardM(t, τ, T, Δt, fData, pDistr, bData, dData, prevSol, cutDict);
                println("Time $(t) Trial $(n) Passed");
            end
        end
    end
    return cutDict;
end
