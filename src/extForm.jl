function extForm(mp, td, ωd, inheritData, baseProb, τ, Δt, T, fData, bData, dData, pDistr)
    # extensive formulation could not have variable/constraint names
    # inheritData: [1]: sp, [2]: w, [3]: u (only contains the information for the linking time period)

    println("========= Disruption time $(td), scenario $(ωd) modeling =========");
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

    # set up the variables
    spDict = Dict();
    sqDict = Dict();
    for i in fData.genIDList
        for t in (td - 1):T
            spDict[i,t] = @variable(mp, lower_bound = fData.Pmin[i],  upper_bound = fData.Pmax[i]);
            sqDict[i,t] = @variable(mp, lower_bound = fData.Qmin[i],  upper_bound = fData.Qmax[i]);
        end
    end
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += spDict[j,t];
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
                    sqhatsum[i,t] += sqDict[j,t];
                end
            end
        end
    end

    pDict = Dict();
    qDict = Dict();
    vDict = Dict();
    wDict = Dict();
    yDict = Dict();
    zpDict = Dict();
    zqDict = Dict();
    lpDict = Dict();
    lqDict = Dict();
    uDict = Dict();
    for k in fData.brList
        for t in td:T
            pDict[k,t] = @variable(mp);
            qDict[k,t] = @variable(mp);
        end
    end
    for i in fData.IDList
        for t in td:T
            vDict[i,t] = @variable(mp, lower_bound = fData.Vmin[i]^2, upper_bound = fData.Vmax[i]^2);
            lpDict[i,t] = @variable(mp, lower_bound = 0);
            lqDict[i,t] = @variable(mp, lower_bound = 0);
        end
    end
    for i in bData.IDList
        wDict[i,td - 1] = @variable(mp, lower_bound = 0, upper_bound = bData.cap[i]);
        uDict[i] = @variable(mp, lower_bound = 0);
        for t in td:T
            wDict[i,t] = @variable(mp, lower_bound = 0, upper_bound = bData.cap[i]);
            yDict[i,t] = @variable(mp);
            zpDict[i,t] = @variable(mp);
            zqDict[i,t] = @variable(mp);
        end
    end

    # set up the constraints
    bigM = 1000;
    for i in fData.IDList
        for t in td:T
            @constraint(mp, sum(zpDict[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpDict[i,t] +
                sphatsum[i,t] - dData.pd[i][t] == sum(pDict[k,t] for k in fData.branchDict1[i]));
            @constraint(mp, sum(zqDict[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqDict[i,t] +
                sqhatsum[i,t] - dData.qd[i][t] == sum(qDict[k,t] for k in fData.branchDict1[i]));
        end
    end
    for i in fData.genIDList
        if td != 1
            @constraint(mp, spDict[i,td - 1] == inheritData[1][i]);
        end
        for t in td:T
            if t != 1
                @constraint(mp, spDict[i,t] - spDict[i,t - 1] <= fData.RU[i] + bigM*(1 - Bparams[i,t]));
                @constraint(mp, spDict[i,t] - spDict[i,t - 1] >= fData.RD[i] - bigM*(1 - Bparams[i,t]));
            end
        end
    end
    for k in fData.brList
        for t in td:T
            @constraint(mp, pDict[k,t] == -pDict[(k[2],k[1],k[3]),t]);
            @constraint(mp, qDict[k,t] == -qDict[(k[2],k[1],k[3]),t]);
            @constraint(mp, pDict[k,t]^2 + qDict[k,t]^2 <= fData.rateA[k]^2*Bparams[k,t]);
            @constraint(mp, vDict[k[2],t] <= vDict[k[1],t] - 2*(Rdict[k]*pDict[k,t] + Xdict[k]*qDict[k,t]) + (1 - Bparams[k,t])*bigM);
            @constraint(mp, vDict[k[2],t] >= vDict[k[1],t] - 2*(Rdict[k]*pDict[k,t] + Xdict[k]*qDict[k,t]) - (1 - Bparams[k,t])*bigM);
        end
    end
    for i in bData.IDList
        if td != 1
            @constraint(mp, uDict[i] == inheritData[3][i]);
        end
        @constraint(mp, wDict[i,td - 1] == inheritData[2][i]);
        for t in td:T
            @constraint(mp, wDict[i,t] == wDict[i,t - 1] - yDict[i,t]*Δt);
            @constraint(mp, zpDict[i,t]^2 + zqDict[i,t]^2 <= uDict[i]^2);
            for l in 1:length(bData.ηα[i])
                @constraint(mp, zpDict[i,t] <= bData.ηα[i][l]*yDict[i,t] + bData.ηβ[i][l]);
            end
            @constraint(mp, wDict[i,t] <= bData.cap[i]);
        end
    end

    # recursion through the possible scenario
    tList = sort([t for t in keys(pDistr.tDistrn)]);
    objExpr = objective_function(mp);
    for tp in tList
        if td == 1
            objExpr += sum(bData.cost[i]*uDict[i] for i in bData.IDList);
            if td + tp > T
                # if the next disruption is over the time horizon
                dExpr = fData.cz*sum(sum(lpDict[i,t] + lqDict[i,t] for i in fData.IDList) for t in td:T);
                for t in td:T
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                objExpr += baseProb*pDistr.tDistrn[tp]*dExpr;
            else
                # if not
                dExpr = fData.cz*sum(sum(lpDict[i,t] + lqDict[i,t] for i in fData.IDList) for t in td:(td + tp - 1));
                for t in td:(td + tp - 1)
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                for ω in Ω
                    objExpr += baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω]*dExpr;
                    spInherit = Dict();
                    wInherit = Dict();
                    uInherit = Dict();
                    for i in fData.genIDList
                        spInherit[i] = spDict[i,td + tp - 1];
                    end
                    for i in bData.IDList
                        wInherit[i] = wDict[i,td + tp - 1];
                        uInherit[i] = uDict[i];
                    end
                    inheritData = [spInherit,wInherit,uInherit];
                    mp = extForm(mp, td + tp, ω, inheritData, baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω], τ, Δt, T, fData, bData, dData, pDistr);
                end
            end
        else
            if td + τ + tp > T
                # if the next disruption is over the time horizon
                # add to the objective function
                dExpr = fData.cz*sum(sum(lpDict[i,t] + lqDict[i,t] for i in fData.IDList) for t in td:T);
                for t in td:T
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                objExpr += baseProb*pDistr.tDistrn[tp]*dExpr;
            else
                # if not
                # add the current part to the objective function
                # recursion to the next disruption
                dExpr = fData.cz*sum(sum(lpDict[i,t] + lqDict[i,t] for i in fData.IDList) for t in td:(td + tp + τ - 1));
                for t in td:(td + tp + τ - 1)
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                for ω in Ω
                    objExpr += baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω]*dExpr;
                    spInherit = Dict();
                    wInherit = Dict();
                    uInherit = Dict();
                    for i in fData.genIDList
                        spInherit[i] = spDict[i,td + tp + τ - 1];
                    end
                    for i in bData.IDList
                        wInherit[i] = wDict[i,td + tp + τ - 1];
                        uInherit[i] = uDict[i];
                    end
                    inheritData = [spInherit,wInherit,uInherit];
                    mp = extForm(mp, td + tp + τ, ω, inheritData, baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω], τ, Δt, T, fData, bData, dData, pDistr);
                end
            end
        end
    end
    return mp;
end
