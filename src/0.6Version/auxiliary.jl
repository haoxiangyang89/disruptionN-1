# auxiliary functions

using Distributions;
import Base.rand;

struct CategoricalSamplerNew <: Sampleable{Univariate,Discrete}
    mass::Vector{Float64}
    category::Vector{Float64}
end

function rand(s::CategoricalSamplerNew)
    s.category[rand(Categorical(s.mass))]
end

function genScenario(pDistr)
    # generate a disruption time
    tSupport = [i for i in keys(pDistr.tDistrn)];
    tProb = [pDistr.tDistrn[i] for i in tSupport];
    tDistrObj = Categorical(tProb);
    t = rand(tDistrObj);

    # generate a disruption location
    ωSupport = [i for i in keys(pDistr.ωDistrn)];
    ωProb = [pDistr.ωDistrn[i] for i in ωSupport];
    ωDistrObj = Categorical(ωProb);
    ω = rand(ωDistrObj);
    return tSupport[t],ωSupport[ω];
end

function modifyOmega(pDistr,hardComp)
    ωDistrNew = Dict();
    releaseProb = pDistr.ωDistrn[hardComp];
    avgNo = length(values(pDistr.ωDistrn)) - 1;
    for i in keys(pDistr.ωDistrn)
        if i != hardComp
            # if it is not the hardened component
            ωDistrNew[i] = pDistr.ωDistrn[i] + releaseProb/avgNo;
        end
    end
    pDistrNew = probDistrn(pDistr.tDistrn,ωDistrNew);
    return pDistrNew;
end

function modifyT(pDistr,λD,T)
    # modify the time distribution using the given λD
    tDistrnNew = Dict();
    for t in 1:(T-1)
        tDistrnNew[t] = exp(-λD*(t - 1)) - exp(-λD*t);
    end
    # probability of the time T+1
    tDistrnNew[T] = 1 - sum([tDistrnNew[t] for t in 1:(T-1)]);
    pDistrNew = probDistrn(tDistrnNew,pDistr.ωDistrn);
    return pDistrNew;
end

function calObj(mp, td, τ, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ)
    # function to formulate the objective function for f_D
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    objExpr = @expression(mp, 0);
    for tp in 1:maximum(keys(pDistr.tDistrn))
        dExpr = @expression(mp, 0);
        if tp <= T - (td + τ)
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
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*(dExpr + sum(pDistr.ωDistrn[ω]*θ[tp + td + τ,ω] for ω in Ω));
        else
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
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*dExpr;
        end
    end
    return objExpr;
end

function calObjDet(mp, T, fData, bData, sp, lpp, lqp, lpm, lqm, u)
    # function to formulate the objective function for deterministic case

    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        fData.cz*sum(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList) for t in 1:T));
    for t in 1:T
        for i in fData.genIDList
            # add generator cost
            if fData.cp[i].n == 3
                objExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
            elseif fData.cp[i].n == 2
                objExpr += fData.cp[i].params[1]*sp[i,t];
            end
        end
    end
    return objExpr;
end

function calObj1(mp, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ, u)
    # function to formulate the objective function for f_1
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

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
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
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
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*dExpr;
        end
    end
    return objExpr;
end

function calDualC1(T, fData, pDistr)
    # function to formulate the objective function for f_D
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    fDict = Dict();
    lDict = Dict();
    θDict = Dict();
    for tp in 1:maximum(keys(pDistr.tDistrn))
        if tp < T
            for t in 1:tp
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
            for ω in Ω
                θDict[tp + 1,ω] = pDistr.tDistrn[tp]*pDistr.ωDistrn[ω];
            end
        else
            for t in 1:T
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
        end
    end
    return fDict,lDict,θDict;
end

function calCostF(costn, currentSol, T, fData, nowT, disT)
    # function to calculate costs for forward pass
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
    return costn;
end

function calDualC(td, τ, T, fData, pDistr)
    # function to formulate the objective function for f_D
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    fDict = Dict();
    lDict = Dict();
    θDict = Dict();
    for tp in 1:maximum(keys(pDistr.tDistrn))
        if tp <= T - (td + τ)
            for t in td:(tp + td + τ - 1)
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
            for ω in Ω
                θDict[tp + td + τ,ω] = pDistr.tDistrn[tp]*pDistr.ωDistrn[ω];
            end
        else
            for t in td:T
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
        end
    end
    return fDict,lDict,θDict;
end

function simuPath(τ,T,pDistr)
    pathList = [];
    nowT = 1;
    while nowT <= T
        tp,ωd = genScenario(pDistr);
        push!(pathList, (tp,ωd));
        if nowT == 1
            nowT += tp;
            nowT = min(nowT, T + 1);
        else
            nowT += tp + τ;
            nowT = min(nowT, T + 1);
        end
    end
    return pathList;
end

function pathSimu_cover(pathHist, N)
    # maximize the coverage from the sample
end
