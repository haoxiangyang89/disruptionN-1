# auxiliary functions

using Distributions;
import Base.rand;

function rand(b :: DiscreteNonParametric)
    n = rand();
    i = 0;
    while n >= 0
        i += 1;
        n -= b.p[i];
    end
    x = b.support[i];
    return x;
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
    ωDistrNew = copy(pDistr.ωDistrn);
    releaseProb = ωDistrNew[hardComp];
    avgNo = length(values(pDistr.ωDistrn)) - 1;
    for i in keys(ωDistrNew)
        if i == hardComp
            # if it is the hardened component
            ωDistrNew[i] = 0;
        else
            ωDistrNew[i] += releaseProb/avgNo;
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
