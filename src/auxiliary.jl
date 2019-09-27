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
