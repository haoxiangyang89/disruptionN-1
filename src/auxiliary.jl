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
