# test 3 hardening tests
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;

for ci in 1:length(caseList)
    dataList = Dict();
    τ = Int64(1/6*T);
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
    end

    for ω in keys(pDistr.ωDistrn)
        cutDictPG = preGen(τ, T, Δt, N, iterMax, false, Dict(), [ω]);

        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false,
            max(Int64(round(300/N)),20), max(Int64(round(300/N)),20), cutDictPG, false, 0, [ω]);
        dataList[ω] = [LBHist,UBHist,UBuHist,UBlHist,timeHist];
    end

    save("hardResults_$(ci).jld","data",dataList);
end
