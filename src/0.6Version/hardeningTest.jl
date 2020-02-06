# test 3 hardening tests
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;
pathTrain = load("pathHist_600.jld");

for ci in 1:length(caseList)
    pathDict = pathTrain["pathDict"][ci][T];
    dataList = Dict();
    τ = Int64(1/6*T);
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T,1e4);
        #remotecall_fetch(readInData_old,j,T,[(2,9),(8,12),(10,13)],10000,0);
    end
    cutDictPG = preGen(τ, T, Δt, N, iterMax, false, Dict(), []);
    # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
    #     20, 20, cutDictPG, false, 0, [], 0, pathDict);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, false, 0, []);

    dataList["NoD"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist];

    for ω in keys(pDistr.ωDistrn)
    #for ω in [(2,9),(8,12),(10,13)]
        cutDictPG = preGen(τ, T, Δt, N, iterMax, false, Dict(), [ω]);
        # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
        #     20, 20, cutDictPG, false, 0, [ω], 0, pathDict);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false,
            max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, false, 0, [ω]);
        dataList[ω] = [LBHist,UBHist,UBuHist,UBlHist,timeHist];
    end

    save("hardResults_$(ci).jld","data",dataList);
end
