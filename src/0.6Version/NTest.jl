# test 3: number of trial paths tests
# test 2: deterministic vs. stochastic
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

NList = [1,5,10,20,30];
dataList = Dict();

caseList = [13,33,123];
T = 96;
τ = Int64(1/6*T);
Δt = 0.25;

for ci in 1:length(caseList)
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, 3, false, 2, 2);

    for N in NList
        startT = time();
        cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, N, false, max(Int64(round(500/i)),20), max(Int64(round(500/i)),20));
        elapsedT = time() - startT;
        dataList[N] = [cutDict,LBHist,UBHist,UBuHist,UBlHist,elapsedT];

        save("NResults_$(ci).jld","NOut",dataList);
    end
end
