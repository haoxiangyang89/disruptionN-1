# test to compare allGen and only disruption gen
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
T = 24;
τ = Int64(1/6*T);
Δt = 0.25;

for ci in 1:length(caseList)
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, 30, false,false, 2, 2);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(200/N)),20), max(Int64(round(200/N)),20));
    elapsedT = time() - startT;
    dataList["dOnly"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("dOnlyResults_$(ci).jld","data",dataList);


    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, max(Int64(round(200/N)),20), max(Int64(round(200/N)),20));
    elapsedT = time() - startT;
    dataList["allGen"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("dOnlyResults_$(ci).jld","data",dataList);
end
