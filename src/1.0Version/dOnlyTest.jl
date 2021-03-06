# test to compare allGen and only disruption gen
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
T = 96;
τ = Int64(1/6*T);
Δt = 0.25;
pathTrain = load("pathHist_600.jld");

for ci in 1:length(caseList)
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T,τ);
    end

    pathDict = pathTrain["pathDict"][ci][T];
    pathDict = reverseScen(pathDict,τ,pDistr);

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, 30, false,false, 2, 2);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), Dict(), false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["dOnly"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("dOnlyResults_$(ci).jld","data",dataList);
end
