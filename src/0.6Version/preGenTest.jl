# test to check if preGen helps reduce the solution time
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
T = 96;
τ = Int64(1/6*T);
Δt = 0.25;
iterMax = 20;
pathTrain = load("pathHist_600.jld");

for ci in 1:length(caseList)
    pathDict = pathTrain["pathDict"][ci][T];
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, 30, false,false, 2, 2);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), Dict(),
        false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["dOnly"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("pgResults_$(ci).jld","data",dataList);

    startPGT = time();
    cutDictPG = preGen(τ, T, Δt, N, iterMax);
    preGenT = time() - startPGT;

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG,
        false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["preGen"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT,preGenT];

    save("pgResults_$(ci).jld","data",dataList);
end
