# test the relationship between the cost coefficients and the run time performance
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
T = 24;
τ = Int64(1/6*T);
Δt = 0.25;
pathTrain = load("pathHist_600.jld");

cList = [0.01,0.05,0.1,0.2,0.3,0.5];
ci = 2;
pathDict = pathTrain["pathDict"][ci][T];


for cmulti in cList
    # adjust the cost
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
        remotecall_fetch(changeCost,j,fData,bData,cmulti);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, 30, false,false, 2, 2);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(300/N)),20), max(Int64(round(300/N)),20), Dict(),
        false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["dOnly"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("cchangeResults_$(ci).jld","data",dataList);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 20, 20, Dict(),
        false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["allGen"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("cchangeResults_$(ci).jld","data",dataList);

end
