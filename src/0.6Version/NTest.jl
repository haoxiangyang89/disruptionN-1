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
pDistr = modifyT(pDistr,4/T,T);
Δt = 0.25;

for ci in 1:length(caseList)
    fileAdd = "case$(caseList[ci])_ieee.m";
    fData = readStatic(fileAdd,10000);
    disAdd = "testProbRead_$(caseList[ci]).csv";
    pDistr = readDisruption(disAdd,"csv");
    pAdd = "testDataP_$(caseList[ci]).csv";
    qAdd = "testDataQ_$(caseList[ci]).csv";
    dData = readDemand(pAdd,qAdd,"csv");
    bAdd = "testDataB_$(caseList[ci]).csv";
    bData = readBattery(bAdd,"csv");
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, 30, false, 2, 2);

    for N in NList
        startT = time();
        cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 50, 50);
        elapsedT = time() - startT;
        dataList[N] = [cutDict,LBHist,UBHist,UBuHist,UBlHist,elapsedT];

        save("NResults_$(ci).jld","NOut",dataList);
    end
end
