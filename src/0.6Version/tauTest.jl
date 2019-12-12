# test effectiveness of different τ
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;
NN = 1000;

τList = [2,4,6,8,10];

for ci in 1:length(caseList)
    dataList[ci] = Dict();
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
    end

    for τ in τList
        # train the stochastic strategy
        cutDictPG = preGen(τ, T, Δt, N, iterMax);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 10, 10, cutDictPG);
        # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(300/N)),20), max(Int64(round(300/N)),20),cutDictPG);

        pathListData = pmap(i -> simuPath(τ,T,pDistr), 1:NN);
        pathDict = Dict();
        for i in 1:NN
            pathDict[i] = pathListData[i];
        end
        solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, NN, false, pathDict);
        listSDDP = [costSDDP[i] for i in 1:NN];
        meanSDDP = mean(listSDDP);
        sigmaSDDP = std(listSDDP);
        dataList[ci][τ] = [LBHist,LBSDDP,listSDDP,meanSDDP,sigmaSDDP];
        save("tauResults.jld","data",dataList);
    end
end
