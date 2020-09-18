# test 5 lines with scenario specific tau
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

caseList = [13,33,123];
T = 24;
Δt = 0.25;
NN = 1000;
N = 5;

dataList = Dict();
for ci in 1:length(caseList)
    for j in procs()
        remotecall_fetch(readInData_tau,j,ci,caseList,T);
        #remotecall_fetch(readInData_old,j,T,ωSet0,10000,0);
    end

    # train the stochastic strategy
    # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 20, 20, cutDictPG);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, true,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20),Dict());

    pathListData = pmap(i -> simuPath(T,pDistr), 1:NN);
    pathDict = Dict();
    for i in 1:NN
        pathDict[i] = pathListData[i];
    end
    solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, true, pathDict);
    listSDDP = [costSDDP[i] for i in 1:NN];
    meanSDDP = mean(listSDDP);
    sigmaSDDP = std(listSDDP);

    dataList[ci] = [LBHist,LBSDDP,listSDDP,meanSDDP,sigmaSDDP];

    save("tauResults_specific.jld","data",dataList);
end

dataListSame = Dict();
for ci in 1:length(caseList)
    for j in procs()
        remotecall_fetch(readInData_tau,j,ci,caseList,T,1e4,0,"Same");
        #remotecall_fetch(readInData_old,j,T,ωSet0,10000,0);
    end

    # train the stochastic strategy
    # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 20, 20, cutDictPG);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20),Dict());

    pathListData = pmap(i -> simuPath(T,pDistr), 1:NN);
    pathDict = Dict();
    for i in 1:NN
        pathDict[i] = pathListData[i];
    end
    solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, false, pathDict);
    listSDDP = [costSDDP[i] for i in 1:NN];
    meanSDDP = mean(listSDDP);
    sigmaSDDP = std(listSDDP);

    dataListSame[ci] = [LBHist,LBSDDP,listSDDP,meanSDDP,sigmaSDDP];

    save("tauResults_specificSame.jld","data",dataListSame);
end
