# test 4a: effectiveness of different τ
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

fileAdd = "case13_ieee.m";
fData = readStatic(fileAdd,500);
disAdd = "testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "testDataP_96_Fixed.csv";
qAdd = "testDataQ_96_Fixed.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "testDataB.csv";
bData = readBattery(bAdd,"csv");

T = 24;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,4/T,T);

τList = [2,4,6,8,10];
data = Dict();
for τ in τList
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20);
    pathListData = pmap(i -> simuPath(τ,T,pDistr), 1:NN);
    pathDict = Dict();
    for i in 1:NN
        pathDict[i] = pathListData[i];
    end
    solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, fData, bData, dData, pDistr, NN, cutDict, pathDict);
    listSDDP = [costSDDP[i] for i in 1:NN];
    meanSDDP = mean(listSDDP);
    sigmaSDDP = std(listSDDP);
    data[τ] = [LBHist,LBSDDP,listSDDP,meanSDDP,sigmaSDDP];
end

save("tauResults.jld","data",data);
