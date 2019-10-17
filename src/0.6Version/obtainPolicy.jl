# obtain the policies for different time horizon
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

Δt = 0.25;
N = 30;
NN = 1000;
fileAdd = "case13_ieee.m";
fData = readStatic(fileAdd,500);
disAdd = "testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "testDataP_96_Fixed.csv";
qAdd = "testDataQ_96_Fixed.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "testDataB.csv";
bData = readBattery(bAdd,"csv");

pathDict = Dict();
TList = [12,24,48,72,96];
for T in TList
    τ = Int64(1/6*T);
    pDistr = modifyT(pDistr,4/T,T);

    pathListData = pmap(i -> simuPath(τ,T,pDistr), 1:NN);
    pathDict[T] = Dict();
    for i in 1:length(pathListData)
        pathDict[T][i] = pathListData[i];
    end

end
save("pathHist.jld","pathDict",pathDict);

policyDict = Dict();
for T in TList
    τ = Int64(1/6*T);
    pDistr = modifyT(pDistr,4/T,T);
    tic();
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, 20, 20);
    runTime = toc();
    policyDict[T] = [cutDict,LBHist,UBHist,UBuHist,UBlHist,runTime];
    save("policy_$(T).jld","policyDict",policyDict);
end
