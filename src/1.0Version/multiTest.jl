# test multiple components broken
# upper bound tests
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
#pmap(i -> importIpopt(),1:30);

dataList = Dict();
iterMax = 20;

caseList = [13,33,123];
T = 96;
Δt = 0.25;
N = 5;
pathTrain = load("pathHist_600.jld");
ci = 3;

for j in procs()
    remotecall_fetch(readInData_tau,j,ci,caseList,T,1e4,0,"Multi");
end

cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, 30, false,false, 2, 2);

startPGT = time();
cutDictPG = preGen(T, Δt, N, iterMax);
preGenT = time() - startPGT;
cutDictPGOri = deepcopy(cutDictPG);

startT = time();
# cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
#     Int64(round(100/N)),Int64(round(100/N)), cutDictPG, false, 200, [], 0, pathDict);
cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
    max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, true, 500, [], 0, pathDict);
elapsedT = time() - startT;
dataList = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT,preGenT];

save("MultiResults.jld","MultiOut",dataList);
