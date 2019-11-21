# running scripts
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
pmap(i -> importIpopt(),1:30);

caseList = [13,33,123];

i = 3;
fileAdd = "case$(caseList[i])_ieee.m";
fData = readStatic(fileAdd,10000);
disAdd = "testProbRead_$(caseList[i]).csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "testDataP_$(caseList[i]).csv";
qAdd = "testDataQ_$(caseList[i]).csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "testDataB_$(caseList[i]).csv";
bData = readBattery(bAdd,"csv");

τ = 4;
T = 12;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,1/4,T);

τ = 4;
T = 24;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,4/T,T);
