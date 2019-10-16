# running scripts
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
pmap(i -> importIpopt(),1:30);

fileAdd = "case13_ieee.m";
fData = readStatic(fileAdd,500);
disAdd = "testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "testDataP_96_Fixed.csv";
qAdd = "testDataQ_96_Fixed.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "testDataB.csv";
bData = readBattery(bAdd,"csv");

τ = 4;
T = 16;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,1/4,T);

τ = 4;
T = 24;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,4/T,T);
