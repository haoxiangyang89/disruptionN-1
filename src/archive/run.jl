# running scripts
using Distributed;
addprocs(30);
@everywhere include("./src/loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
pmap(i -> importIpopt(),1:30);

fileAdd = "./test/case13_ieee_linCost.m";
fData = readStatic(fileAdd,10000);
disAdd = "./test/testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "./test/testDataP_96_Fixed.csv";
qAdd = "./test/testDataQ_96_Fixed.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "./test/testDataB.csv";
bData = readBattery(bAdd,"csv");

τ = 4;
T = 16;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,1/4,T);
