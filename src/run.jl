# running scripts
using Distributed;
addprocs(30);
@everywhere include("./src/loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

fileAdd = "./test/case13_ieee.m";
fData = readStatic(fileAdd,10000);
disAdd = "./test/testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "./test/testDataP_96.csv";
qAdd = "./test/testDataQ_96.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "./test/testDataB.csv";
bData = readBattery(bAdd,"csv");

τ = 16;
T = 96;
Δt = 0.25;
N = 50;
