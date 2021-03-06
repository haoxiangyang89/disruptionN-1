# running scripts for a smaller example
using Distributed;
addprocs(30);
@everywhere include("./src/loadMod.jl");

fileAdd = "./test/case13_ieee.m";
fData = readStatic(fileAdd);
disAdd = "./test/testProbRead.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "./test/testDataP.csv";
qAdd = "./test/testDataQ.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "./test/testDataB.csv";
bData = readBattery(bAdd,"csv");

τ = 4;
T = 16;
Δt = 0.25;
N = 50;
