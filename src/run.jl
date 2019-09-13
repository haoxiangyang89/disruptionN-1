# running scripts
include("loadMod.jl");

fileAdd = "./test/case13_ieee.m";
fData = readStatic(fileAdd);
disAdd = "./test/testProbRead.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "./test/testDataP.csv";
qAdd = "./test/testDataQ.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "./test/testDataB.csv";
bData = readBattery(bAdd,"csv");
