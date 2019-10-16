# test 3 hardening tests
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
T = 24;
Δt = 0.25;
N = 30;

# obtain the no hardening costs for comparison
ΩSet = [3,5,(2,9),(8,12),(10,13)];
pDistr = readDisruption(disAdd,"csv");
pDistr = modifyT(pDistr,1/2,T);

data = [];
for ω in keys(pDistr.ωDistrn)
    if !(ω in ωSet[i])
        pDistr = modifyOmega(pDistr,ω);
    end
end
cutDict0,LBHist0,UBHist0,UBuHist0,UBlHist0 = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 25, 25, cutDict0);
push!(data,(cutDict0,LBHist0,UBHist0,UBuHist0,UBlHist0));
for ω in ΩSet
    pDistrNew = modifyOmega(pDistr,ω);
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistrNew, bData, dData, N, false, 20, 20);
    push!(data,(cutDict,LBHist,UBHist,UBuHist,UBlHist));
end

save("hardOut.jld","data",outputData);
