# test 3 hardening tests
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

fileAdd = "case13_ieee.m";
fData = readStatic(fileAdd,10000);
disAdd = "testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "testDataP_96_Fixed.csv";
qAdd = "testDataQ_96_Fixed.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "testDataB.csv";
bData = readBattery(bAdd,"csv");

Δt = 0.25;
N = 30;
T = 12;
τ = Int64(1/6*T);
pDistr = modifyT(pDistr,4/T,T);

# obtain the no hardening costs for comparison
ΩSet = [3,5,(2,9),(8,12),(10,13)];
for ω in keys(pDistr.ωDistrn)
    if !(ω in ΩSet)
        pDistr = modifyOmega(pDistr,ω);
    end
end
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20);
data = Dict();
data[0] = [cutDict,LBHist,UBHist,UBuHist,UBlHist];

for ω in ΩSet
    pDistrNew = modifyOmega(pDistr,ω);
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistrNew, bData, dData, N, false, 20, 20);
    # cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20, Dict(), 0, [ω]);
    data[ω] = [cutDict,LBHist,UBHist,UBuHist,UBlHist];
    save("hardOut_10000.jld","data",data);
end
