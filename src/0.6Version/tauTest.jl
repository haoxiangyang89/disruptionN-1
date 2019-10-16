# test 4a: effectiveness of different τ
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

fileAdd = "case13_ieee.m";
fData = readStatic(fileAdd,500);
disAdd = "testProbRead_96.csv"
pDistr = readDisruption(disAdd,"csv");
pAdd = "testDataP_96_Fixed.csv";
qAdd = "testDataQ_96_Fixed.csv";
dData = readDemand(pAdd,qAdd,"csv");
bAdd = "testDataB.csv";
bData = readBattery(bAdd,"csv");

T = 24;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,4/T,T);

τList = [2,4,6,8,10];
data = [];
for τ in τList
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, 20, 20);
    push!(data,(cutDict,LBHist,UBHist,UBuHist,UBlHist));
end

cutDictT = Dict();
LBHistT = Dict();
UBHistT = Dict();
UBuHistT = Dict();
UBlHistT = Dict();
for τ in 0:Int64(round(T/2))
    cutDictT[τ],LBHistT[τ],UBHistT[τ],UBuHistT[τ],UBlHistT[τ] =
        solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
    println("Disruption Length $(τ): LB = $(LBHistT[τ][length(LBHistT[τ])])");
end

outputData = [cutDictT,LBHistT,UBHistT,UBuHistT,UBlHistT];
save("tauOut.jld","data",outputData);
