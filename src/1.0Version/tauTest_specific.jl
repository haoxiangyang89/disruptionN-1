# test 5 lines with scenario specific tau
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

dataList = Dict();

caseList = [13,33,123];
T = 24;
Δt = 0.25;
pathTrain = load("pathHist_600.jld");
iterMax = 20;
NN = 1000;

τList = [2,4,6];
for ci in 1:length(caseList)
    dataList[ci] = Dict();
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
        #remotecall_fetch(readInData_old,j,T,ωSet0,10000,0);
    end
end
