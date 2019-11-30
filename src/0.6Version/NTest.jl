# test 3: number of trial paths tests
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
#pmap(i -> importIpopt(),1:30);

NList = [1,5,10,20,30];
dataList = Dict();
iterMax = 20;
NN = 20;

caseList = [13,33,123];
T = 96;
τ = Int64(1/6*T);
Δt = 0.25;

for ci in 1:length(caseList)
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, 30, false,false, 2, 2);

    for N in NList
        startPGT = time();
        cutDictPG = preGen(τ, T, Δt, N, iterMax);
        preGenT = time() - startPGT;

        startT = time();
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(200/N)),20), max(Int64(round(200/N)),20), cutDictPG);
        elapsedT = time() - startT;
        dataList[N] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT,preGenT];

        save("NResults_$(ci).jld","NOut",dataList);
    end
end
