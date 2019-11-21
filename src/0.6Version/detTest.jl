# test 2: deterministic vs. stochastic
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

caseList = [13,33,123];

for ci in 1:length(caseList)
    fileAdd = "case$(caseList[ci])_ieee.m";
    fData = readStatic(fileAdd,10000);
    disAdd = "testProbRead_$(caseList[ci]).csv"
    pDistr = readDisruption(disAdd,"csv");
    pAdd = "testDataP_$(caseList[ci]).csv";
    qAdd = "testDataQ_$(caseList[ci]).csv";
    dData = readDemand(pAdd,qAdd,"csv");
    bAdd = "testDataB_$(caseList[ci]).csv";
    bData = readBattery(bAdd,"csv");

    NN = 1000;
    Δt = 0.25;
    TList = [12,24,48,72,96];
    pathListDRaw = load("pathHist.jld");
    pathDictA = pathListDRaw["pathDict"];
    stochDataRaw = load("policy_trained.jld");
    detOut = Dict();
    stochOut = Dict();
    for T in TList
        τ = Int64(1/6*T);
        pDistr = modifyT(pDistr,4/T,T);

        # select a preset pathDict
        pathDict = pathDictA[T];
        solDet,costDet = exeDet(τ, T, Δt, fData, bData, dData, pDistr, NN, pathDict);
        listDet = [costDet[i] for i in 1:NN];
        meanDet = mean(listDet);
        sigmaDet = std(listDet);
        println(round(meanDet,2)," ",round(meanDet - 1.96*sigmaDet,2)," ",round(meanDet + 1.96*sigmaDet,2));
        detOut[T] = [costDet,listDet,meanDet,sigmaDet];

        cutDict = stochDataRaw["policyDict"][T][1];
        solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, fData, bData, dData, pDistr, NN, cutDict, pathDict);
        listSDDP = [costSDDP[i] for i in 1:NN];
        meanSDDP = mean(listSDDP);
        sigmaSDDP = std(listSDDP);
        println(round(meanSDDP,2)," ",round(meanSDDP - 1.96*sigmaSDDP,2)," ",round(meanSDDP + 1.96*sigmaSDDP,2));
        stochOut[T] = [LBSDDP, costSDDP, listSDDP, meanSDDP, sigmaSDDP];

        save("detResults_$(ci).jld","detOut",detOut,"stochOut",stochOut);
    end
end

# for T in TList
#     println(T," & ",round(detOut[T][3],2)," & ",round(stochOut[T][4],2)," & ",round((detOut[T][3] - stochOut[T][4])/detOut[T][3]*100,2)," & ");
# end

# pathListData = pmap(i -> simuPath(τ,T,pDistr), 1:NN);
# pathDict = Dict();
# for i in 1:length(pathListData)
#     pathDict[i] = pathListData[i];
# end
# solDet,costDet = exeDet(τ, T, Δt, fData, bData, dData, pDistr, NN, pathDict);
# listDet = [costDet[i] for i in 1:NN];
# meanDet = mean(listDet);
# sigmaDet = std(listDet);
# println(round(meanDet,2)," ",round(meanDet - 1.96*sigmaDet,2)," ",round(meanDet + 1.96*sigmaDet,2));
#
# cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, 20, 20);
# solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, fData, bData, dData, pDistr, NN, cutDict, pathDict);
# listSDDP = [costSDDP[i] for i in 1:NN];
# meanSDDP = mean(listSDDP);
# sigmaSDDP = std(listSDDP);
# println(round(meanSDDP,2)," ",round(meanSDDP - 1.96*sigmaSDDP,2)," ",round(meanSDDP + 1.96*sigmaSDDP,2));
#
# outputData = [solDet,costDet,meanDet,sigmaDet,
#     cutDict,LBHist,UBHist,UBuHist,UBlHist,solSDDP,LBSDDP,costSDDP,meanSDDP,sigmaSDDP];
# save("detOut.jld","data",outputData);
