# test 2: deterministic vs. stochastic
include("run.jl");

NN = 1000;
pathListData = pmap(i -> simuPath(τ,T,pDistr), 1:NN);
pathDict = Dict();
for i in 1:length(pathListData)
    pathDict[i] = pathListData[i];
end
solDet,costDet = exeDet(τ, T, Δt, fData, bData, dData, pDistr, NN, pathDict);
listDet = [costDet[i] for i in 1:NN];
meanDet = mean(listDet);
sigmaDet = std(listDet);
println(round(meanDet,2)," ",round(meanDet - 1.96*sigmaDet,2)," ",round(meanDet + 1.96*sigmaDet,2));

cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, 20, 20);
solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, fData, bData, dData, pDistr, NN, cutDict, pathDict);
listSDDP = [costSDDP[i] for i in 1:NN];
meanSDDP = mean(listSDDP);
sigmaSDDP = std(listSDDP);
println(round(meanSDDP,2),round(meanSDDP - 1.96*sigmaSDDP,2),round(meanSDDP + 1.96*sigmaSDDP,2));

outputData = [solDet,costDet,meanDet,sigmaDet,
    cutDict,LBHist,UBHist,UBuHist,UBlHist,solSDDP,LBSDDP,costSDDP,meanSDDP,sigmaSDDP];
save("detOut.jld","data",outputData);
