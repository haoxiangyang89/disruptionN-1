# test 2: deterministic vs. stochastic
include("run.jl");

NN = 1000;
solDet,costDet = exeDet(τ, T, Δt, fData, bData, dData, pDistr, NN);
listDet = [costDet[i] for i in 1:NN];
meanDet = mean(listDet);
sigmaDet = std(listDet);
println(round(meanDet,digits = 2),round(meanDet - 1.96*sigmaDet,digits = 2),round(meanDet + 1.96*sigmaDet,digits = 2));

cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, fData, bData, dData, pDistr, NN, cutDict);
listSDDP = [costSDDP[i] for i in 1:NN];
meanSDDP = mean(listSDDP);
sigmaSDDP = std(listSDDP);
println(round(meanSDDP,digits = 2),round(meanSDDP - 1.96*sigmaSDDP,digits = 2),round(meanSDDP + 1.96*sigmaSDDP,digits = 2));

outputData = [solDet,costDet,meanDet,sigmaDet,
    cutDict,LBHist,UBHist,UBuHist,UBlHist,solSDDP,LBSDDP,costSDDP,meanSDDP,sigmaSDDP];
save("detOut.jld","data",outputData);
