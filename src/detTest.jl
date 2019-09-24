# test 2: deterministic vs. stochastic
startTD = time();
solOut,objV = detBuild(Δt, T, fData, bData, dData);
elapsedTD = time() - startTD;

N = 10;
startT = time();
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
elapsedT = time() - startT;
