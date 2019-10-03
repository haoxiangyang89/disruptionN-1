# test 3: number of trial paths tests
NList = [5,10,20,50,75,100,150,200];
dataList = Dict();

for N in NList
    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
    elapsedT = time() - startT;
    dataList[N] = [cutDict,LBHist,UBHist,UBuHist,UBlHist,elapsedT];
end
