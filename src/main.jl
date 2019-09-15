# main program structure construction

function solveMain(τ, T, fData, pDistr, bData, dData, N, printLog = false, iterMin = 100, iterMax = 1000)
    # readin data and execute the SDDP algorithm

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    cutDict = Dict();
    keepIter = true;
    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        # forward pass: obtain the trial paths
        trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict);
        currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
        currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
        currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
        if (currentLB <= currentUBu)&(currentLB >= currentUBl)&(iterNo >= iterMin)
            keepIter = false;
        else
            cutDict = exeBackward(τ, T, Δt, fData, pDistr, bData, dData, trialPaths, cutDict);
        end
    end
end
