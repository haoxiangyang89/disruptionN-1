# main program structure construction

function solveMain(τ, T, fData, pDistr, bData, dData, N, printLog = false, iterMin = 1000, ϵ = 1e-2)
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
        trialPaths,currentLB,currentUBList = exeForward(τ, T, fData, pDistr, bData, dData, N, cutDict);
        currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
        currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
        if (currentLB <= currentUBu)&(currentLB >= currentUBl)&(iterNo >= iterMin)
            keepIter = false;
        else
            cutDict = exeBackward(τ, T, fData, pDistr, bData, dData, trialPaths, cutDict);
        end
    end
end
