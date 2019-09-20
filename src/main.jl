# main program structure construction

function solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, printLog = false, iterMin = 100, iterMax = 1000)
    # readin data and execute the SDDP algorithm

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    cutDict = Dict();
    keepIter = true;
    LBHist = [];
    UBHist = [];
    UBuHist = [];
    UBlHist = [];
    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        # forward pass: obtain the trial paths
        trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict);
        push!(LBHist,currentLB);
        currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
        push!(UBHist,mean(currentUBList));
        currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
        push!(UBuHist,currentUBu);
        currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
        push!(UBlHist,currentUBl);
        if (currentLB <= currentUBu)&(currentLB >= currentUBl)&(iterNo >= iterMin)
            keepIter = false;
        else
            cutDict = exeBackward(τ, T, Δt, fData, pDistr, bData, dData, trialPaths, cutDict);
        end
        println("========= Iteration $(iterNo) Finished, LB = $(round(currentLB,2)), UB = $(round(mean(currentUBList),2)) +/- $(round(1.96*std(currentUBList),2)) =========")
    end
    return cutDict,LBHist,UBHist,UBuHist,UBlHist;
end
