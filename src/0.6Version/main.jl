# main program structure construction

function solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, ubGen = false, iterMin = 100, iterMax = 1000, cutDict = Dict(), ubM = 200, hardened = [])
    # readin data and execute the SDDP algorithm

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    keepIter = true;
    LBHist = [];
    UBHist = [];
    UBuHist = [];
    UBlHist = [];
    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        # forward pass: obtain the trial paths
        if !(ubGen)
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict, Dict(), hardened);
            push!(LBHist,currentLB);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
            push!(UBHist,mean(currentUBList));
            currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
            push!(UBuHist,currentUBu);
            currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
            push!(UBlHist,currentUBl);
        else
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict, Dict(), hardened);
            push!(LBHist,currentLB);
            trialPathsUB,currentLBUB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, ubM, cutDict, Dict(), hardened);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
            push!(UBHist,mean(currentUBList));
            currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
            push!(UBuHist,currentUBu);
            currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
            push!(UBlHist,currentUBl);
        end
        if (currentLB <= currentUBu)&(currentLB >= currentUBl)&(iterNo >= iterMin)
            keepIter = false;
        else
            cutDict = exeBackwardAll(τ, T, Δt, fData, pDistr, bData, dData, trialPaths, cutDict, hardened);
        end
        println("========= Iteration $(iterNo) Finished, LB = $(round(currentLB,2)), UB = [$(round(currentUBl,2)),$(round(currentUBu,2))] =========")
    end
    return cutDict,LBHist,UBHist,UBuHist,UBlHist;
end

function solveMain_simuRules(τ, T, Δt, fData, pDistr, bData, dData, N, simuRule, ubGen = false, iterMin = 100, iterMax = 1000, cutDict = Dict(), ubM = 200, hardened = [])
    # readin data and execute the SDDP algorithm
    # simulation rule 0: MC simulation
    # simulation rule 1: to cover as many time period as it can

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    keepIter = true;
    LBHist = [];
    UBHist = [];
    UBuHist = [];
    UBlHist = [];
    pathHist = [];
    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        if simuRule == 0
            pathDict = Dict();
        elseif simuRule == 1
            pathDict = pathSimu_cover(pathHist, N);
        end
        # forward pass: obtain the trial paths
        if !(ubGen)
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict, pathDict, hardened);
            push!(LBHist,currentLB);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
            push!(UBHist,mean(currentUBList));
            currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
            push!(UBuHist,currentUBu);
            currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
            push!(UBlHist,currentUBl);
        else
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, N, cutDict, pathDict, hardened);
            push!(LBHist,currentLB);
            trialPathsUB,currentLBUB,currentUBDict = exeForward(τ, T, Δt, fData, bData, dData, pDistr, ubM, cutDict, Dict(), hardened);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
            push!(UBHist,mean(currentUBList));
            currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
            push!(UBuHist,currentUBu);
            currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
            push!(UBlHist,currentUBl);
        end
        if (currentLB <= currentUBu)&(currentLB >= currentUBl)&(iterNo >= iterMin)
            keepIter = false;
        else
            cutDict = exeBackward(τ, T, Δt, fData, pDistr, bData, dData, trialPaths, cutDict, hardened);
        end
        println("========= Iteration $(iterNo) Finished, LB = $(round(currentLB,2)), UB = [$(round(currentUBl,2)),$(round(currentUBu,2))] =========")
    end
    return cutDict,LBHist,UBHist,UBuHist,UBlHist;
end
