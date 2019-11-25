# main program structure construction

function solveMain(τ, T, Δt, N, allGen = false, qpopt = false, iterMin = 100, iterMax = 1000, cutDict = Dict(), ubGen = false, ubM = 200, hardened = [])
    # readin data and execute the SDDP algorithm

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    keepIter = true;
    LBHist = [];
    UBHist = [];
    UBuHist = [];
    UBlHist = [];
    timeHist = [];
    # initialize the cutDict
    for j in procs()
        remotecall_fetch(cutIni,j,cutDict);
    end

    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        iterStart = time();
        # forward pass: obtain the trial paths
        if !(ubGen)
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, N, qpopt, Dict(), hardened);
            push!(LBHist,currentLB);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
            push!(UBHist,mean(currentUBList));
            currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
            push!(UBuHist,currentUBu);
            currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
            push!(UBlHist,currentUBl);
        else
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, N, qpopt, Dict(), hardened);
            push!(LBHist,currentLB);
            trialPathsUB,currentLBUB,currentUBDict = exeForward(τ, T, Δt, ubM, qpopt, Dict(), hardened);
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
            if allGen
                exeBackwardAll(τ, T, Δt, trialPaths, qpopt, hardened);
            else
                exeBackward(τ, T, Δt, trialPaths, qpopt, hardened);
            end
        end
        iterElapsed = time() - iterStart;
        push!(timeHist,iterElapsed);
        println("========= Iteration $(iterNo) Finished, LB = $(round(currentLB,2)), UB = [$(round(currentUBl,2)),$(round(currentUBu,2))], Time = $(iterElapsed) sec. =========")
    end
    return cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist;
end

function solveMain_simuRules(τ, T, Δt, N, simuRule, ubGen = false, iterMin = 100, iterMax = 1000, cutDict = Dict(), ubM = 200, hardened = [])
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
    # initialize the cutDict
    for j in procs()
        remotecall_fetch(cutIni,j,cutDict);
    end

    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        if simuRule == 0
            pathDict = Dict();
        elseif simuRule == 1
            pathDict = pathSimu_cover(pathHist, N);
        end
        iterStart = time();
        # forward pass: obtain the trial paths
        if !(ubGen)
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, N, pathDict, hardened);
            push!(LBHist,currentLB);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
            push!(UBHist,mean(currentUBList));
            currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
            push!(UBuHist,currentUBu);
            currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
            push!(UBlHist,currentUBl);
        else
            trialPaths,currentLB,currentUBDict = exeForward(τ, T, Δt, N, pathDict, hardened);
            push!(LBHist,currentLB);
            trialPathsUB,currentLBUB,currentUBDict = exeForward(τ, T, Δt, ubM, Dict(), hardened);
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
            exeBackward(τ, T, Δt, trialPaths, qpopt, hardened);
        end
        iterElapsed = time() - iterStart;
        push!(timeHist,iterElapsed);
        println("========= Iteration $(iterNo) Finished, LB = $(round(currentLB,2)), UB = [$(round(currentUBl,2)),$(round(currentUBu,2))], Time = $(iterElapsed) sec. =========")
    end
    return cutDict,LBHist,UBHist,UBuHist,UBlHist;
end
