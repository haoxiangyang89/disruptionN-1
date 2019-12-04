# test the utilization rate of batteries
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;
NN = 1000;
τ = 4;
pathListDRaw = load("pathHist_1000.jld");

for ci in 1:length(caseList)
    TList = [24,36,48,72,96];

    pathDictA = pathListDRaw["pathDict"][ci];

    for T in TList
        dataList[ci] = Dict();
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T);
        end

        # select a preset pathDict
        pathDict = pathDictA[T];
        solDet,costDet = exeDet(τ, T, Δt, fData, bData, dData, pDistr, NN, pathDict);

        # train the stochastic strategy
        cutDictPG = preGen(τ, T, Δt, N, iterMax);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(200/N)),20), max(Int64(round(200/N)),20),cutDictPG);
        for j in procs()
            remotecall_fetch(cutIni,j,cutDict);
        end

        solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, NN, false, pathDict);

        # read the battery utilization data from solDet and solSDDP
        uDet = solDet[1][1][1].u;
        uSDDP = solSDDP[1][1][1].u;
        zDet = Dict();
        zSDDP = Dict();
        wDet = Dict();
        wSDDP = Dict();
        for i in 1:NN
            zDet[i] = Dict();
            zSDDP[i] = Dict();
            wDet[i] = Dict();
            wSDDP[i] = Dict();
            for b in bData.IDList
                zDet[i][b] = Dict();
                zSDDP[i][b] = Dict();
                wDet[i][b] = Dict();
                wSDDP[i][b] = Dict();
            end

            # for each scenario, obtain the deterministic solution and the stochastic solution
            dDetList = [solDet[i][j][2] for j in 1:length(solDet[i])];
            push!(dDetList,T+1);
            for j in 1:length(solDet[i])
                # obtain the inventory & actual injection to network until the next solution
                for t in dDetList[j]:(dDetList[j+1] - 1)
                    for b in bData.IDList
                        zDet[i][b][t] = solDet[i][j][1].zp;
                        wDet[i][b][t] = solDet[i][j][1].w;
                    end
                end
            end

            dSDDPList = [solSDDP[i][j][2] for j in 1:length(solSDDP[i])];
            push!(dSDDPList,T+1);
            for j in 1:length(solSDDP[i])
                # obtain the inventory & actual injection to network until the next solution
                for t in dSDDPList[j]:(dSDDPList[j+1] - 1)
                    for b in bData.IDList
                        zSDDP[i][b][t] = solSDDP[i][j][1].zp;
                        wSDDP[i][b][t] = solSDDP[i][j][1].w;
                    end
                end
            end
        end
        save("butilResults.jld","detOut",[uDet,zDet,wDet],"stochOut",[uSDDP,zSDDP,wSDDP]);
    end
end
