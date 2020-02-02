# test 1: compare SDDP and extensive formulation
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
pmap(i -> importIpopt(),1:30);

τ = 2;
T = 12;
Δt = 0.25;
N = 30;
###########################################################################
data = [];
for T in [8,12,16]
    for j in procs()
        remotecall_fetch(readInData_old,j,T);
    end

    # global mExt = Model(solver = IpoptSolver(linear_solver = "ma27"));
    global mExt = Model(solver = GurobiSolver());
    startTE = time();
    global mExt = extForm(1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr, true);
    solve(mExt);
    elapsedTE = time() - startTE;
    mExtObj = getobjectivevalue(mExt);

    startT = time();
    #cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, 20, 20);
    elapsedT = time() - startT;
    push!(data,(T,elapsedTE,mExtObj,elapsedT,LBHist));
end
save("extensiveSDDP.jld","data",data);

τ = 2;
T = 8;
Δt = 0.25;
N = 30;
ωSet = [[5,(2,9),(8,12)],
        [3,5,(2,9),(8,12)],
        [3,5,(2,9),(8,12),(10,13)]];

for i in 1:3
    pDistr = readDisruption(disAdd,"csv");
    pDistr = modifyT(pDistr,1/2,T);

    for ω in keys(pDistr.ωDistrn)
        if !(ω in ωSet[i])
            pDistr = modifyOmega(pDistr,ω);
        end
    end

    pDistr = modifyT(pDistr,1/2,T);
    mp = Model(solver = IpoptSolver(linear_solver = "ma27"));
    startTE = time();
    mp = extForm(mp, 1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr);
    solve(mp);
    elapsedTE = time() - startTE;
    mpObj = getobjectivevalue(mp);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20);
    elapsedT = time() - startT;
    push!(data,(T,elapsedT,cutDict,LBHist,UBHist,UBuHist,UBlHist));
end
