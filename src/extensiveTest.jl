# test 1: compare SDDP and extensive formulation
include("run_small.jl");

mp = Model(with_optimizer(Ipopt.Optimizer, linear_solver = "ma27"));
startTE = time();
mp = extForm(mp, 1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr);
optimize!(mp);
elapsedTE = time() - startTE;
mpObj = objective_value(mp);

startT = time();
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 100);
elapsedT = time() - startT;

outputData = [elapsedTE,mpObj,cutDict,LBHist,UBHist,UBuHist,UBlHist];
save("extOut.jld","data",outputData);
