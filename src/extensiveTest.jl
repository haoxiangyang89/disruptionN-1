# test 1: compare SDDP and extensive formulation
include("run.jl");

mp = Model();
startTE = time();
mp = extForm(mp, 1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr);
# optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
#     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma97", acceptable_tol = 1e-8));

elapsedTE = time() - startTE;
mpObj = objective_value(mp);

startT = time();
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 100);
elapsedT = time() - startT;

outputData = [elapsedTE,mpObj,cutDict,LBHist,UBHist,UBuHist,UBlHist];
save("extOut.jld","data",outputData);
