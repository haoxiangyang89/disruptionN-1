# a diagnostic small test for extensive formulation
τ = 1;
T = 3;
Δt = 0.25;
N = 10;
pDistr = modifyT(pDistr,1/3,T);
ΩSet = [1];
for ω in keys(pDistr.ωDistrn)
    if !(ω in ΩSet)
        pDistr = modifyOmega(pDistr,ω);
    end
end

mExt = Model(solver = IpoptSolver(linear_solver = "ma27"));
mExt = extForm(mExt, 1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr);
solve(mExt);
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20);
