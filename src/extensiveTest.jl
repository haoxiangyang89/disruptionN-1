mp = Model(with_optimizer(Ipopt.Optimizer, linear_solver = "ma27"));
startTE = time();
mp = extForm(mp, 1, 0, [], 1, τ, Δt, T, fData, bData, dData, pDistr);
optimize!(mp);
elapsedTE = time() - startTE;
mpObj = objective_value(mp);

N = 10;
startT = time();
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
elapsedT = time() - startT;
