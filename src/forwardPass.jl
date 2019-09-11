using JuMP, Gurobi, Ipopt;

# forward pass of the SDDP algorithm
function noDisruptionBuild(T, fData, bData, dData, cutDict)
    # construct the first stage without disruption occurring
    mp = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "ma27"));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end

    sqhatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sqhatsum[i] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in 1:T]);
    @variable(mp, q[k in fData.brList, t in 1:T]);
    @variable(mp, fData.Vmin[i] <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]);
    @variable(mp, 0 <= w[i in bData.IDList, t in 1:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lq[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, u[i in bData.IDList] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) +
        sphatsum[i,t] - )

    # set up the cuts

    # set up the objective function

    # solve the problem

    # obtain the solutions

    return sol;
end

function fBuild(td, ωd, sol, τ, T, fData, bData, dData, cutDict)
    mp = Model(with_optimizer(Ipopt.Optimizer, print_level = 0, linear_solver = "ma27"));
    return sol;
end

function constructForwardM(td, ωd, sol, τ, T, fData, bData, dData, N, cutDict)
    # construct the math program given the state variables and current stage
    if td == 0
        # if it is the no-disruption problem
        sol = noDisruptionBuild(T, fData, pDistr, bData, dData, cutDict);
    else
        # if it is f_{kt}^ω
        sol = fBuild(td, ωd, sol, τ, T, fData, bData, dData, cutDict);
    end
end

function exeForward(τ, T, fData, pDistr, bData, dData, N, cutDict)
    # execution of forward pass
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts
    for n in 1:N
        # for each trial path
        disT = 0;
        currentω = 0;
        solPack =
        while disT > T
            # solve the current stage problem, state variables are passed
            solPack = constructForwardM(disT,currentω,solPack);

            # generate disruption
            t,ω = genScenario(pDistr);
            disT += t;
        end
    end
end
