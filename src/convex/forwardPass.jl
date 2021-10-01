# forward pass of the SDDP algorithm
# QC relaxation
function build1_QC(Δt, T)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    # construct the first stage without disruption occurring
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0 ,"Threads" => 1));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
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
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    θu = pi/30;

    @variable(mp, p[k in fData.brList, t in 1:T]);
    @variable(mp, q[k in fData.brList, t in 1:T]);
    @variable(mp, fData.Vmin[i] <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]);

    @variable(mp, fData.Vmin[i]^2 <= vhat[i in fData.IDList, t in 1:T] <= fData.Vmax[i]^2);
    @constraint(mp, vhatMock1[i in fData.IDList, t in 1:T], [vhat[i,t] + 1/4; [vhat[i,t] - 1/4, v[i,t]]] in SecondOrderCone());
    @constraint(mp, vhatMock2[i in fData.IDList, t in 1:T], vhat[i,t] - (fData.Vmax[i] + fData.Vmin[i])*v[i,t] <= -fData.Vmax[i]*fData.Vmin[i]);

    @variable(mp, vv[k in fData.brList, t in 1:T]);
    @constraint(mp, vvEquality[k in fData.brList, t in 1:T], vv[k,t] == vv[(k[2],k[1],k[3]),t]);

    @variable(mp, va[i in fData.IDList, t in 1:T]);
    @variable(mp, vad_var[k in fData.brList, t in 1:T]);
    @constraint(mp, vadCon[k in fData.brList, t in 1:T], vad_var[k,t] == va[k[1],t] - va[k[2],t]);
    # set up the reference bus angle = 0
      for i in fData.IDList
          for t in 1:T
              if fData.bType[i] == 3
                  @constraint(mp, va[i,t] == 0);
              end
          end
      end
    @constraint(mp, angleDiff[k in fData.brList, t in 1:T], va[k[1],t] - va[k[2],t] <= θu);

    @variable(mp, cos(θu) <= cs[k in fData.brList, t in 1:T] <= 1);
    @constraint(mp, csMock[k in fData.brList, t in 1:T],
                    [5/4 - cs[k,t]; [sqrt((1-cos(θu))/(θu)^2)*(va[k[1],t] - va[k[2],t]), cs[k,t] - 3/4]] in SecondOrderCone());
    @constraint(mp, csEquality[k in fData.brList, t in 1:T], cs[k,t] == cs[(k[2],k[1],k[3]),t]);

    @variable(mp, -sin(θu) <= ss[k in fData.brList, t in 1:T] <= sin(θu));
    @constraint(mp, ssMock1[k in fData.brList, t in 1:T], ss[k,t] <= cos(θu/2)*(va[k[1],t] - va[k[2],t] - θu/2) + sin(θu/2));
    @constraint(mp, ssMock2[k in fData.brList, t in 1:T], -ss[k,t] <= -(cos(θu/2)*(va[k[1],t] - va[k[2],t] + θu/2) - sin(θu/2)));
    @constraint(mp, ssEquality[k in fData.brList, t in 1:T], ss[k,t] == -ss[(k[2],k[1],k[3]),t]);

    @variable(mp, wc[k in fData.brList, t in 1:T]);
    @variable(mp, ws[k in fData.brList, t in 1:T]);
    @constraint(mp, wcEquality[k in fData.brList, t in 1:T], wc[k,t] == wc[(k[2],k[1],k[3]),t]);
    @constraint(mp, wsEquality[k in fData.brList, t in 1:T], ws[k,t] == -ws[(k[2],k[1],k[3]),t]);

    # McCormick Relaxation
    @constraint(mp, vvMcCormick1[k in fData.brList, t in 1:T],
                vv[k,t] >= fData.Vmin[k[1]]*v[k[2],t] + fData.Vmin[k[2]]*v[k[1],t] - fData.Vmin[k[1]]*fData.Vmin[k[2]]);
    @constraint(mp, vvMcCormick2[k in fData.brList, t in 1:T],
                vv[k,t] >= fData.Vmax[k[1]]*v[k[2],t] + fData.Vmax[k[2]]*v[k[1],t] - fData.Vmax[k[1]]*fData.Vmax[k[2]]);
    @constraint(mp, vvMcCormick3[k in fData.brList, t in 1:T],
                vv[k,t] <= fData.Vmin[k[1]]*v[k[2],t] + fData.Vmax[k[2]]*v[k[1],t] - fData.Vmin[k[1]]*fData.Vmax[k[2]]);
    @constraint(mp, vvMcCormick4[k in fData.brList, t in 1:T],
                vv[k,t] <= fData.Vmax[k[1]]*v[k[2],t] + fData.Vmin[k[2]]*v[k[1],t] - fData.Vmax[k[1]]*fData.Vmin[k[2]]);

    @constraint(mp, wcMcCormick1[k in fData.brList, t in 1:T],
                wc[k,t] >= fData.Vmin[k[1]]*fData.Vmin[k[2]]*cs[k,t] + cos(θu)*vv[k,t] - fData.Vmin[k[1]]*fData.Vmin[k[2]]*cos(θu));
    @constraint(mp, wcMcCormick2[k in fData.brList, t in 1:T],
                wc[k,t] >= fData.Vmax[k[1]]*fData.Vmax[k[2]]*cs[k,t] + 1*vv[k,t] - fData.Vmax[k[1]]*fData.Vmax[k[2]]*1);
    @constraint(mp, wcMcCormick3[k in fData.brList, t in 1:T],
                wc[k,t] <= fData.Vmin[k[1]]*fData.Vmin[k[2]]*cs[k,t] + vv[k,t]*1 - fData.Vmin[k[1]]*fData.Vmin[k[2]]*1);
    @constraint(mp, wcMcCormick4[k in fData.brList, t in 1:T],
                wc[k,t] <= fData.Vmax[k[1]]*fData.Vmax[k[2]]*cs[k,t] + vv[k,t]*cos(θu) - fData.Vmax[k[1]]*fData.Vmax[k[2]]*cos(θu));

    @constraint(mp, wsMcCormick1[k in fData.brList, t in 1:T],
                ws[k,t] >= fData.Vmin[k[1]]*fData.Vmin[k[2]]*ss[k,t] + (-sin(θu))*vv[k,t] - fData.Vmin[k[1]]*fData.Vmin[k[2]]*(-sin(θu)));
    @constraint(mp, wsMcCormick2[k in fData.brList, t in 1:T],
                ws[k,t] >= fData.Vmax[k[1]]*fData.Vmax[k[2]]*ss[k,t] + sin(θu)*vv[k,t] - fData.Vmax[k[1]]*fData.Vmax[k[2]]*sin(θu));
    @constraint(mp, wsMcCormick3[k in fData.brList, t in 1:T],
                ws[k,t] <= fData.Vmin[k[1]]*fData.Vmin[k[2]]*ss[k,t] + vv[k,t]*sin(θu) - fData.Vmin[k[1]]*fData.Vmin[k[2]]*sin(θu));
    @constraint(mp, wsMcCormick4[k in fData.brList, t in 1:T],
                ws[k,t] <= fData.Vmax[k[1]]*fData.Vmax[k[2]]*ss[k,t] + vv[k,t]*(-sin(θu)) - fData.Vmax[k[1]]*fData.Vmax[k[2]]*(-sin(θu)));

    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lpp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, lineThermal[k in fData.brList, t in 1:T], [fData.rateA[k];[p[k,t],q[k,t]]] in SecondOrderCone());
    @constraint(mp, bThermal[i in bData.IDList, t in 1:T], [u[i];[zp[i,t],zq[i,t]]] in SecondOrderCone());
    @constraint(mp, powerflow_p[k in fData.brList, t in 1:T], p[k,t] == fData.g[k]*vhat[k[1],t] - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]);
    @constraint(mp, powerflow_q[k in fData.brList, t in 1:T], q[k,t] == -fData.b[k]*vhat[k[1],t] + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]);

    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # set up the cuts, here tp is the disruption time
    for tp in 2:T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in 1:T]);
    @variable(mp,tAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3] >= 0);
    @variable(mp,tAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);
    @variable(mp,tAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);

    @constraint(mp,gcAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
        (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
        (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
        fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    @constraint(mp, genCost1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3],
        [tAux1[i,t];[tAux2[i,t],tAux3[i,t]]] in SecondOrderCone());
    @constraint(mp, genCost2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC1(T, fData, pDistr);
    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in 1:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in 2:T));

    @objective(mp, Min, objExpr);
    return mp;
end

# SOCP relaxation
function build1_SOCP(Δt, T)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    # construct the first stage without disruption occurring
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0 ,"Threads" => 1));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
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
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in 1:T]);
    @variable(mp, q[k in fData.brList, t in 1:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]^2);

    @variable(mp, vcos[k in fData.brList, t in 1:T] >= 0);
    @variable(mp, vsin[k in fData.brList, t in 1:T]);
    @constraint(mp, vcSymm[k in fData.brList, t in 1:T], vcos[k,t] == vcos[(k[2],k[1],k[3]),t]);
    @constraint(mp, vsSymm[k in fData.brList, t in 1:T], vsin[k,t] == -vsin[(k[2],k[1],k[3]),t]);
    @constraint(mp, socpCon[k in fData.brList, t in 1:T], [1/sqrt(2)*(v[k[1],t] + v[k[2],t]); [1/sqrt(2)*v[k[1],t], 1/sqrt(2)*v[k[2],t],
                                                            vcos[k,t], vsin[k,t]]] in SecondOrderCone());

    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lpp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, lineThermal[k in fData.brList, t in 1:T], [fData.rateA[k];[p[k,t],q[k,t]]] in SecondOrderCone());
    @constraint(mp, bThermal[i in bData.IDList, t in 1:T], [u[i];[zp[i,t],zq[i,t]]] in SecondOrderCone());
    @constraint(mp, powerflow_p[k in fData.brList, t in 1:T], p[k,t] == fData.g[k]*v[k[1],t] - fData.g[k]*vcos[k,t] - fData.b[k]*vsin[k,t]);
    @constraint(mp, powerflow_q[k in fData.brList, t in 1:T], q[k,t] == -fData.b[k]*v[k[1],t] + fData.b[k]*vcos[k,t] - fData.g[k]*vsin[k,t]);

    #@constraint(mp, powerflow[k in fData.brList, t in 1:T], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # set up the cuts, here tp is the disruption time
    for tp in 2:T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in 1:T]);
    @variable(mp,tAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3] >= 0);
    @variable(mp,tAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);
    @variable(mp,tAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);

    @constraint(mp,gcAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
        (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
        (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
        fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    @constraint(mp, genCost1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3],
        [tAux1[i,t];[tAux2[i,t],tAux3[i,t]]] in SecondOrderCone());
    @constraint(mp, genCost2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC1(T, fData, pDistr);
    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in 1:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in 2:T));

    @objective(mp, Min, objExpr);
    return mp;
end

# LinDistFlow relaxation
function build1_LF(Δt, T)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    # construct the first stage without disruption occurring
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0 ,"Threads" => 1));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
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
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in 1:T]);
    @variable(mp, q[k in fData.brList, t in 1:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]^2);
    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lpp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, pequal[k in fData.brList, t in 1:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in 1:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    @constraint(mp, lineThermal[k in fData.brList, t in 1:T], [fData.rateA[k];[p[k,t],q[k,t]]] in SecondOrderCone());
    @constraint(mp, bThermal[i in bData.IDList, t in 1:T], [u[i];[zp[i,t],zq[i,t]]] in SecondOrderCone());
    @constraint(mp, powerflow[k in fData.brList, t in 1:T], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # set up the cuts, here tp is the disruption time
    for tp in 2:T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in 1:T]);
    @variable(mp,tAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3] >= 0);
    @variable(mp,tAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);
    @variable(mp,tAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);

    @constraint(mp,gcAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
        (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
        (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
        fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    @constraint(mp, genCost1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3],
        [tAux1[i,t];[tAux2[i,t],tAux3[i,t]]] in SecondOrderCone());
    @constraint(mp, genCost2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC1(T, fData, pDistr);
    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in 1:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in 2:T));

    @objective(mp, Min, objExpr);
    return mp;
end

# original nonconvex QCQP
function build1_NC(Δt, T)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    # construct the first stage without disruption occurring
    mp = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0 ,"linear_solver" => "ma27"));

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
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
            sqhatsum[i,t] = @expression(mp,0.0);
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

    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lpp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    @variable(mp, va[i in fData.IDList, t in 1:T]);
    for i in fData.IDList
        for t in 1:T
            if fData.bType[i] == 3
                @constraint(mp, va[i,t] == 0);
            end
        end
    end

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @NLconstraint(mp, lineThermal[k in fData.brList, t in 1:T], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    @NLconstraint(mp, bThermal[i in bData.IDList, t in 1:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    @NLconstraint(mp, powerflow_p[k in fData.brList, t in 1:T], p[k,t] == fData.g[k]*v[k[1],t]^2 - fData.g[k]*v[k[1],t]*v[k[2],t]*cos(va[k[1],t] - va[k[2],t]) -
        fData.b[k]*v[k[1],t]*v[k[2],t]*sin(va[k[1],t] - va[k[2],t]));
    @NLconstraint(mp, powerflow_q[k in fData.brList, t in 1:T], q[k,t] == -fData.b[k]*v[k[1],t]^2 + fData.b[k]*v[k[1],t]*v[k[2],t]*cos(va[k[1],t] - va[k[2],t]) -
        fData.g[k]*v[k[1],t]*v[k[2],t]*sin(va[k[1],t] - va[k[2],t]));

    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # set up the cuts, here tp is the disruption time
    for tp in 2:T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in 1:T]);
    @NLconstraint(mp, genCost1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]^2 + fData.cp[i].params[2]*sp[i,t]);
    @NLconstraint(mp, genCost2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC1(T, fData, pDistr);
    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in 1:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in 2:T));

    @objective(mp, Min, objExpr);
    return mp;
end

# build the no disruption nominal case and solve it
function noDisruptionBuild(Δt, T, form_mode, solveOpt = true)
    # precalculate data
    if form_mode == "LinDistFlow"
        mp = build1_LF(Δt, T);
    elseif form_mode == "QC"
        mp = build1_QC(Δt, T);
    elseif form_mode == "SOCP"
        mp = build1_SOCP(Δt, T);
    else
        error("Wrong formulation category.")
    end

    if solveOpt
        # solve the problem
        stopIter = false;
        while !(stopIter)
            optimize!(mp);
            statusMp = JuMP.termination_status(mp);
            println("First stage, solving status $(statusMp)");
            if statusMp != MOI.OPTIMAL && statusMp != MOI.OTHER_LIMIT
                set_optimizer(mp, optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV),
                                                            "OutputFlag" => 0,
                                                            "NumericFocus" => 3,
                                                            "Threads" => 1));
            else
                stopIter = true;
            end
        end
        #with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", acceptable_tol = 1e-8, print_level = 0, max_iter = 10000));
        mpObj = objective_value(mp);
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        solzp = Dict();
        for i in fData.genIDList
            for t in 1:T
                if abs(value(mp[:sp][i,t])) > 1e-5
                    solSp[i,t] = value(mp[:sp][i,t]);
                else
                    solSp[i,t] = 0;
                end
                if abs(value(mp[:sq][i,t])) > 1e-5
                    solSq[i,t] = value(mp[:sq][i,t]);
                else
                    solSq[i,t] = 0;
                end
            end
        end
        for i in bData.IDList
            if abs(value(mp[:u][i])) > 1e-5
                solu[i] = value(mp[:u][i]);
                for t in 1:T
                    if abs(value(mp[:w][i,t])) > 1e-5
                        solw[i,t] = value(mp[:w][i,t]);
                    else
                        solw[i,t] = 0;
                    end
                    if abs(value(mp[:zp][i,t])) > 1e-6
                        solzp[i,t] = value(mp[:zp][i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            else
                solu[i] = 0;
                for t in 1:T
                    solw[i,t] = 0;
                    if abs(value(mp[:zp][i,t])) > 1e-6
                        solzp[i,t] = value(mp[:zp][i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            end
        end
        for i in fData.IDList
            for t in 1:T
                if (abs(value(mp[:lpp][i,t])) > 1e-8)|(abs(value(mp[:lpm][i,t])) > 1e-8)
                    solLp[i,t] = value(mp[:lpp][i,t]) - value(mp[:lpm][i,t]);
                else
                    solLp[i,t] = 0;
                end
                if (abs(value(mp[:lqp][i,t])) > 1e-8)|(abs(value(mp[:lqm][i,t])) > 1e-8)
                    solLq[i,t] = value(mp[:lqp][i,t]) - value(mp[:lqm][i,t]);
                else
                    solLq[i,t] = 0;
                end
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq,solzp);
        return sol,mpObj;
    else
        return mp;
    end
end

# QC relaxation for later time periods
function buildt_QC(td, ωd, currentSol, τ, Δt, T, hardened)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    Bparams = Dict();
    for t in td:T
        # create B parameters
        for k in fData.brList
            # if the line is disrupted and it is within disruption time
            if (((k[1],k[2]) in ωd)||((k[2],k[1]) in ωd))&&(t <= td + τ)
                if ωd in hardened
                    Bparams[k,t] = 1;
                else
                    Bparams[k,t] = 0;
                end
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i in ωd)&&(t <= td + τ)
                if i in hardened
                    Bparams[i,t] = 1;
                else
                    Bparams[i,t] = 0;
                end
            else
                Bparams[i,t] = 1;
            end
        end
    end
    for i in fData.genIDList
        Bparams[i,td - 1] = 1;
    end

    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0 ,"Threads" => 1));

    # set up the variables
    @variable(mp, fData.Pmin[i]*Bparams[i,t] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]*Bparams[i,t]);
    @variable(mp, fData.Qmin[i]*Bparams[i,t] <= sq[i in fData.genIDList,t in td:T] <= fData.Qmax[i]*Bparams[i,t]);
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    θu = pi/10;

    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i] <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]);

    @variable(mp, fData.Vmin[i]^2 <= vhat[i in fData.IDList, t in td:T] <= fData.Vmax[i]^2);
    @constraint(mp, vhatMock1[i in fData.IDList, t in td:T], [vhat[i,t] + 1/4; [vhat[i,t] - 1/4, v[i,t]]] in SecondOrderCone());
    @constraint(mp, vhatMock2[i in fData.IDList, t in td:T], vhat[i,t] - (fData.Vmax[i] + fData.Vmin[i])*v[i,t] <= -fData.Vmax[i]*fData.Vmin[i]);

    @variable(mp, vv[k in fData.brList, t in td:T]);
    @constraint(mp, vvEquality[k in fData.brList, t in td:T], vv[k,t] == vv[(k[2],k[1],k[3]),t]);

    @variable(mp, va[i in fData.IDList, t in td:T]);
    # set up the reference bus angle = 0
      for i in fData.IDList
          for t in td:T
              if fData.bType[i] == 3
                  @constraint(mp, va[i,t] == 0);
              end
          end
      end
    @constraint(mp, angleDiff[k in fData.brList, t in td:T], va[k[1],t] - va[k[2],t] <= θu);

    @variable(mp, cos(θu) <= cs[k in fData.brList, t in td:T] <= 1);
    @constraint(mp, csMock[k in fData.brList, t in td:T],
                    [5/4 - cs[k,t]; [sqrt((1-cos(θu))/(θu)^2)*(va[k[1],t] - va[k[2],t]), cs[k,t] - 3/4]] in SecondOrderCone());
    @constraint(mp, csEquality[k in fData.brList, t in td:T], cs[k,t] == cs[(k[2],k[1],k[3]),t]);

    @variable(mp, -sin(θu) <= ss[k in fData.brList, t in td:T] <= sin(θu));
    @constraint(mp, ssMock1[k in fData.brList, t in td:T], ss[k,t] <= cos(θu/2)*(va[k[1],t] - va[k[2],t] - θu/2) + sin(θu/2));
    @constraint(mp, ssMock2[k in fData.brList, t in td:T], -ss[k,t] <= -(cos(θu/2)*(va[k[1],t] - va[k[2],t] + θu/2) - sin(θu/2)));
    @constraint(mp, ssEquality[k in fData.brList, t in td:T], ss[k,t] == -ss[(k[2],k[1],k[3]),t]);

    @variable(mp, wc[k in fData.brList, t in td:T]);
    @variable(mp, ws[k in fData.brList, t in td:T]);
    @constraint(mp, wcEquality[k in fData.brList, t in td:T], wc[k,t] == wc[(k[2],k[1],k[3]),t]);
    @constraint(mp, wsEquality[k in fData.brList, t in td:T], ws[k,t] == -ws[(k[2],k[1],k[3]),t]);

    # McCormick Relaxation
    @constraint(mp, vvMcCormick1[k in fData.brList, t in td:T],
                vv[k,t] >= fData.Vmin[k[1]]*v[k[2],t] + fData.Vmin[k[2]]*v[k[1],t] - fData.Vmin[k[1]]*fData.Vmin[k[2]]);
    @constraint(mp, vvMcCormick2[k in fData.brList, t in td:T],
                vv[k,t] >= fData.Vmax[k[1]]*v[k[2],t] + fData.Vmax[k[2]]*v[k[1],t] - fData.Vmax[k[1]]*fData.Vmax[k[2]]);
    @constraint(mp, vvMcCormick3[k in fData.brList, t in td:T],
                vv[k,t] <= fData.Vmin[k[1]]*v[k[2],t] + fData.Vmax[k[2]]*v[k[1],t] - fData.Vmin[k[1]]*fData.Vmax[k[2]]);
    @constraint(mp, vvMcCormick4[k in fData.brList, t in td:T],
                vv[k,t] <= fData.Vmax[k[1]]*v[k[2],t] + fData.Vmin[k[2]]*v[k[1],t] - fData.Vmax[k[1]]*fData.Vmin[k[2]]);

    @constraint(mp, wcMcCormick1[k in fData.brList, t in td:T],
                wc[k,t] >= fData.Vmin[k[1]]*fData.Vmin[k[2]]*cs[k,t] + cos(θu)*vv[k,t] - fData.Vmin[k[1]]*fData.Vmin[k[2]]*cos(θu));
    @constraint(mp, wcMcCormick2[k in fData.brList, t in td:T],
                wc[k,t] >= fData.Vmax[k[1]]*fData.Vmax[k[2]]*cs[k,t] + 1*vv[k,t] - fData.Vmax[k[1]]*fData.Vmax[k[2]]*1);
    @constraint(mp, wcMcCormick3[k in fData.brList, t in td:T],
                wc[k,t] <= fData.Vmin[k[1]]*fData.Vmin[k[2]]*cs[k,t] + vv[k,t]*1 - fData.Vmin[k[1]]*fData.Vmin[k[2]]*1);
    @constraint(mp, wcMcCormick4[k in fData.brList, t in td:T],
                wc[k,t] <= fData.Vmax[k[1]]*fData.Vmax[k[2]]*cs[k,t] + vv[k,t]*cos(θu) - fData.Vmax[k[1]]*fData.Vmax[k[2]]*cos(θu));

    @constraint(mp, wsMcCormick1[k in fData.brList, t in td:T],
                ws[k,t] >= fData.Vmin[k[1]]*fData.Vmin[k[2]]*ss[k,t] + (-sin(θu))*vv[k,t] - fData.Vmin[k[1]]*fData.Vmin[k[2]]*(-sin(θu)));
    @constraint(mp, wsMcCormick2[k in fData.brList, t in td:T],
                ws[k,t] >= fData.Vmax[k[1]]*fData.Vmax[k[2]]*ss[k,t] + sin(θu)*vv[k,t] - fData.Vmax[k[1]]*fData.Vmax[k[2]]*sin(θu));
    @constraint(mp, wsMcCormick3[k in fData.brList, t in td:T],
                ws[k,t] <= fData.Vmin[k[1]]*fData.Vmin[k[2]]*ss[k,t] + vv[k,t]*sin(θu) - fData.Vmin[k[1]]*fData.Vmin[k[2]]*sin(θu));
    @constraint(mp, wsMcCormick4[k in fData.brList, t in td:T],
                ws[k,t] <= fData.Vmax[k[1]]*fData.Vmax[k[2]]*ss[k,t] + vv[k,t]*(-sin(θu)) - fData.Vmax[k[1]]*fData.Vmax[k[2]]*(-sin(θu)));

    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lpp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], [fData.rateA[k],p[k,t],q[k,t]] in SecondOrderCone());
    @constraint(mp, bThermal[i in bData.IDList, t in td:T], [u[i],zp[i,t],zq[i,t]] in SecondOrderCone());
    @constraint(mp, lineThermal2[k in fData.brList, t in td:T;Bparams[k,t] == 0], p[k,t] == 0);
    @constraint(mp, lineThermal3[k in fData.brList, t in td:T;Bparams[k,t] == 0], q[k,t] == 0);
    @constraint(mp, powerflow_p[k in fData.brList, t in td:T; Bparams[k,t] == 1], p[k,t] == fData.g[k]*vhat[k[1],t] - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]);
    @constraint(mp, powerflow_q[k in fData.brList, t in td:T; Bparams[k,t] == 1], q[k,t] == -fData.b[k]*vhat[k[1],t] + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]);

    # @constraint(mp, powerflow1[k in fData.brList, t in td:T;Bparams[k,t] == 1], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the cuts, here tp is the disruption time
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in td:T]);
    @variable(mp,tAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3] >= 0);
    @variable(mp,tAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);
    @variable(mp,tAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);

    @constraint(mp,gcAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
        (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
        (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
        fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    @constraint(mp, genCost1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3],
        [tAux1[i,t],tAux2[i,t],tAux3[i,t]] in SecondOrderCone());
    @constraint(mp, genCost2[i in fData.genIDList, t in td:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC(td, τ, T, fData, pDistr);
    objExpr = @expression(mp, sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in td:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in (td + τ + 1):T));
    @objective(mp, Min, objExpr);
    return mp;
end

# SOCP relaxation for later time periods
function buildt_SOCP(td, ωd, currentSol, τ, Δt, T, hardened)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    Bparams = Dict();
    for t in td:T
        # create B parameters
        for k in fData.brList
            # if the line is disrupted and it is within disruption time
            if (((k[1],k[2]) in ωd)||((k[2],k[1]) in ωd))&&(t <= td + τ)
                if ωd in hardened
                    Bparams[k,t] = 1;
                else
                    Bparams[k,t] = 0;
                end
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i in ωd)&&(t <= td + τ)
                if i in hardened
                    Bparams[i,t] = 1;
                else
                    Bparams[i,t] = 0;
                end
            else
                Bparams[i,t] = 1;
            end
        end
    end
    for i in fData.genIDList
        Bparams[i,td - 1] = 1;
    end

    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0 ,"Threads" => 1));

    # set up the variables
    @variable(mp, fData.Pmin[i]*Bparams[i,t] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]*Bparams[i,t]);
    @variable(mp, fData.Qmin[i]*Bparams[i,t] <= sq[i in fData.genIDList,t in td:T] <= fData.Qmax[i]*Bparams[i,t]);
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]^2);

    @variable(mp, vcos[k in fData.brList, t in td:T] >= 0);
    @variable(mp, vsin[k in fData.brList, t in td:T]);
    @constraint(mp, vcSymm[k in fData.brList, t in td:T], vcos[k,t] == vcos[(k[2],k[1],k[3]),t]);
    @constraint(mp, vsSymm[k in fData.brList, t in td:T], vsin[k,t] == -vsin[(k[2],k[1],k[3]),t]);
    @constraint(mp, socpCon[k in fData.brList, t in td:T], [1/sqrt(2)*(v[k[1],t] + v[k[2],t]); [1/sqrt(2)*v[k[1],t], 1/sqrt(2)*v[k[2],t],
                                                            vcos[k,t], vsin[k,t]]] in SecondOrderCone());

    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lpp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], [fData.rateA[k],p[k,t],q[k,t]] in SecondOrderCone());
    @constraint(mp, bThermal[i in bData.IDList, t in td:T], [u[i],zp[i,t],zq[i,t]] in SecondOrderCone());
    @constraint(mp, lineThermal2[k in fData.brList, t in td:T;Bparams[k,t] == 0], p[k,t] == 0);
    @constraint(mp, lineThermal3[k in fData.brList, t in td:T;Bparams[k,t] == 0], q[k,t] == 0);
    @constraint(mp, powerflow_p[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t] == fData.g[k]*v[k[1],t] - fData.g[k]*vcos[k,t] - fData.b[k]*vsin[k,t]);
    @constraint(mp, powerflow_q[k in fData.brList, t in td:T;Bparams[k,t] == 1], q[k,t] == -fData.b[k]*v[k[1],t] + fData.b[k]*vcos[k,t] - fData.g[k]*vsin[k,t]);

    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the cuts, here tp is the disruption time
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in td:T]);
    @variable(mp,tAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3] >= 0);
    @variable(mp,tAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);
    @variable(mp,tAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);

    @constraint(mp,gcAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
        (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
        (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
        fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    @constraint(mp, genCost1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3],
        [tAux1[i,t],tAux2[i,t],tAux3[i,t]] in SecondOrderCone());
    @constraint(mp, genCost2[i in fData.genIDList, t in td:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC(td, τ, T, fData, pDistr);
    objExpr = @expression(mp, sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in td:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in (td + τ + 1):T));
    @objective(mp, Min, objExpr);
    return mp;
end

# LinDistFlow relaxation for later time periods
function buildt_LF(td, ωd, currentSol, τ, Δt, T, hardened)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    Bparams = Dict();
    for t in td:T
        # create B parameters
        for k in fData.brList
            # if the line is disrupted and it is within disruption time
            if (((k[1],k[2]) in ωd)||((k[2],k[1]) in ωd))&&(t <= td + τ)
                if ωd in hardened
                    Bparams[k,t] = 1;
                else
                    Bparams[k,t] = 0;
                end
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i in ωd)&&(t <= td + τ)
                if i in hardened
                    Bparams[i,t] = 1;
                else
                    Bparams[i,t] = 0;
                end
            else
                Bparams[i,t] = 1;
            end
        end
    end
    for i in fData.genIDList
        Bparams[i,td - 1] = 1;
    end

    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0 ,"Threads" => 1));

    # set up the variables
    @variable(mp, fData.Pmin[i]*Bparams[i,t] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]*Bparams[i,t]);
    @variable(mp, fData.Qmin[i]*Bparams[i,t] <= sq[i in fData.genIDList,t in td:T] <= fData.Qmax[i]*Bparams[i,t]);
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]^2);
    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lpp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, pequal[k in fData.brList, t in td:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in td:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], [fData.rateA[k],p[k,t],q[k,t]] in SecondOrderCone());
    @constraint(mp, bThermal[i in bData.IDList, t in td:T], [u[i],zp[i,t],zq[i,t]] in SecondOrderCone());
    @constraint(mp, lineThermal2[k in fData.brList, t in td:T;Bparams[k,t] == 0], p[k,t] == 0);
    @constraint(mp, lineThermal3[k in fData.brList, t in td:T;Bparams[k,t] == 0], q[k,t] == 0);
    @constraint(mp, powerflow1[k in fData.brList, t in td:T;Bparams[k,t] == 1], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the cuts, here tp is the disruption time
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in td:T]);
    @variable(mp,tAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3] >= 0);
    @variable(mp,tAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);
    @variable(mp,tAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);

    @constraint(mp,gcAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
        (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
        (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    @constraint(mp,gcAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
        fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    @constraint(mp, genCost1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3],
        [tAux1[i,t],tAux2[i,t],tAux3[i,t]] in SecondOrderCone());
    @constraint(mp, genCost2[i in fData.genIDList, t in td:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC(td, τ, T, fData, pDistr);
    objExpr = @expression(mp, sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in td:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in (td + τ + 1):T));
    @objective(mp, Min, objExpr);
    return mp;
end

# SOCP relaxation for later time periods
function buildt_NC(td, ωd, currentSol, τ, Δt, T, hardened)
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    Bparams = Dict();
    for t in td:T
        # create B parameters
        for k in fData.brList
            # if the line is disrupted and it is within disruption time
            if (((k[1],k[2]) in ωd)||((k[2],k[1]) in ωd))&&(t <= td + τ)
                if ωd in hardened
                    Bparams[k,t] = 1;
                else
                    Bparams[k,t] = 0;
                end
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i in ωd)&&(t <= td + τ)
                if i in hardened
                    Bparams[i,t] = 1;
                else
                    Bparams[i,t] = 0;
                end
            else
                Bparams[i,t] = 1;
            end
        end
    end
    for i in fData.genIDList
        Bparams[i,td - 1] = 1;
    end

    mp = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0 ,"linear_solver" => "ma27"));

    # set up the variables
    @variable(mp, fData.Pmin[i]*Bparams[i,t] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]*Bparams[i,t]);
    @variable(mp, fData.Qmin[i]*Bparams[i,t] <= sq[i in fData.genIDList,t in td:T] <= fData.Qmax[i]*Bparams[i,t]);
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i] <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]);

    @variable(mp, va[i in fData.IDList, t in 1:T]);
    for i in fData.IDList
        for t in 1:T
            if fData.bType[i] == 3
                @constraint(mp, va[i,t] == 0);
            end
        end
    end

    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lpp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @NLconstraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    @NLconstraint(mp, bThermal[i in bData.IDList, t in td:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    @constraint(mp, lineThermal2[k in fData.brList, t in td:T;Bparams[k,t] == 0], p[k,t] == 0);
    @constraint(mp, lineThermal3[k in fData.brList, t in td:T;Bparams[k,t] == 0], q[k,t] == 0);
    @NLconstraint(mp, powerflow_p[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t] == fData.g[k]*v[k[1],t]^2 - fData.g[k]*v[k[1],t]*v[k[2],t]*cos(va[k[1],t] - va[k[2],t]) -
        fData.b[k]*v[k[1],t]*v[k[2],t]*sin(va[k[1],t] - va[k[2],t]));
    @NLconstraint(mp, powerflow_q[k in fData.brList, t in td:T;Bparams[k,t] == 1], q[k,t] == -fData.b[k]*v[k[1],t]^2 + fData.b[k]*v[k[1],t]*v[k[2],t]*cos(va[k[1],t] - va[k[2],t]) -
        fData.g[k]*v[k[1],t]*v[k[2],t]*sin(va[k[1],t] - va[k[2],t]));

    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the cuts, here tp is the disruption time
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in td:T]);
    @NLconstraint(mp, genCost1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]^2 + fData.cp[i].params[2]*sp[i,t]);
    @constraint(mp, genCost2[i in fData.genIDList, t in td:T; fData.cp[i].n == 2],
        fs[i,t] == fData.cp[i].params[1]*sp[i,t]);

    fDict,lDict,θDict = calDualC(td, τ, T, fData, pDistr);
    objExpr = @expression(mp, sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in td:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in (td + τ + 1):T));
    @objective(mp, Min, objExpr);
    return mp;
end

function fBuild(td, ωd, currentSol, τ, Δt, T, form_mode, solveOpt = true, hardened = [])
    # precalculate data
    if form_mode == "LinDistFlow"
        mp = buildt_LF(td, ωd, currentSol, τ, Δt, T, hardened);
    elseif form_mode == "QC"
        mp = buildt_QC(td, ωd, currentSol, τ, Δt, T, hardened);
    elseif form_mode == "SOCP"
        mp = buildt_SOCP(td, ωd, currentSol, τ, Δt, T, hardened);
    else
        error("Wrong formulation category.")
    end

    if solveOpt
        stopIter = false;
        while !(stopIter)
            optimize!(mp);
            statusMp = JuMP.termination_status(mp);
            println("Disruption time $(td), scenario $(ωd), solving status $(statusMp)");
            if statusMp != MOI.OPTIMAL && statusMp != MOI.OTHER_LIMIT
                set_optimizer(mp, optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV),
                                                            "OutputFlag" => 0,
                                                            "NumericFocus" => 3,
                                                            "Threads" => 1));
            else
                stopIter = true;
            end
        end
        #optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", acceptable_tol = 1e-8, print_level = 0, max_iter = 10000));
        mpObj = objective_value(mp);
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        solzp = Dict();
        for i in fData.genIDList
            for t in td:T
                if abs(value(mp[:sp][i,t])) > 1e-5
                    solSp[i,t] = value(mp[:sp][i,t]);
                else
                    solSp[i,t] = 0;
                end
                if abs(value(mp[:sq][i,t])) > 1e-5
                    solSq[i,t] = value(mp[:sq][i,t]);
                else
                    solSq[i,t] = 0;
                end
            end
        end
        for i in bData.IDList
            if abs(value(mp[:u][i])) > 1e-5
                solu[i] = value(mp[:u][i]);
                for t in td:T
                    if abs(value(mp[:w][i,t])) > 1e-5
                        solw[i,t] = value(mp[:w][i,t]);
                    else
                        solw[i,t] = 0;
                    end
                    if abs(value(mp[:zp][i,t])) > 1e-6
                        solzp[i,t] = value(mp[:zp][i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            else
                solu[i] = 0;
                for t in td:T
                    solw[i,t] = 0;
                    if abs(value(mp[:zp][i,t])) > 1e-6
                        solzp[i,t] = value(mp[:zp][i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            end
        end
        for i in fData.IDList
            for t in td:T
                if (abs(value(mp[:lpp][i,t])) > 1e-8)|(abs(value(mp[:lpm][i,t])) > 1e-8)
                    solLp[i,t] = value(mp[:lpp][i,t]) - value(mp[:lpm][i,t]);
                else
                    solLp[i,t] = 0;
                end
                if (abs(value(mp[:lqp][i,t])) > 1e-8)|(abs(value(mp[:lqm][i,t]))> 1e-8)
                    solLq[i,t] = value(mp[:lqp][i,t]) - value(mp[:lqm][i,t]);
                else
                    solLq[i,t] = 0;
                end
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq,solzp);
        return sol,mpObj;
    else
        return mp;
    end
end

function constructForwardM(td, ωd, sol, Δt, T, τ, form_mode, hardenend = [])
    # construct the math program given the state variables and current stage
    if td == 1
        # if it is the no-disruption problem
        sol,objV = noDisruptionBuild(Δt, T, form_mode);
    else
        # if it is f_{ht}^ω
        sol,objV = fBuild(td, ωd, sol, τ, Δt, T, form_mode, true, hardenend);
    end
    return sol,objV;
end

function buildPath(T, Δt, form_mode, pathList = [], hardened = [])
    disT = 1;
    ωd = 0;
    τω = 0;
    costn = 0;
    solHist = [];
    currentLB = 0;
    currentSol = solData(Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict());
    iter = 1;
    while disT <= T
        # solve the current stage problem, state variables are passed
        nowT = disT;
        currentSol,objV = constructForwardM(disT, ωd, currentSol, Δt, T, τω, form_mode, hardened);
        push!(solHist,(currentSol,nowT,ωd,τω));

        # generate disruption
        if pathList == []
            tp,ωd,τω = genScenario(pDistr);
        else
            tp,ωd,τω = pathList[iter];
        end
        iter += 1;
        if nowT == 1
            currentLB = objV;
            disT += tp;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn += sum(currentSol.u[i]*bData.cost[i] for i in bData.IDList);
            costn = calCostF(costn, currentSol, T, fData, nowT, disT);
        else
            disT += tp + τω;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn = calCostF(costn, currentSol, T, fData, nowT, disT);
        end
    end
    println("Path Built!");
    return [solHist,currentLB,costn];
end

function exeForward(T, Δt, N, form_mode, pathDict = Dict(), hardened = [])
    # execution of forward pass
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts (global in every core)
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    objV = 0;
    currentLB = 0;
    # for n in 1:N
    #     # for each trial path
    #     returnData = buildPath(τ, T, Δt, fData, bData, dData, pDistr, cutDict);
    #     solDict[n] = returnData[1];
    #     costDict[n] = returnData[2];
    # end
    if pathDict == Dict()
        for i in 1:N
            pathDict[i] = [];
        end
    end
    returnData = pmap(i -> buildPath(T, Δt, form_mode, pathDict[i], hardened), 1:N);
    for n in 1:N
        solDict[n] = returnData[n][1];
        costDict[n] = returnData[n][3];
    end
    currentLB = returnData[1][2];
    return solDict, currentLB, costDict;
end
