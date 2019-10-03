dData = readDemand(pAdd,qAdd,"csv");

dData.pd[2] = dData.pd[2]*0.5;
dData.qd[2] = dData.qd[2]*0.5;

dData.pd[4] = dData.pd[4]*0.3;
dData.qd[4] = dData.qd[4]*0.3;

dData.pd[8] = dData.pd[8]*0.25;
dData.qd[8] = dData.qd[8]*0.3;

dData.pd[10] = dData.pd[10]*0.3;
dData.qd[10] = dData.qd[10]*0.3;

pdNew = zeros(length(fData.IDList),length(dData.pd[1])+1);
qdNew = zeros(length(fData.IDList),length(dData.pd[1])+1);
for i in fData.IDList
   pdNew[i,1] = i;
   qdNew[i,1] = i;
   for t in 1:length(dData.pd[i])
      pdNew[i,t+1] = round(dData.pd[i][t],digits = 4);
      qdNew[i,t+1] = round(dData.qd[i][t],digits = 4);
   end
end
writedlm("./test/testDataP_96_Fixed.csv",pdNew,',');
writedlm("./test/testDataQ_96_Fixed.csv",qdNew,',')

solDet,vDet = detBuild(Δt, T, fData, bData, dData);
println(vDet);

for i in 1:13
   for t in 1:96
      if (solDet.lp[i,t] > 0)|(solDet.lq[i,t] > 0)
         println(i," ",t," ",solDet.lp[i,t]," ",solDet.lq[i,t]);
      end
   end
end

mpDet = detBuild(Δt, T, fData, bData, dData, false);
optimize!(mpDet,with_optimizer(Gurobi.Optimizer));
# examine the line load
bigLoad = [];
for t in 1:T
   for k in fData.brList
      pflow = value(mpDet[:p][k,t]);
      qflow = value(mpDet[:q][k,t]);
      if pflow^2 + qflow^2 >= 0.8*fData.rateA[k]^2
         println(k," ",t," ",(pflow^2 + qflow^2)/fData.rateA[k]^2);
         if !(((k[1],k[2]) in bigLoad)|((k[2],k[1]) in bigLoad))
            push!(bigLoad,(k[1],k[2]));
         end
      end
   end
end
