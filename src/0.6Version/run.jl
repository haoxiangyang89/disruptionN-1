# running scripts
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
#pmap(i -> importIpopt(),1:30);

caseList = [13,33,123];

T = 96;
τ = T/6;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,1/4,T);

τ = 4;
T = 24;
Δt = 0.25;
N = 30;
pDistr = modifyT(pDistr,4/T,T);
