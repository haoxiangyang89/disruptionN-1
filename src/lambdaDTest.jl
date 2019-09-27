# test 4b: effectiveness of different λ_D
include("run.jl");

cutDictλ = Dict();
LBHistλ = Dict();
UBHistλ = Dict();
UBuHistλ = Dict();
UBlHistλ = Dict();
λList = 1 ./ ([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]*T);
for λD in λList
    # generate the new tDistrn
    pDistrNew = modifyT(pDistr,λD,T);

    # run the test for the current lambda
    cutDictλ[λD],LBHistλ[λD],UBHistλ[λD],UBuHistλ[λD],UBlHistλ[λD] =
        solveMain(τ, T, Δt, fData, pDistrNew, bData, dData, N);
    println("Disruption Frequency $(λD): LB = $(LBHistλ[λD][length(LBHistλ[λD])])");
end

outputData = [cutDictλ,LBHistλ,UBHistλ,UBuHistλ,UBlHistλ];
save("tauOut.jld","data",outputData);
