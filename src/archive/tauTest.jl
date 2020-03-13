# test 4a: effectiveness of different τ
include("run.jl");

cutDictT = Dict();
LBHistT = Dict();
UBHistT = Dict();
UBuHistT = Dict();
UBlHistT = Dict();
for τ in 0:Int64(round(T/2))
    cutDictT[τ],LBHistT[τ],UBHistT[τ],UBuHistT[τ],UBlHistT[τ] =
        solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
    println("Disruption Length $(τ): LB = $(LBHistT[τ][length(LBHistT[τ])])");
end

outputData = [cutDictT,LBHistT,UBHistT,UBuHistT,UBlHistT];
save("tauOut.jld","data",outputData);
