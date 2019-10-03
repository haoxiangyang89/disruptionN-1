# test 3 hardening tests
include("run.jl");

# obtain the no hardening costs for comparison
cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N);
println("No Hardening: LB = $(LBHist[length(LBHist)])");

cutDictH = Dict();
LBHistH = Dict();
UBHistH = Dict();
UBuHistH = Dict();
UBlHistH = Dict();
for comp in keys(pDistr.ωDistrn)
    pDistrNew = modifyOmega(pDistr,comp);
    cutDictH[comp],LBHistH[comp],UBHistH[comp],UBuHistH[comp],UBlHistH[comp] =
        solveMain(τ, T, Δt, fData, pDistrNew, bData, dData, N);
    println("Hardening Component $(comp): LB = $(LBHistH[comp][length(LBHistH[comp])])");
end

outputData = [cutDict,LBHist,UBHist,UBuHist,UBlHist,
    cutDictH,LBHistH,UBHistH,UBuHistH,UBlHistH];
save("hardOut.jld","data",outputData);
