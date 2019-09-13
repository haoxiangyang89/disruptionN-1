# define the type of data

struct fixedData
    # static network data
    baseMVA :: Float64
    bType :: Dict{Int64,Any}

    IDList :: Array{Int64,1}
    genIDList :: Array{Int64,1}
    brList :: Array{Any,1}
    brRev :: Dict{Any,Any}

    Loc :: Dict{Int64,Any}
    LocRev :: Dict{Int64,Any}
    Vmax :: Dict{Int64,Any}
    Vmin :: Dict{Int64,Any}
    Pmax :: Dict{Int64,Any}
    Pmin :: Dict{Int64,Any}
    Qmax :: Dict{Int64,Any}
    Qmin :: Dict{Int64,Any}
    gs :: Dict{Int64,Any}
    bs :: Dict{Int64,Any}
    Vmag :: Dict{Int64,Any}
    Vang :: Dict{Int64,Any}
    Pd :: Dict{Int64,Any}
    Qd :: Dict{Int64,Any}
    Pg :: Dict{Int64,Any}
    Qg :: Dict{Int64,Any}
    RU :: Dict{Int64,Any}
    RD :: Dict{Int64,Any}

    g :: Dict{Tuple{Int64,Int64,Int64},Any}
    b :: Dict{Tuple{Int64,Int64,Int64},Any}
    bc :: Dict{Tuple{Int64,Int64,Int64},Any}
    θmax :: Dict{Tuple{Int64,Int64,Int64},Any}
    θmin :: Dict{Tuple{Int64,Int64,Int64},Any}
    rateA :: Dict{Tuple{Int64,Int64,Int64},Any}
    τ1 :: Dict{Tuple{Int64,Int64,Int64},Any}
    τ2 :: Dict{Tuple{Int64,Int64,Int64},Any}
    σ :: Dict{Tuple{Int64,Int64,Int64},Any}

    cp :: Dict{Any,Any}
    cq :: Dict{Any,Any}
    cz :: Float64

    busInd :: Dict{Any,Any}
    branchDict1 :: Dict{Any,Any}
    branchDict2 :: Dict{Any,Any}
    connectPair :: Array{Any,1}
    connectDict :: Dict{Any,Any}
    kpDict :: Dict{Any,Any}

end

struct costDataType
    model :: Int64
    upCost :: Float64
    downCost :: Float64
    n :: Int64
    params :: Array{Float64,1}
end

struct probDistrn
    # probability distribution of disruption timing
    tDistrn :: Dict{Int64,Any}

    # probability distribution of disruption magnitude
    ωDistrn :: Dict{Any,Any}
end

struct batteryData
    # battery information: charging/discharging factor, capacity
    IDList :: Array{Any,1}
    Loc :: Dict{Int64,Any}
    bInv :: Dict{Int64,Any}
    ηα :: Dict{Int64,Any}
    ηβ :: Dict{Int64,Any}
    cap :: Dict{Int64,Any}
    cost :: Dict{Int64,Any}
end

struct solData
    # solution information: generaion level, storage level, maximum generation output
    sp :: Dict{Any,Any}
    sq :: Dict{Any,Any}
    w :: Dict{Any,Any}
    u :: Dict{Any,Any}
    lp :: Dict{Any,Any}
    lq :: Dict{Any,Any}
end

struct demandData
    # demand data: active demand and reactive demand
    T :: Int64
    pd :: Dict{Any,Any}
    qd :: Dict{Any,Any}
end

struct cutData
    # cut data: slope for sp, w and u and the intercept
    λ :: Dict{Any,Any}
    γ :: Dict{Any,Any}
    μ :: Dict{Any,Any}
    vhat :: Float64
    sphat :: Dict{Any,Any}
    what :: Dict{Any,Any}
    uhat :: Dict{Any,Any}
end
