# define the type of data

type fixedData
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

    busInd :: Dict{Any,Any}
end

type probDistrn
    # probability distribution of disruption timing
    tDistrn :: Dict{Int64,Any}

    # probability distribution of disruption magnitude
    ωDistrn :: Dict{Int64,Any}
end

type batteryData
    # battery information: charging/discharging factor, capacity
    ηd :: Dict{Int64,Any}
    ηc :: Dict{Int64,Any}
    cap :: Dict{Int64,Any}
    cost :: Dict{Int64,Any}
end
