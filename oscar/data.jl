###
### Utils for reading in the problem data.
###

to_ω(x::Int) = x

function to_ω(x::AbstractString)
    m = match(r"\((\d+).+?(\d+)\)", x)
    @assert m !== nothing
    return parse(Int, m[1]) => parse(Int, m[2])
end

function build_data(ieee_file, B_file, P_file, Q_file, prob_file)
    A = DelimitedFiles.readdlm(prob_file, ',')
    T = size(A, 2)
    ρ = Statistics.mean(A[2, 1:(T - 2)] ./ A[2, 2:(T - 1)])
    support = Union{Int, Pair{Int, Int}}[]
    nominal_probability = Float64[]
    for t = 1:T
        if isempty(A[3, t])
            break
        end
        push!(support, to_ω(A[3, t]))
        push!(nominal_probability, A[4, t])
    end
    return (
        network_data = PowerModels.parse_file(ieee_file),
        B = DelimitedFiles.readdlm(B_file, ','),
        P = DelimitedFiles.readdlm(P_file, ',')[:, 2:end],
        Q = DelimitedFiles.readdlm(Q_file, ',')[:, 2:end],
        T = T,
        ρ = ρ,
        support = support,
        nominal_probability = nominal_probability,
    )
end

###
### Construct the structure of the policy graph.
###

struct DisruptionNode
    stage::Int
    disruption_length::Int
end

function SDDP.sort_nodes(x::Vector{DisruptionNode})
    return sort!(x; by = (y) -> (y.stage, y.disruption_length))
end

function build_graph(T, ρ, N)
    graph = SDDP.Graph(DisruptionNode(0, 0))
    for t = 1:T
        SDDP.add_node(graph, DisruptionNode(t, 0))
        nt = min(N, T - t + 1)
        SDDP.add_node(graph, DisruptionNode(t, nt))
        SDDP.add_edge(
            graph, DisruptionNode(t - 1, 0) => DisruptionNode(t, 0), 1 - ρ
        )
        SDDP.add_edge(
            graph, DisruptionNode(t - 1, 0) => DisruptionNode(t, nt), ρ
        )
        if t > N
            SDDP.add_edge(
                graph, DisruptionNode(t - N, N) => DisruptionNode(t, 0), 1.0
            )
        end
    end
    return graph
end

function readBattery(dataRaw,baseMVA = 100)
    # read in the battery information: charging/discharging factor, capacity
    # csv file format:
    # First column: node ID
    # Second column: ηd
    # Third column: ηc
    # Fourth column: capacity
    # Fifth column: cost (optional)
    mb,nb = size(dataRaw);
    capacity = Dict();
    cost = Dict();
    ηα = Dict();
    ηβ = Dict();
    IDList = [];
    LocDict = Dict();
    bInv = Dict();
    uCap = Dict();

    for i in 1:mb
        ID = Int64(dataRaw[i,1]);
        push!(IDList,ID);
        loc = Int64(dataRaw[i,2]);
        LocDict[ID] = loc;
        capacity[ID] = dataRaw[i,3]/baseMVA;
        cost[ID] = dataRaw[i,4]*baseMVA;
        bInv[ID] = dataRaw[i,5]/baseMVA;
        uCap[ID] = dataRaw[i,6]/baseMVA;
        ηparams = [j for j in dataRaw[i,7:nb] if j != ""];
        ηα[ID] = [];
        ηβ[ID] = [];
        for j in 1:2:length(ηparams)
            push!(ηα[ID],ηparams[j]);
            # ηα must be monotonously decreasing
            push!(ηβ[ID],ηparams[j + 1]);
        end
    end
    bData = (IDList = IDList,
            LocDict = LocDict,
            bInv = bInv,
            ηα = ηα,
            ηβ = ηβ,
            capacity = capacity,
            cost = cost,
            uCap = uCap);
    return bData;
end
