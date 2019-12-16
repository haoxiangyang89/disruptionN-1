# output csv files
@everywhere include("loadMod.jl");
caseList = [13,33,123];

# output the csv file for GenAll and dOnly
for ci in 1:length(caseList)
    datad = load("dOnlyResults_$(ci)_96.jld");
    datapg = load("pgResults_$(ci).jld");
    dataag = [datad["data"]["allGen"][1][1:20];zeros(80,1)];
    lbOutList = [dataag datapg["data"]["dOnly"][1][1:100] datapg["data"]["preGen"][1][1:100]];
    timeOutList = zeros(101,3);
    cT1 = 0;
    cT2 = 0;
    cT3 = 0;
    for j in 1:100
        cT2 += datapg["data"]["dOnly"][5][j];
        cT3 += datapg["data"]["preGen"][5][j];
        if j <= 20
            cT1 += datad["data"]["allGen"][5][j];
        end
        timeOutList[j,1] = cT1;
        timeOutList[j,2] = cT2;
        timeOutList[j,3] = cT3;
    end
    writedlm("./csvOut/dOnly_lb_$(ci).csv",lbOutList,',');
    writedlm("./csvOut/dOnly_time_$(ci).csv",timeOutList,',');
end

# output the csv file for preGen
for ci in 1:length(caseList)
    datapg = load("pgResults_$(ci).jld");
    lbOutList = [datapg["data"]["dOnly"][1][1:60] datapg["data"]["preGen"][1][1:60]];
    timeOutList = zeros(61,2);
    cT1 = 0;
    cT2 = datapg["data"]["preGen"][5][7];
    timeOutList[1,2] = cT2;
    for j in 1:60
        cT1 += datapg["data"]["dOnly"][5][j];
        cT2 += datapg["data"]["preGen"][5][j];
        timeOutList[j+1,1] = cT1;
        timeOutList[j+1,2] = cT2;
    end
    writedlm("./csvOut/pg_lb_$(ci).csv",lbOutList,',');
    writedlm("./csvOut/pg_time_$(ci).csv",timeOutList,',');
end

# output the csv file for NTest
NList = [1,2,3,4,5,10,20];
NNoList = [100,50,33,25,20,10,5];
lbDict = Dict();
timeDict = Dict();
dataN = Dict();
for ci in 1:length(caseList)
    dataN[ci] = Dict();
    for Nind in 1:length(NList)
        dataN[ci][NList[Nind]] = load("NResults_$(ci)_AG_$(NList[Nind]).jld");
    end
end
for Nind in 1:length(NList)
    N = NList[Nind];
    lbDict[N] = [];
    timeDict[N] = [];
    NNo = NNoList[Nind];
    for ci in 1:length(caseList)
        println(N, " ", ci);
        if ci != 1
            lbDict[N] = [lbDict[N] dataN[ci][N]["NOut"][N][1][1:NNo]];
            tdTemp = zeros(NNo+1);
            CT = dataN[ci][N]["NOut"][N][7];
            for j in 1:NNo
                CT += dataN[ci][N]["NOut"][N][5][j];
                tdTemp[j+1] = CT;
            end
            timeDict[N] = [timeDict[N] tdTemp];
        else
            lbDict[N] = dataN[ci][N]["NOut"][N][1][1:NNo];
            tdTemp = zeros(NNo+1);
            CT = dataN[ci][N]["NOut"][N][7];
            for j in 1:NNo
                CT += dataN[ci][N]["NOut"][N][5][j];
                tdTemp[j+1] = CT;
            end
            timeDict[N] = tdTemp;
        end
    end
    writedlm("./csvOut/N_lb_$(Nind).csv",lbDict[N],',');
    writedlm("./csvOut/N_time_$(Nind).csv",timeDict[N],',');
end

# print out the datadet
TList = [24,36,48,72,96];
for i in 1:3
    datadet[i] = load("detResults_$(i)_AG.jld");
end
for i in 1:3
    for T in TList
        print(" & ",round(datadet[i]["detOut"][T][3],1));
    end
    println("\\\\");
end
for i in 1:3
    for T in TList
        print(" & ",round(datadet[i]["stochOut"][T][4],1));
    end
    println("\\\\");
end
for i in 1:3
    for T in TList
        print(" & ",round(datadet[i]["nomOut"][T]));
    end
    println("\\\\");
end

# print out the dataStab
dataStab = Dict();
for i in 1:20
    dataStab[i] = load("stabilityResults_$(i)_AG.jld");
end
for i in 1:3
    LBMean = [];
    LBStdEv = [];
    UBMean = [];
    UBStdEv = [];
    for T in [24,36,48,72]
        LBList = zeros(20);
        UBList = zeros(20);
        for j in 1:20
            LBList[j] = dataStab[j]["data"][i][T][j][1];
            UBList[j] = dataStab[j]["data"][i][T][j][4];
        end
        push!(LBMean,mean(LBList));
        push!(LBStdEv,std(LBList)*1.96);
        push!(UBMean,mean(UBList));
        push!(UBStdEv,std(UBList)*1.96);
    end
    for tint in 1:length(LBMean)
        print(" & \$", round(LBMean[tint],1), " \\pm ", round(LBStdEv[tint],1), "\$");
    end
    println("\\\\");
    for tint in 1:length(LBMean)
        print(" & \$", round(UBMean[tint],1), " \\pm ", round(UBStdEv[tint],1), "\$");
    end
    println("\\\\");
end

# print out battery utilization results
databutil = Dict();
for i in 1:3
    databutil[i] = load("butilResults_$(i)_AG.jld");
end
for i in 1:3
    println("udet: ",sum(values(databutil[i]["detOut"][1])));
    println("ustoch: ",sum(values(databutil[i]["stochOut"][1])));
end


# print out the hardening results
datahard = Dict();
for i in 1:3
    datahard[i] = load("hardResults_$(i)_AG.jld");
end
for i in 1:3
    println("Null ",datadet[i]["stochOut"][24][1]);
    for ω in keys(datahard[i]["data"])
        println(ω,"    ",datahard[i]["data"][ω][1][length(datahard[i]["data"][ω][1])]);
    end
    println("+++++++++++++++++++++++++++");
end
