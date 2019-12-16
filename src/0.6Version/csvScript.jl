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
