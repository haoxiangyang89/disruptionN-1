# output csv files
@everywhere include("loadMod.jl");
caseList = [13,33,123];

# output the csv file for GenAll and dOnly
for ci in 1:length(caseList)
    datad = load("dOnlyResults_$(ci).jld");
    lbOutList = [datad["data"]["dOnly"][1][1:40] datad["data"]["allGen"][1][1:40]];
    timeOutList = zeros(41,2);
    cT1 = 0;
    cT2 = 0;
    for j in 1:40
        cT1 += datad["data"]["dOnly"][5][j];
        cT2 += datad["data"]["allGen"][5][j];
        timeOutList[j+1,1] = cT1;
        timeOutList[j+1,2] = cT2;
    end
    writedlm("./csvOut/dOnly_lb_$(ci).csv",lbOutList,',');
    writedlm("./csvOut/dOnly_time_$(ci).csv",timeOutList,',');
end

# output the csv file for preGen
for ci in 1:length(caseList)
    datapg = load("pgResults_$(ci).jld");
    lbOutList = [datapg["data"]["dOnly"][1][1:40] datapg["data"]["preGen"][1][1:40]];
    timeOutList = zeros(41,2);
    cT1 = 0;
    cT2 = datapg["data"]["preGen"][5][7];
    timeOutList[1,2] = cT2;
    for j in 1:40
        cT1 += datapg["data"]["dOnly"][5][j];
        cT2 += datapg["data"]["preGen"][5][j];
        timeOutList[j+1,1] = cT1;
        timeOutList[j+1,2] = cT2;
    end
    writedlm("./csvOut/pg_lb_$(ci).csv",lbOutList,',');
    writedlm("./csvOut/pg_time_$(ci).csv",timeOutList,',');
end

# output the csv file for NTest
NList = [1,5,10,20,30];
NNoList = [200,40,20,20,20];
lbDict = Dict();
timeDict = Dict();
dataN = Dict();
for ci in 1:length(caseList)
    dataN[ci] = load("NResults_$(ci).jld");
end
for Nind in 1:length(NList)
    N = NList[Nind];
    lbDict[N] = [];
    timeDict[N] = [];
    NNo = NNoList[Nind];
    for ci in 1:length(caseList)
        if ci != 1
            lbDict[N] = [lbDict[N] dataN[ci]["NOut"][N][1][1:NNo]];
            tdTemp = zeros(NNo+1);
            CT = 0;
            for j in 1:NNo
                CT += dataN[ci]["NOut"][N][5][j];
                tdTemp[j+1] = CT;
            end
            timeDict[N] = [timeDict[N] tdTemp];
        else
            lbDict[N] = dataN[ci]["NOut"][N][1][1:NNo];
            tdTemp = zeros(NNo+1);
            CT = 0;
            for j in 1:NNo
                CT += dataN[ci]["NOut"][N][5][j];
                tdTemp[j+1] = CT;
            end
            timeDict[N] = tdTemp;
        end
    end
    writedlm("./csvOut/N_lb_$(Nind).csv",lbDict[N],',');
    writedlm("./csvOut/N_time_$(Nind).csv",timeDict[N],',');
end
