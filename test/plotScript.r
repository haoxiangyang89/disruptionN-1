LBList <- read.csv("./Desktop/Git/disruptionN-1/test/LBList.csv")
png(file = "./Desktop/Git/disruptionN-1/test/LBList.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(LBList$LB96,type = "l",ylim=range(c(0,14000)),xlab = "Iterations",ylab = "Obj Value", 
     main = "Lower Bound", col = "#377EB8", lwd = 5, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(LBList$LB72, col = "#E41A1C", lwd = 5)
lines(LBList$LB48, col = "#4DAF4A", lwd = 5)
lines(LBList$LB24, col = "#984EA3", lwd = 5)
lines(LBList$LB12, col = "#FF7F00", lwd = 5)
legend("topleft",c("T = 96","T = 72","T = 48","T = 24","T = 12"), col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5)
dev.off();

tauList <- read.csv("./Desktop/Git/disruptionN-1/test/tauList.csv")
png(file = "./Desktop/Git/disruptionN-1/test/tauList.png", width= 10,height=6,units = 'in',res = 300);
par(mar = c(5,5,2,2));
plot(tauList$tau,tauList$LB,type = "l",ylim=range(c(1580,2800)),xlab = "Recovery Length",ylab = "Obj Value", 
     main = "Cost vs. Recovery Length", col = "#377EB8", lwd = 5, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(tauList$tau, tauList$X1000_UB, col = "#E41A1C", lwd = 5)
legend("topleft",c("Lower bound","Simulated Estimation"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5)

dev.off();

#=============================================================================================
# plot 1: GenAll vs. DOnly
dOnlylb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_1.csv", header = FALSE)
dOnlylb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_2.csv", header = FALSE)
dOnlylb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_lb_3.csv", header = FALSE)

dOnlyt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_1.csv", header = FALSE)
dOnlyt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_2.csv", header = FALSE)
dOnlyt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/dOnly_time_3.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/GenAll.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:20, dOnlylb1$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,40000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:100, dOnlylb1$V2[1:100], col = "#E41A1C", lwd = 3, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlylb2$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,38000)), xlab = "Iteration", ylab = "LB", main = "Case 33",
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:100, dOnlylb2$V2[1:100], col = "#E41A1C", lwd = 3, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlylb3$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,50000)), xlab = "Iteration", ylab = "LB", main = "Case 123",
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:100, dOnlylb3$V2[1:100], col = "#E41A1C", lwd = 3, log="x");
legend("bottomright",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(1:20, dOnlyt1$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,15000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:100, dOnlyt1$V2[1:100], col = "#E41A1C", lwd = 3, log="x");
legend("topleft",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlyt2$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,20000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:100, dOnlyt2$V2[1:100], col = "#E41A1C", lwd = 3, log="x");
legend("topleft",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:20, dOnlyt3$V1[1:20], type = "l", xlim = range(c(1,100)), ylim=range(c(0,40000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:100, dOnlyt3$V2[1:100], col = "#E41A1C", lwd = 3, log="x");
legend("topleft",c("GenAll","DOnly"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 2: preGen
pglb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/pg_lb_1.csv", header = FALSE)
pglb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/pg_lb_2.csv", header = FALSE)
pglb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/pg_lb_3.csv", header = FALSE)

pgt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/pg_time_1.csv", header = FALSE)
pgt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/pg_time_2.csv", header = FALSE)
pgt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/pg_time_3.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/pgFig.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:60, pglb1$V1[1:60], type = "l", ylim=range(c(0,18000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:60, pglb1$V2[1:60], col = "#E41A1C", lwd = 3);
legend("bottomright",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:60, pglb2$V1[1:60], type = "l", ylim=range(c(0,3600)), xlab = "Iteration", ylab = "LB", main = "Case 33",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:60, pglb2$V2[1:60], col = "#E41A1C", lwd = 3);
legend("bottomright",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:60, pglb3$V1[1:60], type = "l", ylim=range(c(0,48000)), xlab = "Iteration", ylab = "LB", main = "Case 123",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:60, pglb3$V2[1:60], col = "#E41A1C", lwd = 3);
legend("bottomright",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(pgt1$V1[1:60], type = "l", ylim=range(c(0,2500)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(pgt1$V2[1:60], col = "#E41A1C", lwd = 3);
legend("topleft",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(pgt2$V1[1:60], type = "l", ylim=range(c(0,4500)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(pgt2$V2[1:60], col = "#E41A1C", lwd = 3);
legend("topleft",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(pgt3$V1[1:60], type = "l", ylim=range(c(0,9000)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(pgt3$V2[1:60], col = "#E41A1C", lwd = 3);
legend("topleft",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 3: NTest
Nlb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_1.csv", header = FALSE)
Nlb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_2.csv", header = FALSE)
Nlb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_3.csv", header = FALSE)
Nlb4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_4.csv", header = FALSE)
Nlb5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_5.csv", header = FALSE)
Nlb6 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_6.csv", header = FALSE)
Nlb7 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_7.csv", header = FALSE)

Nt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_1.csv", header = FALSE)
Nt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_2.csv", header = FALSE)
Nt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_3.csv", header = FALSE)
Nt4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_4.csv", header = FALSE)
Nt5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_5.csv", header = FALSE)
Nt6 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_6.csv", header = FALSE)
Nt7 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_7.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/NFig.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:60, Nlb1$V1[1:60], type = "l", ylim=range(c(0,40000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
#lines(1:50, Nlb2$V1, col = "#000000", lwd = 3, log="x");
lines(1:20, Nlb3$V1[1:20], col = "#E41A1C", lwd = 3, log="x");
#lines(1:25, Nlb4$V1, col = "#FF7F00", lwd = 3, log="x");
lines(1:12, Nlb5$V1[1:12], col = "#4DAF4A", lwd = 3, log="x");
lines(1:6, Nlb6$V1[1:6], col = "#984EA3", lwd = 3, log="x");
lines(1:3, Nlb7$V1[1:3], col = "#3CAEA3", lwd = 3, log="x");
#legend("bottomright",c("N = 1","N = 2","N = 3","N = 4","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#000000","#3CAEA3"),pch = 20,cex = 1.5);
legend("bottomright",c("N = 1","N = 3","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(1:60, Nlb1$V2[1:60], type = "l", ylim=range(c(0,38000)), xlab = "Iteration", ylab = "LB", main = "Case 33",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
#lines(1:50, Nlb2$V2, col = "#000000", lwd = 3, log="x");
lines(1:20, Nlb3$V2[1:20], col = "#E41A1C", lwd = 3, log="x");
#lines(1:25, Nlb4$V2, col = "#FF7F00", lwd = 3, log="x");
lines(1:12, Nlb5$V2[1:12], col = "#4DAF4A", lwd = 3, log="x");
lines(1:6, Nlb6$V2[1:6], col = "#984EA3", lwd = 3, log="x");
lines(1:3, Nlb7$V2[1:3], col = "#3CAEA3", lwd = 3, log="x");
#legend("bottomright",c("N = 1","N = 2","N = 3","N = 4","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#000000","#3CAEA3"),pch = 20,cex = 1.5);
legend("bottomright",c("N = 1","N = 3","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(1:60, Nlb1$V3[1:60], type = "l", ylim=range(c(0,50000)), xlab = "Iteration", ylab = "LB", main = "Case 123",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
#lines(1:50, Nlb2$V3, col = "#000000", lwd = 3, log="x");
lines(1:20, Nlb3$V3[1:20], col = "#E41A1C", lwd = 3, log="x");
#lines(1:25, Nlb4$V3, col = "#FF7F00", lwd = 3, log="x");
lines(1:12, Nlb5$V3[1:12], col = "#4DAF4A", lwd = 3, log="x");
lines(1:6, Nlb6$V3[1:6], col = "#984EA3", lwd = 3, log="x");
lines(1:3, Nlb7$V3[1:3], col = "#3CAEA3", lwd = 3, log="x");
#legend("bottomright",c("N = 1","N = 2","N = 3","N = 4","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#000000","#3CAEA3"),pch = 20,cex = 1.5);
legend("bottomright",c("N = 1","N = 3","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:60, Nt1$V1[1:61], type = "l", ylim=range(c(0,20000)), xlab = "Iteration", ylab = "Time (sec.)",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
#lines(1:50, Nt2$V1, col = "#000000", lwd = 3, log="x");
lines(0:20, Nt3$V1[1:21], col = "#E41A1C", lwd = 3, log="x");
#lines(1:25, Nt4$V1, col = "#FF7F00", lwd = 3, log="x");
lines(0:12, Nt5$V1[1:13], col = "#4DAF4A", lwd = 3, log="x");
lines(0:6, Nt6$V1[1:7], col = "#984EA3", lwd = 3, log="x");
lines(0:3, Nt7$V1[1:4], col = "#3CAEA3", lwd = 3, log="x");
#legend("topleft",c("N = 1","N = 2","N = 3","N = 4","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#000000","#3CAEA3"),pch = 20,cex = 1.5);
legend("topleft",c("N = 1","N = 3","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);

plot(0:60, Nt1$V2[1:61], type = "l", ylim=range(c(0,35000)), xlab = "Iteration", ylab = "Time (sec.)",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
#lines(1:50, Nt2$V2, col = "#000000", lwd = 3, log="x");
lines(0:20, Nt3$V2[1:21], col = "#E41A1C", lwd = 3, log="x");
#lines(1:25, Nt4$V2, col = "#FF7F00", lwd = 3, log="x");
lines(0:12, Nt5$V2[1:13], col = "#4DAF4A", lwd = 3, log="x");
lines(0:6, Nt6$V2[1:7], col = "#984EA3", lwd = 3, log="x");
lines(0:3, Nt7$V2[1:4], col = "#3CAEA3", lwd = 3, log="x");
#legend("topleft",c("N = 1","N = 2","N = 3","N = 4","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#000000","#3CAEA3"),pch = 20,cex = 1.5);
legend("topleft",c("N = 1","N = 3","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);

plot(0:60, Nt1$V3[1:61], type = "l", ylim=range(c(0,60000)), xlab = "Iteration", ylab = "Time (sec.)",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
#lines(1:50, Nt2$V3, col = "#000000", lwd = 3, log="x");
lines(0:20, Nt3$V3[1:21], col = "#E41A1C", lwd = 3, log="x");
#lines(1:25, Nt4$V3, col = "#FF7F00", lwd = 3, log="x");
lines(0:12, Nt5$V3[1:13], col = "#4DAF4A", lwd = 3, log="x");
lines(0:6, Nt6$V3[1:7], col = "#984EA3", lwd = 3, log="x");
lines(0:3, Nt7$V3[1:4], col = "#3CAEA3", lwd = 3, log="x");
#legend("topleft",c("N = 1","N = 2","N = 3","N = 4","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00","#000000","#3CAEA3"),pch = 20,cex = 1.5);
legend("topleft",c("N = 1","N = 3","N = 5","N = 10","N = 20"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#3CAEA3"),pch = 20,cex = 1.5);
dev.off()