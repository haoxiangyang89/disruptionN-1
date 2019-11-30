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
plot(1:40, dOnlylb1$V1, type = "l", ylim=range(c(0,25)), xlab = "Iteration", ylab = "LB", main = "Case 13",
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:40, dOnlylb1$V2, col = "#E41A1C", lwd = 3);
legend("bottomright",c("DOnly","GenAll"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:40, dOnlylb2$V1, type = "l", ylim=range(c(0,25)), xlab = "Iteration", ylab = "LB", main = "Case 33",
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:40, dOnlylb2$V2, col = "#E41A1C", lwd = 3);
legend("bottomright",c("DOnly","GenAll"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:40, dOnlylb3$V1, type = "l", ylim=range(c(0,15)), xlab = "Iteration", ylab = "LB", main = "Case 123",
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:40, dOnlylb3$V2, col = "#E41A1C", lwd = 3);
legend("bottomright",c("DOnly","GenAll"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(dOnlyt1$V1, type = "l", ylim=range(c(0,3100)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(dOnlyt1$V2, col = "#E41A1C", lwd = 3);
legend("topleft",c("DOnly","GenAll"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt2$V1, type = "l", ylim=range(c(0,3200)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(dOnlyt2$V2, col = "#E41A1C", lwd = 3);
legend("topleft",c("DOnly","GenAll"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(dOnlyt3$V1, type = "l", ylim=range(c(0,7000)), xlab = "Iteration", ylab = "Time (sec.)", 
      col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(dOnlyt3$V2, col = "#E41A1C", lwd = 3);
legend("topleft",c("DOnly","GenAll"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
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
plot(1:40, pglb1$V1, type = "l", ylim=range(c(0,25)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:40, pglb1$V2, col = "#E41A1C", lwd = 3);
legend("bottomright",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:40, pglb2$V1, type = "l", ylim=range(c(0,25)), xlab = "Iteration", ylab = "LB", main = "Case 33",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:40, pglb2$V2, col = "#E41A1C", lwd = 3);
legend("bottomright",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(1:40, pglb3$V1, type = "l", ylim=range(c(0,15)), xlab = "Iteration", ylab = "LB", main = "Case 123",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(1:40, pglb3$V2, col = "#E41A1C", lwd = 3);
legend("bottomright",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(pgt1$V1, type = "l", ylim=range(c(0,500)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(pgt1$V2, col = "#E41A1C", lwd = 3);
legend("topleft",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(pgt2$V1, type = "l", ylim=range(c(0,700)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(pgt2$V2, col = "#E41A1C", lwd = 3);
legend("topleft",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
plot(pgt3$V1, type = "l", ylim=range(c(0,2000)), xlab = "Iteration", ylab = "Time (sec.)", 
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5);
lines(pgt3$V2, col = "#E41A1C", lwd = 3);
legend("topleft",c("No pre-generated cuts","Pre-generated cuts"), col = c("#377EB8","#E41A1C"),pch = 20,cex = 1.5);
dev.off();

# plot 3: NTest
Nlb1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_1.csv", header = FALSE)
Nlb2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_2.csv", header = FALSE)
Nlb3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_3.csv", header = FALSE)
Nlb4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_4.csv", header = FALSE)
Nlb5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_lb_5.csv", header = FALSE)

Nt1 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_1.csv", header = FALSE)
Nt2 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_2.csv", header = FALSE)
Nt3 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_3.csv", header = FALSE)
Nt4 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_4.csv", header = FALSE)
Nt5 <- read.csv("./Desktop/Git/disruptionN-1/test/csvOut/N_time_5.csv", header = FALSE)

outString = "./Desktop/Git/disruptionN-1/test/csvOut/NFig.png"
png(file = outString, width= 13, height = 8, units = 'in',res = 300);
par(mfrow=c(2,3));
par(mar = c(5,5,2.5,2.5));
plot(1:200, Nlb1$V1, type = "l", ylim=range(c(0,40)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:40, Nlb2$V1, col = "#E41A1C", lwd = 3, log="x");
lines(1:20, Nlb3$V1, col = "#4DAF4A", lwd = 3, log="x");
lines(1:20, Nlb4$V1, col = "#984EA3", lwd = 3, log="x");
lines(1:20, Nlb5$V1, col = "#FF7F00", lwd = 3, log="x");
legend("topleft",c("N = 1","N = 5","N = 10","N = 20","N = 30"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);

plot(1:200, Nlb1$V2, type = "l", ylim=range(c(0,60)), xlab = "Iteration", ylab = "LB", main = "Case 33",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:40, Nlb2$V2, col = "#E41A1C", lwd = 3, log="x");
lines(1:20, Nlb3$V2, col = "#4DAF4A", lwd = 3, log="x");
lines(1:20, Nlb4$V2, col = "#984EA3", lwd = 3, log="x");
lines(1:20, Nlb5$V2, col = "#FF7F00", lwd = 3, log="x");
legend("topleft",c("N = 1","N = 5","N = 10","N = 20","N = 30"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);

plot(1:200, Nlb1$V3, type = "l", ylim=range(c(0,40)), xlab = "Iteration", ylab = "LB", main = "Case 123",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(1:40, Nlb2$V3, col = "#E41A1C", lwd = 3, log="x");
lines(1:20, Nlb3$V3, col = "#4DAF4A", lwd = 3, log="x");
lines(1:20, Nlb4$V3, col = "#984EA3", lwd = 3, log="x");
lines(1:20, Nlb5$V3, col = "#FF7F00", lwd = 3, log="x");
legend("topleft",c("N = 1","N = 5","N = 10","N = 20","N = 30"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);

par(mar = c(5,5,1.5,2.5));
plot(0:200, Nt1$V1, type = "l", ylim=range(c(0,8000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(0:40, Nt2$V1, col = "#E41A1C", lwd = 3, log="x");
lines(0:20, Nt3$V1, col = "#4DAF4A", lwd = 3, log="x");
lines(0:20, Nt4$V1, col = "#984EA3", lwd = 3, log="x");
lines(0:20, Nt5$V1, col = "#FF7F00", lwd = 3, log="x");
legend("topleft",c("N = 1","N = 5","N = 10","N = 20","N = 30"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);

plot(0:200, Nt1$V2, type = "l", ylim=range(c(0,13000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(0:40, Nt2$V2, col = "#E41A1C", lwd = 3, log="x");
lines(0:20, Nt3$V2, col = "#4DAF4A", lwd = 3, log="x");
lines(0:20, Nt4$V2, col = "#984EA3", lwd = 3, log="x");
lines(0:20, Nt5$V2, col = "#FF7F00", lwd = 3, log="x");
legend("topleft",c("N = 1","N = 5","N = 10","N = 20","N = 30"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);

plot(0:200, Nt1$V3, type = "l", ylim=range(c(0,40000)), xlab = "Iteration", ylab = "LB", main = "Case 13",
     col = "#377EB8", lwd = 3, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, log="x");
lines(0:40, Nt2$V3, col = "#E41A1C", lwd = 3, log="x");
lines(0:20, Nt3$V3, col = "#4DAF4A", lwd = 3, log="x");
lines(0:20, Nt4$V3, col = "#984EA3", lwd = 3, log="x");
lines(0:20, Nt5$V3, col = "#FF7F00", lwd = 3, log="x");
legend("topleft",c("N = 1","N = 5","N = 10","N = 20","N = 30"), , col = c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00"),pch = 20,cex = 1.5);
dev.off();