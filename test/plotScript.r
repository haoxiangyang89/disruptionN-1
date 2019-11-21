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