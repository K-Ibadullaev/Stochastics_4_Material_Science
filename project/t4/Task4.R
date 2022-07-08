source("fcts.r")
library(EBImage)
library(spatstat)
library(bitops)

######Generate models ######
m = 40

#CHANGE TO 100 !!
###lambda = 5,R = 0.05####
simL5<- sapply(1:100, function(k) {
  
  
  L5 = rM3.disc.const(lambda = 5,R = 0.05, W = owin(c(0,10),c(0,10)) )
  BW.L5 = digitizeDiscSys(L5,spacing = 0.01)
  display(BW.L5, method = "raster")
  ALX.L5 = estALX(BW.L5 ,spacing=0.01)
})

simL5

# mean and sd for A 
AmeanL5 = mean(simL5[1,])
AmeanL5 



AsdL5 = sd(simL5[1,])
AsdL5

# mean and sd for L 
LmeanL5 = mean(simL5[2,])
LmeanL5 

LsdL5 = sd(simL5[2,])
LsdL5

# mean and sd for Xi 
XimeanL5 = mean(simL5[3,])
XimeanL5 

XisdL5 = sd(simL5[3,])
XisdL5


#####lambda = 10,R = 0.05,#####
simL10<- sapply(1:100, function(k) {
  L10 = rM3.disc.const(lambda = 10,R = 0.05, W = owin(c(0,10),c(0,10)) )
  BW.L10 = digitizeDiscSys(L10,spacing = 0.01)
  display(BW.L10, method = "raster")
  ALX.L10 = estALX(BW.L10 ,spacing=0.01)
  
 
})

simL10

# mean and sd for A 
AmeanL10 = mean(simL10[1,])
AmeanL10

AsdL10 = sd(simL10[1,])
AsdL10

# mean and sd for L 
LmeanL10 = mean(simL10[2,])
LmeanL10 

LsdL10 = sd(simL10[2,])
LsdL10

# mean and sd for Xi 
XimeanL10 = mean(simL10[3,])
XimeanL10 

XisdL10 = sd(simL10[3,])
XisdL10




###### lambda = 15,R = 0.05,#######

simL15<- sapply(1:100, function(k) {
  
  L15 = rM3.disc.const(lambda = 15,R = 0.05, W = owin(c(0,10),c(0,10)) )
  BW.L15 = digitizeDiscSys(L15,spacing = 0.01)
  display(BW.L15, method = "raster")
  ALX.L15 = estALX(BW.L15 ,spacing=0.01)
  
  
})

simL15

# mean and sd for A 
AmeanL15 = mean(simL15[1,]) #monte carlo estimate????
AmeanL15

AsdL15 = sd(simL15[1,])
AsdL15

#boxplot of A
boxplot(simL15[1], ylab=expression(A[A]))

# mean and sd for L 
LmeanL15 = mean(simL15[2,])
LmeanL15 

LsdL15 = sd(simL15[2,])
LsdL15

# mean and sd for Xi 
XimeanL15= mean(simL15[3,])
XimeanL15 

XisdL15= sd(simL15[3,])
XisdL15


#boxplots

#for A
png(file = "boxplots_A.png", width = 480, height = 480)

boxplot(simL5[1,],simL10[1,],simL15[1,], ylab=expression(A[A]), 
        names = c("lambda 5","lambda 10","lambda 15"),main= expression(A[A]))

abline(h=c(AmeanL5,AmeanL10,AmeanL15),col=c(2,3,4)) # mean value
dev.off()

#for L
png(file = "boxplots_L.png", width = 480, height = 480)
boxplot(simL5[2,],simL10[2,],simL15[2,], ylab=expression(L[A]), 
        names = c("lambda 5","lambda 10","lambda 15"),main= expression(L[A]))
abline(h=c(LmeanL5,LmeanL10,LmeanL15),col=c(2,3,4))# mean value
dev.off()

#for X
png(file = "boxplots_X.png", width = 480, height = 480)
boxplot(simL5[3,],simL10[3,],simL15[3,], ylab=expression(X[A]), 
        names = c("lambda 5","lambda 10","lambda 15"),main= expression(X[A]))

abline(h=c(XimeanL5,XimeanL10,XimeanL15),col=c(2,3,4)) # mean value
dev.off()

#total params

total = data.frame(Model=c("lambda 5","lambda 10","lambda 15"),
                   A_mean = c(AmeanL5,AmeanL10,AmeanL15),A_std = c(AsdL5,AsdL10,AsdL15),
                   L_mean = c(LmeanL5,LmeanL10,LmeanL15),L_std = c(LsdL5,LsdL10,LsdL15),
                   X_mean = c(XimeanL5,XimeanL10,XimeanL15),X_std = c(XisdL5,XisdL10,XisdL15)
                   )
total
write.table(x=total, file = "total_tsk4.txt",sep="| ")

