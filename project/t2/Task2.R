source("fcts.r")
library(EBImage)
library(spatstat)
library(bitops)


#### Additional Functions ####
# EDT
estAAbyEDT <- function(BW, m=0, spacing=1) {
  rr <- 0:m
  d <- dim(BW)
  BW.EDT <- EBImage::distmap(1-BW) # EDT
  AA <- sapply(rr, function(r) mean((BW.EDT[(1+r):(d[1]-r),(1+r):(d[2]-r)]<=r)@.Data))
  cbind(r=rr*spacing, AA.EDT=AA, Hs.EDT=(AA-AA[1])/(1-AA[1]))
}

# The following function 'getAA.rect.const' generates 'nrep' repetitions of a 2D Boolean model 
# with some intensity 'lambda' and typical grain equal to a constant rectangle of A x B 
# and each time estimates the first Minkowski function A_A(r).
getAA.rect.const <- function(nrep = 39, lambda, a, b, W = owin(c(0,1),c(0,1)), m=0, spacing=1) {
  sapply(1:nrep, function(k) {
    XYR <- rBM.rect.const(lambda=lambda, a=a,b = b, W=W)
    BW.XYR <- digitizeRectSys(XYR, spacing=spacing) 
    estAAbyEDT(BW=BW.XYR, m=m, spacing=spacing)[,2]
  })
}



#### Part 1####

#load image
img = readImage("image1_03.png")
display(img, method = "raster")



#Morphological opening

pxl = 0.025 #pixel length 

#loop through possible vals of pix size of b
for (px in c(1, 2.1, 2.5,3.3,4.1,4.3,4.7)){
  brush = makeBrush(size= px, shape='disc')
  i = EBImage::opening(img, brush )
  display(i,  method = "raster")
  writeImage(i,paste("img ", px*pxl,".png")) #save each img for each size of b
  print(px)
}

# 4.1 seems to be the right number of pixels
b = pxl *4.1
b
brush = makeBrush(size= b/pxl, shape='disc')
display(img, method = "raster")
display(EBImage::opening(img, brush ), method = "raster")

#Seems like value of b is between 10 and 11.
paste("Size of b is ", b)





####Part 2 ####
alpha = 0.05
m = 40
#generate models with the resized observation window

M1 = rBM.rect.const(lambda = 4.7,a = 1.71, b = 0.15)
img.M1 = digitizeRectSys(M1, spacing=0.025)
display(img.M1,method = "raster")

M2 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
img.M2 = digitizeRectSys( M2, spacing=0.025)
display(img.M2,method = "raster")

# Minkowski functions
M1.ALX = estALXFct(BW=img.M1, m=10, spacing=0.025, ms=T)
M1.ALX

M2.ALX = estALXFct(BW=img.M2, m=10, spacing=0.025, ms=T)
M2.ALX



####* Model 1 ####
M1 = rBM.rect.const(lambda = 4.7,a = 1.71, b = 0.15)
img.M1 = digitizeRectSys(M1, spacing=0.025)
display(img.M1,method = "raster")

# 999 realizations of model 1
M1.model.999 <- getAA.rect.const(nrep=999, lambda=4.7, a = 1.71,b=0.15, W = owin(c(0,10),c(0,10)), m=m, spacing=0.025)

# Estimates of the first Minkowski func
M1.data <- estAAbyEDT(BW=img.M1, m=m, spacing=0.025)[,1:2]

#global envelopes
M1.genv.999 <- globalEnvelopes(M1.model.999, alpha=alpha)
M1.genv.999

#plot
png(file = "Envelopes_model1.png", width = 480, height = 480)

plot(M1.data[,1], M1.data[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")), lwd=2, ylim=c(0.4,1))
lines(M1.data[,1], M1.genv.999[,1], lwd=2, col=2) # red
lines(M1.data[,1], M1.genv.999[,2], lwd=2, col=2) # red

dev.off()



#p values
p.values1 <- sapply(1:1000, function(k) {
  M1 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
  img.M1 = digitizeRectSys(M1, spacing=0.025)
  M1.data <- estAAbyEDT(BW=img.M1, m=m, spacing=0.025)[,1:2]
  globalTest(M1.data[,2], M1.model.999) 
})
p.values1 # shows all p-values
mean(p.values1<alpha)

# The p-value is smaller than alpha in 80.7% rep -> model is unsufficient
#=======================#

####* Model 2 ####
M2 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
img.M2 = digitizeRectSys(M2, spacing=0.025)
display(img.M2,method = "raster")

# 999 realizations of model 2
M2.model.999 <- getAA.rect.const(nrep=999, lambda=4.7, a = 1.71,b=0.15,
                                 W = owin(c(0,10),c(0,10)), m=m, spacing=0.025)

# Estimates of the first Minkowski func
M2.data <- estAAbyEDT(BW=img.M2, m=m, spacing=0.025)[,1:2]

#global envelopes
M2.genv.999 <- globalEnvelopes(M2.model.999, alpha=alpha)
M2.genv.999

#plot
png(file = "Envelopes_model2.png", width = 480, height = 480)

plot(M2.data[,1], M2.data[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")),
     lwd=2, ylim=c(0.4,1),col=4)
lines(M2.data[,1], M2.genv.999[,1], lwd=2, col=2) # red
lines(M2.data[,1], M2.genv.999[,2], lwd=2, col=2) # red


dev.off()


#global test
p.values2 <- sapply(1:1000, function(k) {
  M2 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
  img.M2 = digitizeRectSys(M2, spacing=0.025)
  M2.data <- estAAbyEDT(BW=img.M2, m=m, spacing=0.025)[,1:2]
  globalTest(M2.data[,2], M2.model.999) 
})
p.values2 # shows all p-values
mean(p.values2<alpha)

# The p-value is smaller than alpha  rep -> model is unsufficient





########Try with estALXfct #############

####* mod 1####
M1 = rBM.rect.const(lambda = 4.7,a = 1.71, b = 0.15)
img.M1 = digitizeRectSys(M1, spacing=0.025)
display(img.M1,method = "raster")


# 999 realizations of model 1
M1.model.999 <- getAA.rect.const(nrep=999, lambda=4.7, a = 1.71,b=0.15, W = owin(c(0,10),c(0,10)), m=m, spacing=0.025)

# Estimates of the first Minkowski func
M1.data <-estALXFct(BW=img.M1, m=m, spacing=0.025)

#global envelopes
M1.genv.999 <- globalEnvelopes(M1.model.999, alpha=alpha)
M1.genv.999

#plot
png(file = "Envelopes_model1_estalx.png", width = 480, height = 480)

plot(M1.data[,1], M1.data[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")), lwd=2, ylim=c(0.4,1))
lines(M1.data[,1], M1.genv.999[,1], lwd=2, col=2) # red
lines(M1.data[,1], M1.genv.999[,2], lwd=2, col=2) # red

dev.off()



#p values
p.values1 <- sapply(1:1000, function(k) {
  M1 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
  img.M1 = digitizeRectSys(M1, spacing=0.025)
  M1.data <- estALXFct(BW=img.M1, m=m, spacing=0.025)[,1:2]
  globalTest(M1.data[,2], M1.model.999) 
})
p.values1 # shows all p-values
mean(p.values1<alpha)

# The p-value is smaller than alpha in 94% rep -> model is unsufficient
#=======================#

####* mod 2 ####
M2 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
img.M2 = digitizeRectSys(M2, spacing=0.025)
display(img.M2,method = "raster")

# 999 realizations of model 2
M2.model.999 <- getAA.rect.const(nrep=999, lambda=4.7, a = 1.71,b=0.15,
                                 W = owin(c(0,10),c(0,10)), m=m, spacing=0.025)

# Estimates of the first Minkowski func
M2.data <- estALXFct(BW=img.M2, m=m, spacing=0.025)

#global envelopes
M2.genv.999 <- globalEnvelopes(M2.model.999, alpha=alpha)
M2.genv.999

#plot
png(file = "Envelopes_model2_estalx.png", width = 480, height = 480)

plot(M2.data[,1], M2.data[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")),
     lwd=2, ylim=c(0.4,1),col=4)
lines(M2.data[,1], M2.genv.999[,1], lwd=2, col=2) # red
lines(M2.data[,1], M2.genv.999[,2], lwd=2, col=2) # red


dev.off()


#global test
p.values2 <- sapply(1:1000, function(k) {
  M2 = rBM.rect.const(lambda = 5.1,a = 1.98, b = 0.11)
  img.M2 = digitizeRectSys(M2, spacing=0.025)
  M2.data <- estALXFct(BW=img.M2, m=m, spacing=0.025)[,1:2]
  globalTest(M2.data[,2], M2.model.999) 
})
p.values2 # shows all p-values
mean(p.values2>alpha)
# The p-value is smaller than alpha in 94% rep -> model is unsufficient
pvals_total = c(0.775,0.785,0.961,0.934)
