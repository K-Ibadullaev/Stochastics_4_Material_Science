source("fcts.r")
library(EBImage)
library(spatstat)
library(bitops)

img = readImage("image3_03.png")
display(img, method = "raster")

#### Estimates of Vv and Sv ####

###### Estimate Vv under an assumption of the isotropy ######

img.V = estALX(B = img, spacing = 0.001)[1]
img.V

# The 2nd possibility by the point count method
img.Vfrac.mean = mean(img)
img.Vfrac.mean 
# both estimates of the volume fraction differ slightly
img.Vfrac.mean / img.V

###### Estimate Sv under an assumption of the isotropy ######
img.S = 4*estALX(B = img, spacing = 0.001)[2]/pi
img.S


#### Find diameters and implement segmentation algorithm  ####

# label grains 
img.labelled = bwlabel(img)

# visual check 
cols = c('black', sample(rainbow(max(img.labelled))))
img.labelled.coloured = Image(cols[1+img.labelled], dim=dim(img.labelled))

png(file = "labelled.png", width = 480, height = 480)
display(img.labelled.coloured, method="raster") 
dev.off()


# get features
img.feat <- computeFeatures.moment(img.labelled)

# check names
colnames(img.feat)

# We need the diameter = 'm.majoraxis' and 
# coordinates of their centres (x,y)- "m.cx" and "m.cy"

d = img.feat[,3] #discs diameters
length(d)

centers = img.feat[,1:2] # coordinates of discs centres
nrow(centers)


# Applying an idea of the minus-sampling

centers = img.feat[img.feat[,1]>30 & img.feat[,1]<1371 & img.feat[,2]>30 & img.feat[,2]<1371, 1:2]
nrow(centers) # number of discs with centres inside W

# Discs inside window 1340 x 1340
d.2D <- img.feat[img.feat[,1]>30 & img.feat[,1]<1371 & img.feat[,2]>30 & img.feat[,2]<1371, 3]
head(d.2D)# unit is 1 pixel = 1 micrometre
d.2D <- d.2D/1000 # We switch to unit 1 millimetre.




# number of discs with centre in the reduced window:
n.2D <- length(d.2D)
n.2D



# area of reduced window in square millimeters:
A.W <- prod(dim(img)-2*30)/1000^2 
A.W

# mean and standard deviation
mean(d.2D)
sd(d.2D)

# estimated intensity in 2D:
lambdaA <- n.2D/A.W # unit is number per square millimetre
lambdaA

# range of disc diameters:
range = range(d.2D) # in millimeters

# plot histogram
Delta <- 0.02

png(file = "histogram.png", width = 480, height = 480)
fi.2D <- hist(d.2D, breaks=seq(0, 0.12, by=Delta))$counts
dev.off()



# absolute frequencies of disc diameters in 2D per per square millimetre:
ni.2D <- fi.2D/A.W
ni.2D

total2d = data.frame(Discs_numbers = n.2D,
                     Area_W =A.W,
                     mean_d =mean(d.2D),standard_deviation =sd(d.2D),
                     lambdaA = lambdaA,
                     range = paste(range(d.2D)[1],range(d.2D)[2],sep = "-")
                     )
total2d
write.table(x=total2d, file = "total2d_tsk3.txt",sep="| ")
#### Implementation of the Saltykov's algorithm  ####

Ni.3D <- saltykov(Delta=Delta, ni.2D, nEM=32)
# -> 'f.3D' contains absolute frequencies of ball diameters in 3D per cubic millimetre:
Ni.3D

# comparison of relative frequencies:
png(file = "barplot_3Dvs2D.png", width = 480, height = 480)
barplot(rbind(ni.2D/sum(ni.2D),Ni.3D/sum(Ni.3D)), beside=TRUE, 
        legend.text = c("2D", "3D"), args.legend = list(x = "topleft"),
        names.arg=c("(0,0.02]","(0.02,0.04]","(0.04,0.06]","(0.06,0.08]","(0.08,0.1]","(0.1,0.12]"))
dev.off()


#### Calculate lambda and mean ball diameter by : ####

# estimated intensity and mean ball diameter  in 3D, with formula (6) of the lecture 10
lambdaV.1 <- sum(1/d.2D)*2/pi/A.W 
lambdaV.1

muV.1 <- n.2D*pi/2/sum(1/d.2D)
muV.1

# estimated intensity and mean ball diameter  in 3D, with formula (7) of the lecture 10

lambdaV.2 <- sum(Ni.3D) # estimated intensity in 3D, with formula (7)
lambdaV.2

muV.2 <- n.2D/(lambdaV.2 * A.W) 
muV.2

total3d = data.frame(
                     lambdaV.1=lambdaV.1,
                     muV.1=muV.1,
                     lambdaV.1=lambdaV.1,
                     muV.1=muV.1)
total3d
write.table(x=total3d, file = "total3d_tsk3.txt",sep="| ",)

freq = data.frame(Absolute_freq = Ni.3D)
write.table(x=freq, file = "absfreq3d_tsk3.txt",sep="| ",)
