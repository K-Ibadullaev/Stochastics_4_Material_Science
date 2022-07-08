source("fcts.r")
library(EBImage)
library(spatstat)
library(bitops)

#Load images
img1 = readImage('image1_03.png')
display(img1,method="raster"  )

img2 = readImage('image2_03.png')
display(img2,method="raster"  )

 img1.ALX = estALX(img1, spacing=0.025)
 img1.ALX
 
 img2.ALX = estALX(img2, spacing=0.025)
 img2.ALX

 delta.ALX = data.frame(img2.ALX) - data.frame(img1.ALX)
 delta.ALX
#Estimate Minkowski functions
img1.ALX = estALXFct(BW=img1,m=10, spacing=0.025, ms=FALSE)
img1.ALX

img2.ALX = estALXFct(BW=img2,m=10, spacing=0.025, ms=FALSE)
img2.ALX

#Plot functions
png(file = "plot1_img1.png", width = 480, height = 480)
plotALXFct(BW=img1, ALX=img1.ALX, show.type=0)
dev.off()


png(file = "plot1_img2.png", width = 480, height = 480)
plotALXFct(BW=img2, ALX=img2.ALX, show.type=0)
dev.off()

#Plot jointly each  Minkowski function,resp
png(file = "1st_Minkowski.png", width = 480, height = 480)
plotALXFct(BW=img2, ALX=img1.ALX, ALX2 =img2.ALX ,show.type=1) 
dev.off()

png(file = "2nd_Minkowski.png", width = 480, height = 480)
plotALXFct(BW=img2, ALX=img1.ALX, ALX2 =img2.ALX ,show.type=2)
dev.off()

png(file = "3rd_Minkowski.png", width = 480, height = 480)
plotALXFct(BW=img2, ALX=img1.ALX, ALX2 =img2.ALX ,show.type=3)
dev.off()
