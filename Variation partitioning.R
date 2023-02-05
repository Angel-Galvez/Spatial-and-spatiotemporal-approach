#-----------------------------------------------------------------------------------------------
# Script: gam for multivariate data. Spatial and spatiotemporal approaches: example
# Project: METACOM 1-3
# Date: 05/02/2023
#-----------------------------------------------------------------------------------------------
library(vegan)
library(mgcv)
library(gllvm)
library(doParallel)
library(foreach)
library(spdep)
library(adespatial)
library(dplyr)

setwd("/home/lv31/lv31033")
source("Gam_functions.R") 

# load data
sample <- read.table("samplev.csv", header = TRUE, row.names = 1, sep = ";") #Matrix with sampling campaign and site name for each row
xy <- read.csv2("XYv.csv", header = TRUE, row.names = 1, sep = ";") #XYDays matrix
env <- read.csv2("envbuenov.csv", header = TRUE, row.names = 1, sep = ";") #environmental matrix
sp <- read.csv2("amph_speciesv.csv", header = TRUE, row.names = 1, sep = ";") #Species matrix (Presence/Absence)
#------------------------------------------------------------------------------------------------
# Arrange data
sp <- sp[,colSums(sp) != 0]
gll.hat <- gll.hat.fun(sp = sp, family = "binomial", max.lv = 3) # latent variables

names<-select(sample, site, campaign)
time<-data.frame(xy$Days)
xy2<-xy
xy2$Days<-NULL
coords<-xy2

# PCA environment
pca.env <- as.data.frame(scores(rda(env[,c(1:6)], scale = TRUE), choices = c(1:6), display = "si")) #First 6 PCs

#-------------------------------------------------------------------------------------------------
# Forward selection

# environment
r2.all <- gam.multi.par(sp = gll.hat, pred = pca.env) #R2 with all variables
r2a.thresh <- 1 - (1 - r2.all) * (nrow(sp) - 1)/(nrow(sp) - ncol(pca.env) - 1) #Adjusted R2 with all variables (Threshold)
fw <- fw.sel.gam.par(sp = gll.hat, pred = pca.env, r2a.thresh = r2a.thresh, alpha = 0.05) #Forward selection with double-stopping criterion
fw #selected variables
env.sel <- as.data.frame(pca.env[, fw$variable]) #Selected variables matrix
R2X1<-gam.multi.par(sp = gll.hat, pred = env.sel)#R2 with the selected variables
R2X1
R2aX1<-adjR2X1.3(sp=gll.hat,env=env.sel,xy=xy,names=names,coords=coords)#Adjusted R2 using MSR method
R2aX1


# space
r2.all <- gam.multi.par(sp = gll.hat, pred = xy2) #R2 with all variables
r2a.thresh <- 1 - (1 - r2.all) * (nrow(sp) - 1)/(nrow(sp) - ncol(pca.env) - 1) #Adjusted R2 with all variables (Threshold)
fw <- fw.sel.gam.par(sp = gll.hat, pred = xy2, r2a.thresh = r2a.thresh, alpha = 0.05) #Forward selection with double-stopping criterion
fw #Forward selection with double-stopping criterion. 
spa.sel <- as.data.frame(xy[, fw$variable]) #Selected variables matrix. With very few variables, automatic selection can be troublesome. Manual selection might be needed
R2X2<-gam.multi.par(sp = gll.hat, pred = as.data.frame(spa.sel)) #R2 with the selected variables
R2X2 
R2aX2 <- 1 - (nrow(gll.hat) - 1) / (nrow(gll.hat) - ncol(xy2) - 1) * (1 - R2X2) #Adjusted R2 using Ezequiel formula
R2aX2

# time
r2.all <- gam.multi.par(sp = gll.hat, pred = time) #R2 with all variables
r2a.thresh <- 1 - (1 - r2.all) * (nrow(sp) - 1)/(nrow(sp) - ncol(pca.env) - 1) #Adjusted R2 with all variables (Threshold)
fw <- fw.sel.gam.par(sp = gll.hat, pred = time, r2a.thresh = r2a.thresh, alpha = 0.05) #Forward selection with double-stopping criterion
fw #Forward selection with double-stopping criterion. 
time.sel <- as.data.frame(time[, fw$variable]) #Selected variables matrix. With very few variables, automatic selection can be troublesome. Manual selection might be needed
R2X3<-gam.multi.par(sp = gll.hat, pred = as.data.frame(time.sel)) #R2 with the selected variables
R2X3
R2aX3 <- 1 - (nrow(gll.hat) - 1) / (nrow(gll.hat) - ncol(time) - 1) * (1 - R2X3) #Adjusted R2 using Ezequiel formula
R2aX3

#environment+space
envspa<-cbind(env.sel,spa.sel) #Combination of environmental and spatial variables
R2X1X2 <- gam.multi.par(sp = gll.hat, pred = as.data.frame(envspa)) #R2 with the selected variables
R2X1X2
R2X1X2.msr<-adjR2X1X2.3(sp=gll.hat,env=as.data.frame(env.sel),spa=as.data.frame(spa.sel),xy=xy2,names=names,coords=coords) #Adjusted R2 using MSR method

#space+time
xytime<-cbind(spa.sel,time.sel) #Combination of spatial and temporal variables
R2X2X3 <- gam.multi.par(sp = gll.hat, pred = as.data.frame(xytime)) #R2 with the selected variables
R2X2X3
R2aX2X3<-1 - (nrow(gll.hat) - 1) / (nrow(gll.hat) - ncol(xytime) - 1) * (1 - R2X2X3) #Adjusted R2 using Ezequiel formula

#environment+time
envtime<-cbind(env.sel,time.sel) #Combination of environmental and temporal variables
R2X1X3 <- gam.multi.par(sp = gll.hat, pred = as.data.frame(env.sel)) #R2 with the selected variables
R2X1X3
R2X1X3.msr<-adjR2X1X2.3(sp=gll.hat,env=as.data.frame(env.sel),spa=as.data.frame(time.sel),xy=time,names=names,coords=coords) #Adjusted R2 using MSR method

#environment+space+time
all<-cbind(env.sel,spa.sel,time.sel) #Combination of environmental, spatial and temporal variables
R2X1X2X3 <- gam.multi.par(sp = gll.hat, pred = all) #R2 with the selected variables
R2X1X2X3
R2X1X2X3.msr<-adjR2X1X2.3(sp=gll.hat,env=as.data.frame(env.sel),spa=as.data.frame(xytime),xy=xy,names=names,coords=coords) #Adjusted R2 using MSR method


#Adjusting...
a.W <- 1 - (1-(R2X1X2 - R2X2)) / (1-mean(R2X1X2.msr - R2X2))
b.W <- R2aX1 -  a.W
c.W <- R2aX2 - b.W
R2aX1X2 <- a.W + b.W + c.W

a.Z <- 1 - (1-(R2X1X3 - R2X3)) / (1-mean(R2X1X3.msr - R2X3))
b.Z <- R2aX1 -  a.Z
c.Z <- R2aX3 - b.Z
R2aX1X3 <- a.Z + b.Z + c.Z

#Partitioning
a <- 1 - (1 - (R2X1X2X3 - R2X2X3)) / (1-mean(R2X1X2X3.msr - R2X2X3))
all <- a + R2aX2X3
b <- all - R2aX1X3
c <- all - R2aX1X2 
d <- all - R2aX3 - a - b
e <- all - R2aX1 - b - c
f <- all - R2X2 - a - c
g <- all - a - b - c - d - e - f  
h <- 1 - all 

vp <- c(a, b, c, d, e, f, g, h) #Vector of the fractions
names(vp) <- letters[1:length(vp)]
vp
vpdat<-as.data.frame(vp)
rownames(vpdat)<-c("Pure environment","Pure space", "Pure time", "Environment-Space", "Space-Time", "Environment-Time", "Environment-Space-Time", "Residuals")
vpdat