rm(list=ls())
setwd("D:/Dept.STK/Statistika Spasial (S2)/TA 2017-2018/Praktkum RA/Reference/")
library(rgdal);library (spdep)
data(columbus)

NY8 <- readOGR("spatial-data-with-r-master", "NY8_utm18")
TCE <- readOGR("spatial-data-with-r-master", "TCE")
cities <- readOGR("spatial-data-with-r-master", "NY8cities")
NY_nb <- read.gal("spatial-data-with-r-master/NY_nb.gal", region.id=row.names(NY8))

Modeling single variable and Checking the regression residuals
col.listw <- nb2listw(col.gal.nb)
columbus.lm0<- lm(CRIME ~ 1, data=columbus)
summary(columbus.lm0)

col.moran0 <-  lm.morantest(columbus.lm0,col.listw); col.moran0

col.e <- resid(columbus.lm0)
col.moran_0 <- moran.test(col.e,col.listw, randomisation=FALSE); col.moran_0



