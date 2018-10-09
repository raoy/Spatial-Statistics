require(gstat)
require(sp)

data(meuse)
names(meuse)
summary(meuse)
hist(meuse$zinc)
plot(meuse$x,meuse$y,asp=1)

str(meuse)
coordinates(meuse) = ~x + y # convert data to spatial points data frame
class(meuse)
data(meuse.grid)# this is gridded data for the meuse data set
str(meuse.grid)
coordinates(meuse.grid) = ~x + y # converts to spatial class
gridded(meuse.grid) <- TRUE

#Inverse Distance Weighting Interpolation (IDW)
data(meuse.riv)
meuse.sp <- SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)),
"meuse.riv")))
meuse.lt = list(riv = list("sp.polygons", meuse.sp, fill = "grey"),
pts = list("sp.points", meuse, pch = 3, cex = 0.5, col = "black"))

#Doing IDW interpolation in gstat
zn.idw = idw(log(zinc) ~ 1, meuse, meuse.grid)
spplot(zn.idw, "var1.pred", sp.layout = meuse.lt, main = "log(zinc),
inverse distance interpolation-default")

####silahkan disalin yang di bawah ini:
#evaluating IDW interpolation
zn.idwcv<-krige.cv(log(zinc) ~ 1, meuse)
RMSE.idw<-sqrt(sum(zn.idwcv$residual^2)/length(zn.idwcv$residual))
#######################################################################

#Try it again but set the minimum number of points in the search neighbourhood to 6, a minimum
zn.idw6 <- krige(log(zinc) ~ 1, meuse, meuse.grid,nmax=6)
spplot(zn.idw6, "var1.pred", sp.layout = meuse.lt, main = "log(zinc),
inverse distance interpolation-nmx=6")

#Try it with two different weighting exponents: 0.5,1,5 and 10
meuse.grid$idp05 = idw(log(zinc) ~ 1, meuse, meuse.grid, idp = 0.5)$var1.pred
meuse.grid$idp1 = idw(log(zinc) ~ 1, meuse, meuse.grid, idp = 1)$var1.pred
meuse.grid$idp5 = idw(log(zinc) ~ 1, meuse, meuse.grid, idp = 5)$var1.pred
meuse.grid$idp10 = idw(log(zinc) ~ 1, meuse, meuse.grid, idp = 10)$var1.pred
spplot(meuse.grid, c("idp05", "idp1", "idp5", "idp10"), sp.layout =
 meuse.lt, main = "log(zinc), inverse distance interpolation")

#Thiessen polygons in gstat
zn.tp = krige(log(zinc) ~ 1, meuse, meuse.grid, nmax = 1)#for this search
# the neighborhood is set to nmax=1
image(zn.tp["var1.pred"])
points(meuse, pch = "+", cex = 0.65)
cc = coordinates(meuse)
library(tripack)
plot(voronoi.mosaic(cc[, 1], cc[, 2]), do.points = FALSE, add = TRUE)
title("Thiessen (or Voronoi) polygon interpolation of log(zinc)")



#Kriging in gstat
coordinates(meuse) = ~x + y # convert data to spatial points data frame
class(meuse)
gridded(meuse.grid) = TRUE
class(meuse.grid)
vgm1 <- variogram(log(zinc)~1, meuse)
plot(vgm1, plot.numbers = TRUE, pch = "+")
m <- fit.variogram(vgm1,vgm(.59,"Sph",896,.05))
plot(vgm1, model=m)
show.vgms() # various models for fitting variogram in gstat
m <- vgm(.59, "Sph", 874, .04)

#Perform ordinary kriging:
x.k <- krige(log(zinc)~1, meuse, meuse.grid, model = m)
spplot(x.k["var1.pred"], main = "Meuse zinc ordinary kriging log predictions")
x.k$sek <- sqrt(x.k$var1.var)
spplot(x.k,zcol="sek", main = "Meuse zinc ordinary kriging se")
summary(x.k)

#Evaluating ordinary kriging:
zn.okcv<-krige.cv(log(zinc)~1, meuse, meuse.grid, model = m)
RMSE.ok<-sqrt(sum(zn.okcv$residual^2)/length(zn.okcv$residual))

# Add a covariate
data(meuse.riv)# outline of the river
meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.sr <- SpatialPolygons(meuse.lst)
image(meuse.grid["dist"])# one of the variables in meuse.grid
plot(meuse.sr, add=TRUE)
title("distance to river")
par(mfrow=c(1,2))
plot(zinc ~ dist, meuse)
plot(log(zinc) ~ sqrt(dist), meuse)
abline(lm(log(zinc) ~ sqrt(dist), meuse))

zn.lm <- lm(log(zinc) ~ sqrt(dist), meuse)
meuse$fitted.s <- predict(zn.lm,meuse)-mean(predict(zn.lm,meuse))
meuse$residuals <- residuals(zn.lm)
spplot(meuse,c("fitted.s","residuals"))

# generate a new kriged surface with the covariate - distance to river
vm.uk <- variogram(log(zinc)~sqrt(dist), meuse)
plot(vm.uk, plot.numbers = TRUE, pch = "+")
m.uk <- fit.variogram(vm.uk,vgm(.3,"Sph",800,.06))
plot(vm.uk, model=m.uk)
m.uk <- vgm(.159, "Sph", 880, .05)
ko.uk <- krige(log(zinc)~ x+y, meuse,meuse.grid,model=m.uk)
pts <- list("sp.points",meuse,pch=3,col="black")
meuse.layout <- list(pts)
spplot(ko.uk["var1.pred"], sp.layout=meuse.layout, main = "universal
kriging predictions-Zn/distance river ")
ko.uk$sek <- sqrt(ko.uk$var1.var)
spplot(ko.uk,zcol='sek', sp.layout=meuse.layout, main = "universal
kriging se-Zn(covariate)")
summary(ko.uk)

zn.ukcv<-krige.cv(log(zinc)~ x+y, meuse,meuse.grid,model=m.uk)
RMSE.uk<-sqrt(sum(zn.ukcv$residual^2)/length(zn.ukcv$residual))


#back transform the results
ko.uk$predt <- exp(ko.uk$var1.pred)
spplot(ko.uk["predt"], sp.layout=meuse.layout, main = "UK covariate
predictions-Meuse zinc log/backtrans(Zn)")

#Plot all predictions next to each other:
meuse.grid$zn.id <- zn.idw$var1.pred
meuse.grid$zn.ok <- x.k$var1.pred
meuse.grid$zn.uk <- ko.uk$var1.pred
spplot(meuse.grid[c("zn.ok", "zn.uk", "zn.id")])

#Cross Validation
kcv <- krige.cv(log(zinc) ~ 1, meuse, model=m, nfold=10)
summary(kcv)

#Some statistics:
# mean error, ideally 0:
mean(kcv$residual)
(rmse <- sqrt(mean(kcv$residual^2)))
rmse/sd(meuse$zinc)
#correlation observed and predicted, ideally 1
cor(kcv$observed, kcv$observed - kcv$residual)
#mean(k.cv$residual^2/k.cv$var1.var)
#correlation predicted and residual, ideally 0
cor(kcv$observed - kcv$residual, kcv$residual)
bubble(kcv, "zscore", main = "zinc:10-fold CV-zscores")
#zscores-residual/kriging standard error

#The same for UK with a covariate:
kcv.uk <- krige.cv(log(zinc) ~ sqrt(dist), meuse, model=m.uk, nfold=10)
summary(kcv.uk)
mean(kcv.uk$residual)
mean(kcv.uk$residual^2)
(rmse <- sqrt(mean(kcv.uk$residual^2)))
rmse/sd(meuse$zinc)
cor(kcv.uk$observed, kcv.uk$observed - kcv.uk$residual)
cor(kcv.uk$observed - kcv.uk$residual, kcv.uk$residual)
bubble(kcv.uk, "zscore", main = "zinc with distance to river:10-fold CVzscores")

