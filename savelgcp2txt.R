### This file is to output data from LGCP simulation to csv/txt files

#+ dirNames

#dataDir = '/store/patrick/lgcp'
dataDir = "Y:/temp/personal/tech/R/simulationptoronto"
#'

#+ setupKnitr, include=FALSE
library('knitr')
kFile = 'Y:/temp/personal/tech/R/geostatsp/knitrUtils.R'
#if(!file.exists(kFile)) {
#  download.file("https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/utils/R/knitr.R?root=diseasemapping", kFile)
#}

source(kFile)
knitr::knit_hooks$set(plot=hook_plot_p) 


utilsFile = 'Y:/temp/personal/tech/R/geostatsp/lgcpUtils.R'
#if(!file.exists(utilsFile)) {
#  download.file(
#    'https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/utils/R/lgcp.R?root=diseasemapping', 
#    utilsFile)
#}
source(utilsFile)

#'






#' # Compare lgcp to geostatsp

#+ packages
library('geostatsp')
library("lgcp")
library('mapmisc')
#'


savelgcp2txt = function(datadir, srcdata, surfix="") {
  count = srcdata$inla$.args$data$count[order(srcdata$inla$.args$data$space)]
  x = srcdata$inla$.args$data$x[order(srcdata$inla$.args$space)]
  randomeff = srcdata$raster$random.mean[order(srcdata$inla$.args$space)]
  AmatColumn = c(rep(1, length(x)),x)
  write.table(count, file = paste0(datadir, "/y",surfix,".txt"), sep="\t",col.names=F, row.names=F)
  write.table(AMatColumn, file = paste0(datadir, "/AMatColumn",surfix,".txt"), sep="\t", col.names=F, row.names=F)
  write.table(rep(1, length(count)), file = paste0(datadir,"/O", surfix,".txt"), sep="\t", col.names=F, row.names=F)
  write.table(rep(1, length(count)), file = paste0(datadir, "x", surfix, ".txt"), sep="\t", col.names=F, row.names=F)
  
}


#' 
rasterImage = meuse.raster.aggregate
thebb = bbox(rasterImage)
xseq = seq(thebb[1,1], thebb[1,2], len=dim(rasterImage)[1])
yseq = seq(thebb[2,1], thebb[2,2], len=dim(rasterImage)[2])
xseq = seq(1, 35, len=dim(rasterImage)[1])
yseq = seq(1, 26, len=dim(rasterImage)[2])
Ncoords = prod(dim(rasterImage)[1:2])
Ncoordsx = dim(rasterImage)[1]
Ncoordsy = dim(rasterImage)[2]
coordsMat = matrix(xseq, Ncoordsx, Ncoordsy) + 1i * matrix(yseq, Ncoordsx, Ncoordsy, byrow = T)
coordsVec = as.vector(coordsMat)
distMat = Mod(outer(coordsVec, coordsVec, FUN="-"))

#' 
#' 
gridMat = outer(yseq * 1i, xseq, FUN="+")
coordsGrid = as.vector(gridMat)
distGrid = Mod(outer(coordsGrid, coordsGrid, FUN="-"))
#' 
#' 
meshgrid <- function(a,b) {
  list(
   x = outer(b=0, a, FUN='+')
   ,y = outer(b, a=0, FUN='+')
  )
}
#' 


#' #April 24, write to text file
#' From file lgcpPatrickBen_2.R
DIR="Y:/temp/personal/tech/R/simulationp_p1"
DIR = dataDir
load(file=paste0(DIR, "/pRes.RData"))
fstFit = pRes
count = fstFit$inla$.args$data$count
x = fstFit$inla$.args$data$x
randomeff = fstFit$raster$random.mean[order(fstFit$inla$.args$data$space)]
AMatColumn = c(rep(1, length(x)),x)

Surfix=""
write.table(count,file=paste0(DIR, "/y",Surfix, ".txt"), sep="\t", col.names=F, row.names=F)
write.table(AMatColumn, file = paste0(DIR, "/AMatColumn", Surfix,".txt"), sep="\t", col.names=F, row.names=F)
write.table(rep(1,length(count)), file=paste0(DIR, "/O", Surfix, ".txt"), sep="\t", col.names=F, row.names=F)
write.table(rep(1,length(count)), file=paste0(DIR, "/x", Surfix, ".txt"), sep="\t", col.names=F, row.names=F)
write.table(randomeff, file=paste0(DIR, "/randomeff", Surfix, ".txt"), sep="\t", col.names=F, row.names=F)
#' 
#' 
#' # Write the distance matrix
rasterImage = fstFit$raster
Ncoordsx = dim(rasterImage)[2]
Ncoordsy = dim(rasterImage)[1]
xseq = seq(1,dim(rasterImage)[2], len=dim(rasterImage)[2])
yseq = seq(1,dim(rasterImage)[1], len=dim(rasterImage)[1])
coordsMat = matrix(xseq, Ncoordsx, Ncoordsy) + 1i * matrix(yseq, Ncoordsx, Ncoordsy, byrow = T)
coordsVec = as.vector(coordsMat)
distMat = Mod(outer(coordsVec, coordsVec, FUN="-"))
write.table(distMat, file=paste0(Dir, "/distanceR30C17", Surfix, ".txt"), sep="\t", col.names=F, row.names=F)
#' 
#' 

#' 
#' How to calculate interval of variance
#' The following function shows how to figure out parameters so that
#' it has a similar distribution as other methods 
#' --e.g. One is scale, the other is range --ZF Jan. 14, 2023
#' You can change weight to make sure the distribution covers 97.5%,2.5% quantile
#'  
objFun = function(par) {
     par1 = par[1]
     par2 = par[2]
     qlower = pgamma(slower, shape=par1, rate=par2)
     qupper = pgamma(supper, shape=par1, rate=par2)
     Mode = (par1-1)/par2
     #(qupper-0.975)^2 + (qlower-0.025)^2 + (qupper-qlower-0.95)^2
     4*(qupper-0.975)^2 + (qlower-0.025)^2
}

intSd
slower = (intSd**2)[1]
supper = (intSd**2)[2]
resr = optim(c(slower, supper), objFun);resr$par
#[1] 3.240869 9.377425
#[1] 3.240869 9.377425

intRange
slower = 4/intRange[2]
supper = 4/intRange[1]
resr = optim(c(slower, supper), objFun);resr$par
#[1]     76.05856 320367.91945

resr = optim(resr$par, objFun);resr$par
#' 
#' 
#' 
#' 
setwd(DIR)
load(file='lgcpPatrickBen.RData')
intRangeS = intRange/1422.493   ####xres(fstFit$raster)
intScale = 4/intRangeS
slower = intScale[2]
supper = intScale[1]
resr = optim(c(slower, supper), objFun);resr$par
rgammrange = rgamma(1000, resr$par[1], resr$par[2])
#[1] 11.58874 46.73769
#[1]  6.72168 10.17067 -- for quantile [1]  4565.826 21901.841
rgammarange = 4* 1422.493/rgamma(1000, resr$par[1], resr$par[2]) 
plot(density(rgammarange),xlim=c(1000,45000))
abline(v=intRange[1]);abline(v=intRange[2])
#' 
#' 
intRangeS = intRange/(1422.493/4);intRangeS
intScale = 4/intRangeS
slower = intScale[2];slower
supper = intScale[1];supper
resr = optim(c(slower, supper), objFun);resr$par
#[1]   11.59105 4673.97422 #There were 40 warnings (use warnings() to see them)
rgammarange100 = 1422.493/rgamma(1000, resr$par[1], resr$par[2])
plot(density(rgammarange100))
#' 
#' 
intRangeS = intRange/(1422.493/125);intRangeS
intScale = 4/intRangeS
slower = intScale[2];slower
supper = intScale[1];supper
resr = optim(c(slower, supper), objFun);resr$par
#[1]   11.59181 5844.70457  #There were 12 warnings (use warnings() to see them)
rgammarange125 = rgamma(1000, resr$par[1], resr$par[2])
plot(density(rgammarange125))
#' 
#' 
#' 
intRangeS = intRange/(1422.493/150);intRangeS
intScale = 4/intRangeS
slower = intScale[2];slower
supper = intScale[1];supper
resr = optim(c(slower, supper), objFun);resr$par
#[1]  0.8091434 43.3868189
rgammarange150 = rgamma(1000, resr$par[1], resr$par[2])
plot(density(rgammarange150))
#' 
#' 
#' 
intRangeS = intRange/(1422.493/200);intRangeS
intScale = 4/intRangeS
slower = intScale[2];slower
supper = intScale[1];supper
resr = optim(c(slower, supper), objFun);resr$par
#[1] 1.543168e-02 5.385352e-18
rgammarange200 = rgamma(1000, resr$par[1], resr$par[2])
#' 
#' 
#' 
#comparison of all these range densities
plot(density(4/rgammarange100*(1422.493/100)))
lines(density(4/rgammarange150*(1422.493/150)), col='red')

#' 
#' 
#' 
#' #comparison of range prior settings
# When distance is scaled down by its resolution
rgammarange = rgamma(1000, 100.40604, 87.38268)
rgammarange100 = rgamma(1000, 98.66471, 8536.25592)
rgammarange1000 = rgamma(1000, 98.66471, 85362.5592)
rgammarangeO = rgamma(1000, 98.66471, 8536.25592/100*1422.493)

plot(density(4*1422.493/rgammarange))
lines(density(4*1422.493/(rgammarange100*100)), col='red')
lines(density(4*1422.493/(rgammarange1000*1000)), col='green')
lines(density(4*1422.493/(rgammarangeO*1422.493)), col='blue')
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
# The following is to calculate based on grid--
# we will explore unit grid vs grid 
xres(fstFit$raster)
(fstFit$raster@extent@xmax-fstFit$raster@extent@xmin)/xres(fstFit$raster)
(fstFit$raster@extent@ymax-fstFit$raster@extent@ymin)/xres(fstFit$raster)
intRange
intRangeGrid = intRange/xres(fstFit$raster) 
intRangeGrid

intScale = 4/intRangeGrid
slower = intScale[2];slower
supper = intScale[1];supper
resr = optim(c(slower, supper), objFun);resr$par
#### [1] 100.40606  87.38269
rgammarangeGrid = rgamma(1000, resr$par[1], resr$par[2])
plot(density(rgammarangeGrid))
abline(v=c(slower,supper))
1/30/17   # 0.001960784
#' 
#' 
#' 
#' 
# This is to explore if Unit grid's difference 
intRange
intRangeGridUnit = intRange/xres(fstFit$raster)/30 
intRangeGridUnit

intScale = 4/intRangeGridUnit
intScale
slower = intScale[2];slower
supper = intScale[1];supper
resr = optim(c(slower, supper), objFun);resr$par
#### [1] 100.40606  87.38269
rgammarangeGridUnit = rgamma(1000, resr$par[1], resr$par[2])
plot(density(rgammarangeGridUnit))
abline(v=c(slower,supper))
1/30/17   # 0.001960784
rgammarangeGridUnit = rgamma(1000, 100.40606, 3.869371)
sqrt(1/30/17) # 0.04428074
rgammarangeGridUnit = rgamma(1000, 100.40606, 2.912756)
plot(density(rgammarangeGridUnit))
abline(v=c(slower,supper))
sqrt(1/30/30) # 0.03333333
#' 
#' 
#' 
#' 