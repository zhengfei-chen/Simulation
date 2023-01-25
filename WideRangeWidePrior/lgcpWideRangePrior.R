#' # Compare lgcp to geostatsp
#' #This is to use permissive priors to test Ben's method
#' #Simulate the data with a large range
#' 

#+ packages
library('geostatsp')
library("lgcp")
library('mapmisc')
#'

#+ dirNames
dataDir = "C:/Users/542563036/Documents/personal/R"   
#'

#+ crs, cache=TRUE
#myCrs = omerc(geocode("Toronto, ON"), angle=-18)
load(paste0(dataDir,"/myCrs.RData"))
#'

#+ setupKnitr, include=FALSE
library('knitr')
#'


#+ crs, cache=TRUE
#source("Y:/temp/personal/tech/R/simulationptoronto/lgcpToGeostatsp.R")
#load("Y:/temp/personal/tech/R/simulationp/myCrs.RData")
#'


#+ setupSim, cache=TRUE
data(murder, package='geostatsp')
border = spTransform(torontoBorder, myCrs)

# approximate number of events
border$atrisk = 500

myRaster = squareRaster(border, 800)

myCovariate = myRaster
values(myCovariate) = matrix(seq(-1, 1, len=nrow(myRaster)), 
                             nrow(myRaster), ncol(myRaster))

myOffset = myRaster
values(myOffset) = log(border$atrisk)-sum(log(abs(apply(bbox(myRaster),1,diff))))

covList = list(x = myCovariate, offset = myOffset)
#'


#+ bgMap, cache=TRUE
myBg = openmap(myRaster, path='stamen-toner', buffer=2000, fact=1.5)
myBg = tonerToTrans(myBg)
#'


#+ sim, cache=TRUE, stuff=1
myModel = c(
  mean=0, 
  variance=0.4^2, 
  range=25000, shape=2)

set.seed(0)
myLgcp=simLgcp(param=myModel, 
               betas = c(x=0.4), offset = 'offset', 
               covariates=covList,
               rasterTemplate=myRaster)
eventsSample5k = myLgcp$events[
  sample(length(myLgcp$events), 250)]
#'


#+ plotSim, fig.cap='simulation', fig.height=3, fig.width=5, fig.subcap = c('covariate','events'), fig.ncol=2, echo=FALSE
map.new(myCovariate, TRUE)
plot(myCovariate, add=TRUE)
#plot(myBg, add=TRUE)

map.new(myCovariate, TRUE)
#plot(myBg, add=TRUE)

plot(eventsSample5k, cex=0.5,add=TRUE, col='#FF0000') #### change '#FF000030' to '#FF0000' --April 7, 2022 
#'


#+ priors
cf <- CovFunction(
  RandomFieldsCovFct(
    model="matern",additionalparameters=2
  )
)

priors <- lgcpPrior(
  etaprior=PriorSpec(LogGaussianPrior(
    mean=log(c(0.5, 25000/2)),
    variance=diag(c(0.3, 0.3)))),
  betaprior=PriorSpec(GaussianPrior(
    mean=rep(0,2),
    variance=diag(10^6,2))))

intSd = exp(qnorm(c(0.025, 0.975), 
                  priors$etaprior$mean[1], 
                  priors$etaprior$variance[1,1]))

intRange = exp(qnorm(c(0.025, 0.975), 
                     priors$etaprior$mean[2], 
                     priors$etaprior$variance[2,2]))*2
intSd
intRange
#'



#' #+ resPatrick
#' pResFile = file.path(dataDir,'pRes.RData')
#' if(!file.exists(pResFile)) {
#'   pRes = geostatsp::lgcp(
#'     data=eventsSample5k, 
#'     formula = ~ x,
#'     shape=environment(cf)$additionalparameters,
#'     grid=squareRaster(myLgcp$raster, 30), 
#'     buffer=4000,
#'     covariates=covList,
#'     priorCI=list(
#'       sd = intSd,
#'       range = intRange
#'     )
#'   )
#'   save(pRes, file=pResFile)
#' } else {
#'   load(pResFile)
#' }
#' #'




#+ packageForConverting, echo=TRUE, message=FALSE
library('spatstat')
library('maptools')
#'


#+ objectsBen

myCov2 = aggregate(myCovariate, 4)
spcovs <- as(myCov2,"SpatialPixelsDataFrame")
spcovs@data <- guessinterp(spcovs@data)


boundary = as(extent(myLgcp$raster), 'SpatialPolygons')
boundary = SpatialPolygonsDataFrame(
  boundary, data.frame(offset=rep(1, length(boundary)))
)
boundary@data = assigninterp(df = boundary@data, vars = "offset", 
                             value = "ArealWeightedSum")

eventsPPP = as.ppp(
  eventsSample5k@coords,
  boundary)

#'

#+ overlays
polyolayfl <- file.path(dataDir,"polyolay.RData")
if(!file.exists(polyolayfl)){
  
  polyolay <- getpolyol(
    data=eventsPPP,
    regionalcovariates=boundary,
    pixelcovariates=spcovs,
    #        cellwidth=2000,
    cellwidth=1000,
    ext=2)
  save(polyolay,file=polyolayfl)
  
} else {
  load(polyolayfl)
}

zfl <- file.path(dataDir,"Zmat.RData")
if(!file.exists(zfl)){  # this takes less time than the above, but to speed things up you can save it
  # if you use a model with a different set of covariates you will have to re-run 
  # this bit
  Zmat <- getZmat(
    formula= X ~ layer,
    data=eventsPPP,
    pixelcovariates=spcovs,
    overl=polyolay
  ) 
  save(Zmat,file=zfl)
  
} else {
  load(zfl)
}

#'

#+ testZmat, eval=FALSE, include=FALSE
zArray = array(
  Zmat, c(
    c(polyolay$gridobj$M, polyolay$gridobj$N)/polyolay$ext,	
    ncol(Zmat)	
  )
)
image(zArray[,,1])		
image(zArray[,,2])		
#'



#' fit the model with mcmc
#+ lgcpBen
bResFile = file.path(dataDir, 'bRes.RData')
if(!file.exists(bResFile)) {
  bRes <- lgcpPredictSpatialPlusPars(   
    formula=attributes(Zmat)$FORM,
    sd=eventsPPP, 
    Zmat=Zmat,
    model.priors=priors,
    spatial.covmodel=cf,
    cellwidth=polyolay$cellwidth,
    mcmc.control=mcmcpars(
      mala.length=3000000,burnin=100000, retain=2900,
      #					mala.length=10000,burnin=1000, retain=50,
      adaptivescheme=andrieuthomsh(
        inith=1,alpha=0.5,C=1,
        targetacceptance=0.574)),
    output.control=setoutput(
      gridfunction=
        dump2dir(
          dirname= dataDir,
          forceSave=TRUE)),
    ext=polyolay$ext,
    poisson.offset=NULL)
  
  save(bRes, file=bResFile)
  
} else {
  load(bResFile)
}

warnings()
bRes$mcmcacc
bResP = lgcpToGeostatsp(bRes)
#'

##### comment out the following by Ctrl-Shift-C

#' 
#' #+ parTables, results='asis', echo=FALSE
#' knitr::kable(bResP$summary, digits=1)
#' 
#' knitr::kable(pRes$parameters$summary[,
#'                                      colnames(bResP$summary)], digits=1)
#' #'
#' 
#' #+ tracePlots, fig.cap='trace plots', fig.subcap = c('intercept','beta','sd','range'), fig.height=5, fig.width=5, fig.ncol=2, echo=FALSE
#' 
#' plot(bRes$betarec[,1], xlab='iter', ylab='intercept',type='l')
#' 
#' plot(bRes$betarec[,2], xlab='iter', ylab='beta',type='l')
#' 
#' plot(bResP$samples$sd, xlab='iter', ylab='sd',type='l', log='y')
#' 
#' plot(bResP$samples$range, xlab='iter', ylab='range',type='l', log='y')
#' 
#' #'
#' 
#' #+ priorPost, fig.cap='Posterior distributions', fig.subcap = c('intercept','beta','sd','range'), fig.height=5, fig.width=5, fig.ncol=2, echo=FALSE
#' 
#' 
#' hist(bResP$samples$beta[,'(Intercept)'], xlab='intercept', prob=TRUE, main='')
#' lines(pRes$inla$marginals.fixed[['(Intercept)']], col='red')
#' abline(v=myLgcp$parameters$fixed['intercept'], col='yellow')
#' legend("topright", lty=1, col=c('black','red'), legend=c('b','p'))
#' 
#' 
#' hist(bResP$samples$beta[,2], xlab='beta', prob=TRUE, main='',
#'      xlim=range(pRes$inla$marginals.fixed[['x']][,1]))
#' lines(pRes$inla$marginals.fixed[['x']], col='red')
#' abline(v=myLgcp$parameters$fixed['x'], col='yellow')
#' legend("topright", lty=1, col=c('black','red','yellow'), legend=c('b','p','truth'))
#' 
#' 
#' 
#' plot(bResP$parameters$sd$posterior, xlab='sd', ylab='dens',
#'      xlim = range(c(
#'        range(bResP$parameters$sd$posterior[,'x']),
#'        range(pRes$parameters$sd$posterior[,'x']))),
#'      ylim = range(c(
#'        range(bResP$parameters$sd$posterior[,'y']),
#'        range(pRes$parameters$sd$posterior[,'y']))),
#'      main='', type='l')
#' 
#' lines(pRes$parameters$sd$posterior, col='red', lwd=2)
#' lines(pRes$parameters$sd$posterior[,c('x','prior')], col='red', lty=2)
#' #lines(pRes$parameters$sd$prior, lty=2,col='red')
#' 
#' lines(bResP$parameters$sd$prior,lty=2,col='black')
#' abline(v=sqrt(myLgcp$parameters$random['variance']), col='yellow')
#' legend("topright", lty=c(1,1,1,2,1), col=c('black','red','black','black','yellow'), 
#'        legend=c('b','p', 'post','prior','truth'))
#' 
#' 
#' plot(bResP$parameters$range$posterior, xlab='range', ylab='dens',
#'      xlim = range(c(
#'        range(bResP$parameters$range$posterior[,'x']),
#'        range(pRes$parameters$range$posterior[,'x']))),
#'      ylim = range(c(
#'        range(bResP$parameters$range$posterior[,'y']),
#'        range(pRes$parameters$range$posterior[,'y']))),
#'      main='', type='l')
#' lines(pRes$parameters$range$posterior, col='red', lwd=2)
#' lines(pRes$parameters$range$posterior[,c('x','prior')], lty=2, col='red')
#' #lines(pRes$parameters$range$prior, lty=2,col='red')
#' lines(bResP$parameters$range$prior,
#'       lty=2,col='black')
#' abline(v=myLgcp$parameters$random['range'], col='yellow')
#' legend("topright", lty=c(1,1,1,2,1), col=c('black','red','black','black','yellow'), 
#'        legend=c('b','p', 'post','prior','truth'))
#' 
#' #'
#' 
#' 
#' 
#' #+ twod, eval=FALSE
#' stuff = KernSmooth::bkde2D(cbind(bResP$samples$sd, bResP$samples$range), bandwidth=c(0.01, 100))
#' image(matrix(stuff$fhat, 51, 51))
#' 
#' 
#' #'
#' 
#' #+ plotRandom, fig.cap='residual spatial variation', fig.subcap = c('truth','geostatsp','lgcp'), fig.height=3.5, fig.width=7, echo=FALSE
#' colR = colourScale(bResP$samples$random,
#'                    col='Spectral', breaks=12, dec=-log10(0.25),
#'                    rev=TRUE, style='equal')
#' 
#' map.new(myLgcp$raster,TRUE)
#' plot(myLgcp$raster$random, 
#'      col=colR$col, breaks=colR$breaks, legend=FALSE,
#'      add=TRUE)
#' #plot(myBg, add=TRUE)
#' scaleBar(myBg, 'bottomright')
#' legendBreaks('topright', colR)
#' 
#' map.new(myLgcp$raster,TRUE)
#' plot(pRes$raster$random.mean, 
#'      col=colR$col, breaks=colR$breaks, legend=FALSE,
#'      add=TRUE)
#' #plot(myBg, add=TRUE)
#' scaleBar(myBg, 'bottomright')
#' legendBreaks('topright', colR)
#' 
#' map.new(myLgcp$raster,TRUE)
#' plot(mean(bResP$samples$random), 
#'      col=colR$col, breaks=colR$breaks, legend=FALSE,
#'      add=TRUE)
#' #plot(myBg, add=TRUE)
#' scaleBar(myBg, 'bottomright')
#' legendBreaks('topright', colR)
#' 
#' 
#' #'
#' 
#' #+ randomSamples, fig.cap='random samples', fig.height=3, fig.width=6, fig.ncol=2, echo=FALSE
#' for(D in seq(1, nlayers(bResP$samples$random), len=4)) {
#'   
#'   map.new(bResP$samples$random,TRUE)
#'   plot(bResP$samples$random[[D]], 
#'        col=colR$col, breaks=colR$breaks, legend=FALSE,
#'        add=TRUE)
#'   #	plot(myBg, add=TRUE)
#'   scaleBar(myBg, 'bottomright')
#'   
#'   legendBreaks('topright', colR)
#'   
#'   
#' }
#' 
#' 
#' #'
#' 
#' #+ rescaleModel
#' 
#' offsetTruth = myLgcp$parameters$fixed['intercept']
#' 
#' lambdaTruth = myLgcp$raster[['intensity']] * 
#'   exp(-offsetTruth)
#' 
#' totalArea = exp(sum(log(abs(
#'   apply(bbox(myLgcp$raster),1,diff)
#' ))))
#' 
#' offsetSim = log(length(eventsSample5k)/totalArea)
#' lambdaP = pRes$raster$predict.exp *exp(-offsetSim)
#' lambdaB = mean(bResP$samples$intensity)  * exp(-offsetSim)
#' #'
#' 
#' 
#' #+ testScaling, eval=FALSE, include=FALSE
#' 
#' length(myLgcp$events) / length(eventsSample5k)
#' 
#' mean(values(lambdaTruth))
#' 
#' mean(values(myLgcp$raster[['intensity']])) * totalArea
#' mean(values(lambdaTruth)*exp(offsetTruth))*totalArea
#' 
#' exp(offsetTruth)*totalArea
#' 
#' 
#' mean(values(lambdaP))
#' mean(values(lambdaB))
#' mean(values(pRes$raster$predict.exp )) * totalArea
#' mean(values(lambdaP)*exp(offsetSim))*totalArea
#' 
#' #'
#' 
#' 
#' #+ plotIntensity, fig.cap='intensity', fig.subcap=c('truth','geostatsp','lgcp'), fig.height=3.5, fig.width=7, echo=FALSE
#' colI = colourScale(myLgcp$raster[['relativeIntensity']],
#'                    col='Spectral', breaks=10, dec=-log10(0.5),
#'                    rev=TRUE, style='equal')
#' 
#' map.new(myLgcp$raster,TRUE)
#' plot(myLgcp$raster[['relativeIntensity']], col=colI$col, 
#'      breaks=colI$breaks, legend=FALSE,
#'      add=TRUE)
#' #plot(myBg, add=TRUE)
#' points(eventsSample5k, cex=0.4, col='#FF000030')
#' scaleBar(myBg, 'bottomright')
#' legendBreaks('topright', colI)
#' 
#' map.new(myLgcp$raster,TRUE)
#' plot(lambdaP, col=colI$col, 
#'      breaks=colI$breaks, legend=FALSE,
#'      add=TRUE)
#' #plot(myBg, add=TRUE)
#' scaleBar(myBg, 'bottomright')
#' legendBreaks('topright', colI)
#' 
#' map.new(myLgcp$raster,TRUE)
#' plot(lambdaB, col=colI$col, 
#'      breaks=colI$breaks, legend=FALSE,
#'      add=TRUE)
#' #plot(myBg, add=TRUE)
#' scaleBar(myBg, 'bottomright')
#' legendBreaks('topright', colI)
#' #'
#' 
#' 
#' #+ intensitySamples, fig.cap='intensity samples', fig.height=3, fig.width=6, fig.ncol=2, echo=FALSE
#' breaksScaled = colI$breaks*exp(offsetSim)
#' 
#' for(D in seq(1, nlayers(bResP$samples$intensity), len=4)) {
#'   
#'   map.new(myLgcp$raster,TRUE)
#'   plot(bResP$samples$intensity[[D]], 
#'        col=colI$col, breaks=breaksScaled, legend=FALSE,
#'        add=TRUE)
#'   #	plot(myBg, add=TRUE)
#'   scaleBar(myBg, 'bottomright')
#'   
#'   legendBreaks('topright', colI)
#'   
#'   
#' }
#' 
#' 
#' #'
#' 
