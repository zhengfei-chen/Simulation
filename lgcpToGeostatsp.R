#### Added beta parameter in the return list -- Zhengfei Mar. 2022

#+ convertBres, include=FALSE,echo=TRUE,eval=FALSE
lgcpToGeostatsp = function(x) {
  # convert lgcp results to geostatsp-type stuff
  
  bRes = x
  
  randomFile = file.path(bRes$gridfunction$dirname, "simout.nc")
  fixedFile = file.path(bRes$gridfunction$dirname, 'fixed.nc')
  intensityFile = file.path(bRes$gridfunction$dirname, 'intensity.nc')
  summaryFile = file.path(bRes$gridfunction$dirname, 'summary.grd')
  
  random = brick(randomFile)
  random = setMinMax(random)
  extent(random) = extent(
    c(bRes$xyt$window$xrange,
      bRes$xyt$window$yrange)
  )
  
  smallZ = array(
    bRes$Z,
    c(bRes$N, bRes$ext, bRes$M, bRes$ext, ncol(bRes$Z))
  )
  smallZ = drop(smallZ[,1,,1,])
  smallZ = matrix(smallZ, ncol=ncol(bRes$Z))
  smallZ[as.vector(bRes$cellInside)==0, ] = NA
  
  fixed = tcrossprod(smallZ, bRes$betarec)
  fixed = array(fixed, 
                c(bRes$M, bRes$N, nrow(bRes$betarec)))
  fixed = fixed[,dim(fixed)[2]:1,]
  fixed = aperm(fixed, c(2,1,3))
  fixed = brick(fixed)
  extent(fixed) = extent(random)
  
  fixed = writeRaster(fixed, 
                      filename = fixedFile,
                      overwrite = file.exists(fixedFile)
  )
  
  intensity = exp(fixed + random)
  intensity = setMinMax(intensity)
  
  intensity = writeRaster(intensity,
                          filename = intensityFile,
                          overwrite = file.exists(intensityFile)
  )
  
  summaryFunction = function(x, prefix = '',...) {
    Sprob=c(0.025, 0.5, 0.975)
    res = quantile(x, na.rm=TRUE, prob = Sprob)
    names(res) = paste('xx.',Sprob, "quant", sep='')
    c('xx.mean' = mean(x, na.rm=TRUE), res)
  }				
  
  randomSummary = calc(random, summaryFunction)
  names(randomSummary) = gsub("^xx", 'random', names(randomSummary))
  intensitySummary = calc(intensity, summaryFunction)
  names(intensitySummary) = gsub("^xx", 'intensity', names(intensitySummary))
  
  summaryStack = stack(
    randomSummary, intensitySummary
  )
  summaryStack = writeRaster(
    summaryStack, filename=summaryFile,
    overwrite = file.exists(summaryFile)
  )
  
  samples = list(
    intensity = intensity,
    random = random,
    sd = exp(bRes$etarec[,1]),
    range = 2*exp(bRes$etarec[,2])
  )
  
  colnames(bRes$betarec) = paste('beta', 
                                 seq(0,by=1, len=ncol(bRes$betarec)), sep=''
  )
  colnames(bRes$betarec) = gsub('^beta0$', '(Intercept)', colnames(bRes$betarec))
  
  summaryTable = lapply(
    c(
      as.list(as.data.frame(bRes$betarec)),		
      samples[c('range','sd')]
    ), summaryFunction)
  summaryTable = t(simplify2array(summaryTable))
  colnames(summaryTable) = gsub(
    '^xx\\.', '', colnames(summaryTable)			
  )
  
  xSeq = seq(0,
             exp(bRes$priors$etaprior$mean[1] +
                   4*sqrt(bRes$priors$etaprior$variance[1,1])),
             len=1000)
  
  sd = list(
    posterior = do.call(cbind, 
                        density(samples$sd)[c('x','y')]),
    prior = cbind(
      x = xSeq,
      y = dlnorm(xSeq, 
                 bRes$priors$etaprior$mean[1], 
                 bRes$priors$etaprior$variance[1,1]
      )
    )
  )
  
  
  xSeq = seq(0,
             exp(bRes$priors$etaprior$mean[2] +
                   4*sqrt(bRes$priors$etaprior$variance[2,2])),
             len=1000)
  
  range = list(
    posterior = do.call(cbind, 
                        density(samples$range)[c('x','y')]),
    prior = cbind(
      x = xSeq*2,
      y = dlnorm(xSeq, 
                 bRes$priors$etaprior$mean[2], 
                 bRes$priors$etaprior$variance[2,2]
      )/2
    )
  )
  ### Zhengfei Mar. 17, 2022
  samples = list(
    intensity = intensity,
    random = random,
    sd = exp(bRes$etarec[,1]),
    range = 2*exp(bRes$etarec[,2])
    , beta = bRes$betarec
  )
  
  list(
    summary = summaryTable,
    parameters = list(
      sd = sd,
      range = range),
    raster = summaryStack,
    samples = samples
  )
}		

#'
#'