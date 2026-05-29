
# Acknowledgement: This code was heavily influenced by code developed by Dr Henning Winker (legend that he is!)
# r.pr=c(0.4583,  0.2553, 1); bk.pr=c(0.5,0.3,1); shape.pr=c(2, 0.001, 1); oe=c(0.3, 0.5, 1); pe= c(0.2,0.5,1); fdevs=c(4, 0.5, 1); ce=c(0.05, 0.001, 1); q.pr= c(1,0.5,1); timing=0.625; dteuler=0.25; n_Indices = 1

SPiCT_config = function(Sdata, r.pr=c(0.6,  0.2, 1), bk.pr=c(0.5,0.3,1),
                        shape.pr=c(2, 0.001, 1), oe=c(0.3, 0.5, 1), pe= c(0.4,0.5,1),
                        fdevs=c(4, 0.5, 1), ce=c(0.05, 0.001, 1), q.pr= c(1,0.5,1),
                        timing=0.625, dteuler=0.25){

  config = Sdata
  config$priors = list()
  n_Ind = length(Sdata$obsI)

  config$dteuler = dteuler           # Resolution of continuous SP model calcs
  config$timeI = Map(function(x,y){x+y}, y=as.list(timing),x=config$timeI) # Survey timing
  config$priors$logalpha <- c(0,0,0) # Deactivate ratio of proc to obs
  config$priors$logbeta <- c(0,0,0)  # Deactivate catch to f devs

  if(!is.null(bk.pr))config$priors$logbkfrac <- c(log(bk.pr[1])-0.5*bk.pr[2]^2/2,bk.pr[2],bk.pr[3]) # SET bk prior'

  if(!is.null(oe)){ # Observation error
    config$priors$logsdi = list()
    for(i in 1:length(config$obsI)){
      config$priors$logsdi[[i]] = c(log(oe[1]), oe[2], oe[3])
    }
  }

  if(!is.null(pe)) config$priors$logsdb <- c(log(pe[1])-0.5*pe[2]^2, pe[2], pe[3]) # Process error

  if(!is.null(fdevs)) config$priors$logsdf <- c(log(fdevs[1])-0.5*fdevs[2]^2, fdevs[[2]], fdevs[3]) # F devs: should be higher

  if(!is.null(ce)) config$priors$logsdc <- c(log(ce)[1]-0.5*ce[2]^2, ce[2], ce[3]) # Catch error

  if(!is.null(shape.pr)) config$priors$logn <- c(log(shape.pr[1]),shape.pr[2],shape.pr[3]) # Shape CV

  if(!is.null(r.pr))config$priors$logr <- c(log(r.pr[1])-0.5*r.pr[2]^2,r.pr[2],r.pr[3]) # r prior

  if(!is.null(q.pr)){             # Catchability prior
    config$priors$logq = list()
    for(i in 1:n_Ind){
      config$priors$logq[[i]] = c(log(q.pr[1])-0.5*q.pr[2]^2,q.pr[2],q.pr[3])
    }
  }

  return(config)

}


SPiCT_data = function(sim, simdata, time_step = "year"){

  if(time_step == "year"){

    # Catches
    obsC = apply(simdata$Landings_y[sim,,,drop=F],2,sum) # sum over fleets
    timeC <- as.numeric(names(obsC))

    # Indices
    obsI = timeI = list();
    for(ff in 1:simdata$nFleet){
      obsI[[ff]] = simdata$Survey_y[sim,,ff]
      timeI[[ff]] = as.numeric(names(obsI[[ff]]))
    }

    # Data input list
    outlist = list(obsC=obsC, timeC=timeC, obsI=obsI, timeI=timeI)

  }else if(time_step == "quarter"){
    # Catches
    obsC = apply(simdata$Landings[sim,,,drop=F],2,sum) # sum over fleets
    timeC <- as.numeric(names(obsC))

    # Indices
    obsI = timeI = list();
    for(ff in 1:simdata$nFleet){
      obsI[[ff]] = simdata$Survey[sim,,ff]
      timeI[[ff]] = as.numeric(names(obsI[[ff]]))
    }

    # Data input list
    outlist = list(obsC=obsC, timeC=timeC, obsI=obsI, timeI=timeI)
  }

  outlist
}


getsval = function(fit, nam, ts=T){
  temp = fit$value[names(fit$value)==nam]
  if(ts){
    if(length(temp) == length(fit$inp$time)){
      keep = fit$inp$time %in% fit$inp$timeC
      temp = temp[keep]
      names(temp) = fit$inp$timeC
    }else{
      keep = 1:length(fit$inp$timeC)
      temp = temp[keep]
      names(temp) = fit$inp$timeC
    }
  }
  temp
}

SPiCT_output = function(fit){
  predts = fit$inp$time
  BMSY = getsval(fit, "Bmsy",ts=F)
  MSY = getsval(fit, "MSY",ts=F)
  BMSY_B0 = exp(getsval(fit, "logBpK",ts=F))
  Brel = exp(getsval(fit, nam = "logBBmsy"))
  B = BMSY * Brel
  C = exp(getsval(fit, nam = "logCpred"))
  U = C/B
  list(B = B, C = C, U = U, BMSY = BMSY, MSY = MSY, BMSY_B0 = BMSY_B0, predts = predts)
}


#  sim = 2; timestep="quarter"; r.pr = c(0.8,0.1,1); bk.pr = c(0.5,0.3,1);shape.pr = c(2, 0.001, 1);oe = c(0.5, 0.5, 1);pe = c(0.8,0.5,1);fdevs = c(4, 0.5, 1);ce = c(0.05, 0.001, 1);q.pr = NULL;timing = 0;dteuler = 0.01
do_spict=function(sim, simdata, timestep = 'year',
                  r.pr = c(0.6,0.2,1),
                  bk.pr = c(0.5,0.3,1),
                  shape.pr = c(2, 0.001, 1),
                  oe = c(0.2, 0.5, 1),
                  pe = c(0.5,0.5,1),
                  fdevs = c(4, 0.5, 1),
                  ce = c(0.05, 0.001, 1),
                  q.pr = NULL,
                  timing = 0.01,
                  dteuler = 0.25){

  Sdata = SPiCT_data(sim, simdata, timestep)
  Sinput = SPiCT_config(Sdata, r.pr, bk.pr, shape.pr, oe, pe, fdevs, ce,  q.pr, timing, dteuler)
  check.inp(Sinput)
  fit = fit.spict(Sinput)
  SPiCT_output(fit)

}


#  timestep = "year"; parallel = T; r.pr = c(0.5,0.2,1); bk.pr = c(0.5,0.3,1);shape.pr = c(2, 0.001, 1);oe = c(0.2, 0.5, 1);pe = c(0.2,0.5,1);fdevs = c(4, 0.5, 1);ce = c(0.05, 0.001, 1);q.pr = NULL;timing = 0.625;dteuler = 0.25


SimSam_spict = function(simdata, timestep = "year", parallel =T,
                        r.pr = c(0.5,0.2,1),
                        bk.pr = c(0.5,0.3,1),
                        shape.pr = c(2, 0.001, 1),
                        oe = c(0.2, 0.5, 1),
                        pe = c(0.2,0.5,1),
                        fdevs = c(4, 0.5, 1),
                        ce = c(0.05, 0.001, 1),
                        q.pr = NULL,
                        timing = 0.625,
                        dteuler = 0.25){

  if(parallel){
    setup()
    sfLibrary(slMSE)
    sfLibrary(spict)
    Est = sfLapply(1:nSim,do_spict, simdata = simdata, timestep = timestep, r.pr = r.pr, bk.pr = bk.pr,
                   shape.pr = shape.pr, oe = oe, pe = pe, fdevs = fdevs, ce = ce,
                   q.pr = q.pr, timing = timing, dteuler = dteuler)
  }else{
    Est = lapply(1:nSim,do_spict, simdata = simdata, timestep = timestep, r.pr = r.pr, bk.pr = bk.pr,
                 shape.pr = shape.pr, oe = oe, pe = pe, fdevs = fdevs, ce = ce,
                 q.pr = q.pr, timing = timing, dteuler = dteuler)
  }


  SimSam = list()
  estU = t(sapply(Est,function(x)x$U))
  estB = t(sapply(Est,function(x)x$B))
  estBMSY = sapply(Est,function(x)x$BMSY)
  estMSY = sapply(Est,function(x)x$MSY)
  estUMSY = estMSY / estBMSY
  SimSam$U = list(sim = simdata$U, est = estU)
  SimSam$B = list(sim = simdata$B, est = estB)
  SimSam$BMSY = list(sim = simdata$BMSY, est = estBMSY)
  SimSam$MSY = list(sim = simdata$MSY, est = estMSY)
  SimSam$UMSY = list(sim = simdata$UMSY, est = estUMSY)
  class(SimSam) = "slSimSam"
  SimSam
}


# advice generators
# get.TAC(fit, fractiles = list(catch = 0.35), breakpointB = 0.5, limitB = 0.3).
# fit <- retro(fit)
# fit <- check.ini(fit)). The
# estimates should be the same for all initial values (fit$check.ini$resmat
# inp$optimiser.control = list(iter.max = 1e5, eval.max = 1e5) # increase numerical solve iterations
# inp$ini$logn <- log(2); inp$phases$logn <- -1 # fix to schaefer
# intiial biomass close to carrying capacity inp$priors$logbkfrac <- c(log(0.8),0.5,1)
# fit <- manage(fit, scenarios = "ices").

