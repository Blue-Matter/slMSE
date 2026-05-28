
# r.pr=c(0.4583,  0.2553, 1); bk.pr=c(0.5,0.3,1); shape.pr=c(2, 0.001, 1); oe=c(0.3, 0.5, 1); pe= c(0.2,0.5,1); fdevs=c(4, 0.5, 1); ce=c(0.05, 0.001, 1); q.pr= c(1,0.5,1); timing=0.625; dteuler=0.25; n_Indices = 1

SPiCT_config = function(Sdata, r.pr=c(0.4,  0.2, 1), bk.pr=c(0.5,0.3,1),
                        shape.pr=c(2, 0.001, 1), oe=c(0.3, 0.5, 1), pe= c(0.2,0.5,1),
                        fdevs=c(4, 0.5, 1), ce=c(0.05, 0.001, 1),q.pr= c(1,0.5,1),
                        timing=0.625, dteuler=0.25){

  # vvv ripped from Henning Winker Code vvvvv
  config = Sdata
  config$priors = list()
  n_Ind = length(Sdata$obsI)

  config$dteuler = dteuler           #
  config$timeI = Map(function(x,y){x+y}, y=as.list(timing),x=config$timeI) #A SSIGN survey timing
  config$priors$logalpha <- c(0,0,0) # DEACTIVATE ratio of proc to obs
  config$priors$logbeta <- c(0,0,0)  # DEACTIVATE catch to f devs

  if(!is.null(bk.pr))config$priors$logbkfrac <- c(log(bk.pr[1])-0.5*bk.pr[2]^2/2,bk.pr[2],bk.pr[3]) # SET bk prior'

  if(!is.null(oe)){ # SET obs error
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


SPiCT_data = function(sim, simdata){

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
  list(obsC=obsC, timeC=timeC, obsI=obsI, timeI=timeI)

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
  Brel = exp(getsval(fit, nam = "logBBmsy"))
  B = BMSY * Brel
  C = exp(getsval(fit, nam = "logCpred"))
  U = C/B
  list(B = B, C = C, U = U, BMSY = BMSY, predts = predts)
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

# vvvv --------------- Henning's Code: -------------- vvvvv

# --- spict

# buildTMBinputSPICT {{{

# obsC, timeC, obsI, timeI

buildTMBinputSPICT <- function(catch, indices) {

  # PARSE and CHECK arguments
  ny <- dim(catch)[2]

  # C_hist

  dataC <- as.data.frame(catch)
  obsC <- dataC$data
  timeC <- dataC$year

  # obsI

  dataI <- lapply(indices, function(x) as.data.frame(x))
  obsI <- lapply(dataI, '[[', 'data')
  timeI <- lapply(dataI, '[[', 'year')

  # BUILD input list
  inp <- list(obsC=obsC, timeC=timeC, obsI=obsI, timeI=timeI)

  return(inp)
}
# }}}

# spict.sa {{{

spict.sa <- function(stk, idx, args, tracking,fractile=0.35,qtype=c("B","Fmsy")[1],indicator=c("FLQuant","FLStock")[1]) {

  # EXTRACT inputs
  catch <- catch(stk) # + min(c(catch(stk))[c(catch(stk)) > 0]) * 1
  indices <- lapply(idx, index)

  # OUTPUT
  empty <- catch %=% 0
  out <- vector(mode="list", length=args$it)
  # First to converge
  CONV = FALSE
  conv = rep(0,args$it) # track
  # LOOP over iters
  for (i in seq(args$it)) {

    cat(".")

    # EXTRACT single iter
    cab <- iter(catch, i)
    index <- iter(indices, i)

    inp <- buildTMBinputSPICT(cab, index)

    # SET simple options
    inp <- spict.settings(inp)

    # FIT spict
    # BUG: FAILS at nlminb, so cannot catch
    fit <- tryCatch(fit.spict(inp),
                    # error, RETURN 1 convergence + 0 output
                    error = function(e) {
                      print(i)
                      return(list(
                        opt=list(convergence=3),
                        report=list(MSY=0, Fmsy=0, Bmsy=0),
                        value=c(K=0)))
                    },
                    warning = function(w) {
                      print(i)
                    }
    )


    # convergence flag
    if(!class(fit)=="spictcls"){
      conv[i] <- 1
    } else {


      conv[i] <- fit$opt$convergence



      # EXTRACT output indicators
      res <- FLRef::spict2FLStockR(fit,osa=TRUE)
      # Make B/Bmsy here for HCR metric
      res@stock = ssb(res)/res@refpts["Bmsy"]
      res@name = paste(i)
      # Get quantile of Fmsy

      # Fractile Fmsy
      Fqmsy=exp(qnorm(fractile, get.par("logFmsy", fit)[2], get.par("logFmsy", fit)[4]))
      # Fractile B
      endyr = an(range(res)["maxyear"])
      logB =  get.par("logB", fit)
      logB = logB[which(rownames(logB)==endyr+1)-1,]
      Bq=exp(qnorm(fractile, logB[2], logB[4]))


      if(qtype=="B") Fmsy.star = an(res@refpts[[1]]* Bq/ssb(res[,ac(endyr)]))
      if(qtype=="Fmsy") Fmsy.star = Fqmsy
      res@refpts <- rbind(res@refpts,FLPar(Fmsy.star=Fmsy.star))

      if(!CONV){
        CONV=TRUE # set flag to TRUE
        ind = propagate(res,args$it)
        refs = propagate(res@refpts,args$it)
      }


      iter(ind,i) <- res # Save results
      iter(refs,i)<-res@refpts
    } # converged models
  }

  # TRACK convergence
  #track(tracking, "conv.est", ac(args$ay)) <- unlist(lapply(out, '[[', 'conv'))
  track(tracking, "conv.est", ac(args$ay)) <- conv

  #browser()
  if(indicator=="FLQuant"){
    ind=metrics(ind,metrics=list(Bratio=stock,Biomass=ssb))
    attr(ind,"refpts") = refs

    track(tracking, "Bratio", ac(args$ay)) <- an(ind$Bratio[,ac(args$ay-1)])
    track(tracking, "Bspm", ac(args$ay-1)) <- an(ind$Biomass[,ac(args$ay-1)])
    #track(tracking, "Fspm", ac(args$ay-1)) <- an(ind$F[,ac(args$ay-1)])
  }



  # TODO: IAGO to review
  #ind <- tryCatch(combine(FLStocks(lapply(out,function(x){x$ind}))),error=function(e){})
  #ind <- tryCatch(FLStocks(lapply(setNames(nm=names(inds[[1]])), function(i)
  #  window(Reduce(combine, lapply(inds, "[[", i)), end=args$ay+1))),
  #  error=function(e) {
  #  })

  # USE ind directly as FLStockR
  return(list(stk=stk,ind=ind,refpts=refs,tracking=tracking))
  #return(list(stk=stk, ind=ind, tracking=tracking,rps))
}
# }}}


spict.settings = function(inp, r.pr=c(0.4583,  0.2553, 1), bk.pr=c(0.5,0.3,1),
                          shape.pr=c(2, 0.001, 1), oe=c(0.3, 0.5, 1), pe= c(0.2,0.5,1),
                          fdevs=c(4, 0.5, 1), ce=c(0.05, 0.001, 1),q.pr= c(1,0.5,1),timing=0.625, dteuler=0.25) {

  # ASSIGN arguments to input
  inp$dteuler=dteuler

  #ASSIGN survey timing
  inp$timeI = Map(function(x,y){x+y}, y=as.list(timing),x=inp$timeI)

  # DEACTIVATE ratio of proc to obs
  inp$priors$logalpha <- c(0,0,0)

  # DEACTIVATE catch to f devs
  inp$priors$logbeta <- c(0,0,0)

  # SET bk prior'
  if(!is.null(bk.pr))
    inp$priors$logbkfrac <- c(log(bk.pr[1])-0.5*bk.pr[2]^2/2,bk.pr[2],bk.pr[3])

  # SET obs error
  if(!is.null(oe)){
    inp$priors$logsdi = list()
    for(i in 1:length(inp$obsI)){
      inp$priors$logsdi[[i]] = c(log(oe[1]), oe[2], oe[3])
    }
  }
  # SET process error
  if(!is.null(pe))
    inp$priors$logsdb <- c(log(pe[1])-0.5*pe[2]^2, pe[2], pe[3])

  # SET F devs: should be higher
  if(!is.null(fdevs))
    inp$priors$logsdf <- c(log(fdevs[1])-0.5*fdevs[2]^2, fdevs[[2]], fdevs[3])

  # SET catch error
  if(!is.null(ce))
    inp$priors$logsdc <- c(log(ce)[1]-0.5*ce[2]^2, ce[2], ce[3])

  # REDUCE shape CV
  if(!is.null(shape.pr))
    inp$priors$logn <- c(log(shape.pr[1]),shape.pr[2],shape.pr[3])

  # SET r prior
  if(!is.null(r.pr))
    inp$priors$logr <- c(log(r.pr[1])-0.5*r.pr[2]^2,r.pr[2],r.pr[3])

  if(!is.null(q.pr)){
    inp$priors$logq = list()
    for(i in 1:length(inp$obsI)){
      inp$priors$logq[[i]] = c(log(q.pr[1])-0.5*q.pr[2]^2,q.pr[2],q.pr[3])
    }}


  return(inp)
}

