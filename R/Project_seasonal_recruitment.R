
#' Takes an OpenMSE v1.0 object from RCM and imposes seasonality including estimated recruitment pattern in an OpenMSE v2.0 OM
#'
#' @param OM An object of openMSE class 'OM' (OpenMSE v1.0)
#' @param Seasons Positive integer. The number of seasons within year. Defaults to 4, quarterly.
#' @param trim  Positive integer. The number of recent historical years of estimated recruitment to ignore (due to poor characterization in the lower triangle).
#' @param RecLev Improper fraction. Overall mean recruitment strength relative to historical. Defaults to 1. Can be modified for purposes of OM specification.
#' @param sy Positive integer. The number of initial recruitment years (including age classes in initial year) to impose recruitment. Defaults to 0.
#' @param plot Boolean. Do you want to plot the estmated  / imposed recruitments?
#' @examples
#' Seasonality(RCMfit@OM)
#' @author T. Carruthers
#' @export

Seasonality = function(OM, Seasons = 4, trim = 4, RecLev=1, sy = 1, plot=F){
  # OM = fit@OM; Seasons = 4; trim = 4; RecLev=1; sy = 1; plot=T

  OM2 = ConvertOM(OM, Seasons = 4)        # Convert to seasonal model
  OM3 = make_rec(OM2, RecLev, sy, trim, plot)
  OM3@Interval = 1
  OM3

}

# truncated recruitment sampler (inefficient)
rec_samp=function(n,mu,cv,UB = 6){
  x = trlnorm(n,mu,cv)
  while(any(x>UB)){
    ind = x > UB
    nux = trlnorm(n,mu,cv)
    x[ind]=nux[ind]
  }
  x
}

make_rec= function(OM, RecLev=1, sy = 40, trim = 4, plot=T){ # rec sampled from year 11 onwards up to last four quarters
  # OM = OM2;

  RDI = OM@Stock[[1]]@SRR@RecDevInit
  RDP = OM@Stock[[1]]@SRR@RecDevProj
  RDH = OM@Stock[[1]]@SRR@RecDevHist
  Perr_y = nuPerr = cbind(RDI, RDH, RDP)
  na = nAge(OM,1)[[1]];
  ny = ncol(RDH); py = ncol(RDP); nsim = nSim(OM)
  qs = rep(1:4,1000)[1:ny]

  # Initial years
  #RecYrs = na+1:sy
  #for(i in 1:nsim){
  #  qmu = aggregate(Perr_y[i,RecYrs],by=list(q = qs[RecYrs]),mean)
  # qsd = aggregate(Perr_y[i,RecYrs],by=list(q = qs[RecYrs]),function(x){sd(x)/mean(x)})
  #nuPerr[i,1:na] = rec_samp(na,qmu$x,qsd$x)
  #}
  #nuPerr[,1:na] = 1

  # Projection years
  RecYrs = (na+sy):(ny-trim)
  pord = c((4-trim+1):4,1:(4-trim))[1:4]
  for(i in 1:nsim){
    qmu = aggregate(Perr_y[i,RecYrs],by=list(q = qs[RecYrs]),mean) # debug ; qmu$x = c(0.1,0.5,1,1.5)
    qsd = aggregate(Perr_y[i,RecYrs],by=list(q = qs[RecYrs]),function(x){sd(x)/mean(x)}) # debug ;qsd$x = rep(0.01,4)
    toind = (na+ny-trim):(na+ny+py-1)
    nuPerr[i,toind] = rec_samp(py+trim,qmu$x[pord]*RecLev,qsd$x[pord])
  }

  dat = data.frame(rd = nuPerr[2,toind],q = rep(pord,1000)[1:length(toind)])
  check = cbind(qmu,aggregate(dat$rd,by=list(dat$q),mean))
  names(check) = c("spec q","spec mu","model q","simulated mu")
  print(check)

  if(plot){
    par(mfrow=c(2,2),omi=c(0.4,0.4,0.01,0.01))
    for(i in 1:4){
      plot(nuPerr[i,],col=c("red","green","blue","grey"),pch=19)
      points(Perr_y[i,1:na+ny],,col=c("red","green","blue","grey"),pch=4,lwd=2)
      abline(v=na+ny)
    }
  }

  OM@Stock[[1]]@SRR@RecDevHist[,ny-((trim-1):0)] = nuPerr[,(ny+na-trim):(ny+na-1)]
  OM@Stock[[1]]@SRR@RecDevProj = nuPerr[,(ny+na):(py+ny+na-1)]

  if(plot){
    par(mfrow=c(2,2),mai=c(0.4,0.4,0.05,0.05),omi=c(0.3,0.3,0.01,0.01))
    for(i in 1:4){
      xlab = Years(OM)
      recs = c(OM@Stock[[1]]@SRR@RecDevHist[i,], OM@Stock[[1]]@SRR@RecDevProj[i,])
      plot(xlab,recs,
           ylab=paste0("Recdevs simulation ",i),col=c("red","green","blue","grey"),xlab="",,pch=19)
      grid()
      #points(Perr_y[i,1:na+ny],,col=c("red","green","blue","grey"),pch=4,lwd=2)
      abline(v=OM@CurrentYear+0.5)
      if(i == 1) legend('topleft',c("Winter","Spring","Summer","Autumn"),text.col=c("red","green","blue","grey"),bty="n")
    }
    mtext("Recruitment Deviation",2,outer=T)
    mtext("Year", 1, outer=T)
  }

  OM

}
