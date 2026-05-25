slplot = function(obj){

  if(class(obj)=="om")   slplot.om(obj)
  if(class(obj)=="hist") slplot.hist(obj)
  if(class(obj)=="data") slplot.data(obj)

}


slplot.om = function(obj){


}



slplot.hist = function(obj, scol = c("red","green","blue","black","darkgrey","purple"),labcex = 0.85){

  ylabline = 2.4; xlabline = 2.6
  par(mfrow=c(3,4),mai=c(0.5,0.5,0.2,0.05))
  snames = names(obj@OM@Stock)

  # Life history

  len = sapply(obj@OM@Stock, function(x)x@Length@MeanAtAge[1,,1])
  matplot(as.numeric(rownames(len)), len, type="l", ylim = c(0,max(len)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Length (cm)",2,line=ylabline, cex=labcex)
  legend('bottomright',legend=snames, text.col=scol,bty="n")

  wt = sapply(obj@OM@Stock, function(x)x@Weight@MeanAtAge[1,,1])
  matplot(as.numeric(rownames(wt)), wt, type="l", ylim = c(0,max(wt)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Weight (kg)",2,line=ylabline, cex=labcex)

  mat = sapply(obj@OM@Stock, function(x)x@Maturity@MeanAtAge[1,,1])
  matplot(as.numeric(rownames(mat)), mat, type="l", ylim = c(0,max(mat)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Spawning Fraction",2,line=ylabline, cex=labcex)

  UEN = sapply(obj@Unfished@Equilibrium@Number,function(x) apply(x[,,1,drop=F],2,mean))
  matplot(as.numeric(rownames(UEN)), UEN, type="l", ylim = c(0,max(UEN)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Unfished Eqbm. Num.",2,line=ylabline, cex=labcex)

  # Recruitment

  rec = sapply(obj@Number,function(x)apply(x[,1,,],2,mean))
  matplot(as.numeric(rownames(rec)), rec, type="l", lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Recruitment",2,line=ylabline, cex=labcex)

  rec_s = apply(array(rec,c(obj@OM@Seasons,nYear(obj),nStock(obj))),c(1,3),mean)
  matplot(rec_s, type="l", lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Season",1,line=xlabline, cex=labcex); mtext("Mean Recruitment",2,line=ylabline, cex=labcex)

  # Ontogeny

  Bage = get_spat_age_B(obj)
  ylim=c(0,max(unlist(lapply(Bage,max))))
  xlab = 1:(nAge(obj)[[1]])
  matplot(xlab,Bage[[1]],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nArea(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(xlab, Bage[[ss]],col=scol[ss],type='l',lty=1:nArea(obj),add=T)
  mtext("Age classes",1,line=xlabline, cex=labcex); mtext("Mean Biomass by Area (kg)",2,line=ylabline, cex=labcex)
  legend('topright',legend=paste0("Area",1:nArea(obj)),lty=1:nArea(obj),bty="n")

  # Biomass

  Blist = get_spat_ann_B(obj)
  ylim=c(0,max(unlist(lapply(Blist,max))))
  xlab =  CalcYears(nYear=nYear(obj), pYear = 0, CurrentYear=CurrentYear(obj))
  matplot(xlab,Blist[[1]],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nArea(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(xlab, Blist[[ss]],col=scol[ss],type='l',lty=1:nArea(obj),add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Mean Biomass by Area (kg)",2,line=ylabline, cex=labcex)
  legend('topright',legend=paste0("Area",1:nArea(obj)),lty=1:nArea(obj),bty="n")

  # Spatial Biomass

  Blist = get_spat_B(obj)
  ylim=c(0,max(unlist(lapply(Blist,max))))
  xlab = Years(obj, Period="Historical")
  matplot(xlab,Blist[[1]],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nArea(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(xlab, Blist[[ss]],col=scol[ss],type='l',lty=1:nArea(obj),add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Biomass by Area (kg)",2,line=ylabline, cex=labcex)
  legend('topright',legend=paste0("Area",1:nArea(obj)),lty=1:nArea(obj),bty="n")

  # Exploitation rate

  FbFs = apply(obj@FDead,2:4,mean)
  FbF = apply(array(FbFs,c(nStock(obj),Seasons(obj),nYear(obj),nFleet(obj))),c(1,2,4),mean)
  ylim=c(0,max(FbF))
  matplot(FbF[1,,],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nFleet(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(FbF[ss,,],col=scol[ss],type='l',lty=1:nFleet(obj),add=T)
  mtext("Season",1,line=xlabline, cex=labcex); mtext("Mean Apical Fish. Mort.",2,line=ylabline, cex=labcex)
  legend('top',legend=FleetNames(obj),lty=1:nFleet(obj),bty="n")

  FbFy = apply(array(FbFs,c(nStock(obj),Seasons(obj),nYear(obj),nFleet(obj))),c(1,3,4),sum)
  ylim=c(0,max(FbFy))
  xval = CalcYears(nYear=nYear(obj), pYear = 0, CurrentYear=CurrentYear(obj))
  matplot(xval, FbFy[1,,],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nFleet(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(xval, FbFy[ss,,],col=scol[ss],type='l',lty=1:nFleet(obj),add=T)
  mtext("Season",1,line=xlabline, cex=labcex); mtext("Total Ann. Ap. Fish. Mort.",2,line=ylabline, cex=labcex)
  legend('top',legend=FleetNames(obj),lty=1:nFleet(obj),bty="n")

  Fhist = apply(apply(obj@FDead,2:4,mean),2:1,sum)
  xlab = Years(obj, Period="Historical")
  matplot(xlab, Fhist, type="l",lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Total Apical Fish. Mort.",2,line=ylabline, cex=labcex)

}

# internal functions


get_spat_B = function(hist){

  blist = list()
  for(ss in 1:nStock(hist)){
    N = apply(hist@Number[[ss]],2:4,mean)
    W = hist@OM@Stock[[ss]]@Weight@MeanAtAge[1,,]
    blist[[ss]] = apply(array(N*array(W,dim(N)),dim(N)),2:3,sum) # sum over ages
  }
  blist

}


get_spat_ann_B = function(hist){

  blist = list()
  for(ss in 1:nStock(hist)){
    N = apply(hist@Number[[ss]],2:4,mean)
    W = hist@OM@Stock[[ss]]@Weight@MeanAtAge[1,,]
    Btemp = apply(array(N*array(W,dim(N)),dim(N)),2:3,sum) # sum over ages
    blist[[ss]] = apply(array(Btemp,c(Seasons(hist),nYear(hist),nArea(hist))),c(2,3),mean)
  }
  blist

}

get_spat_age_B = function(hist){

  blist = list()
  for(ss in 1:nStock(hist)){
    N = apply(hist@Number[[ss]],2:4,mean)
    W = hist@OM@Stock[[ss]]@Weight@MeanAtAge[1,,]
    blist[[ss]] = apply(array(N*array(W,dim(N)),dim(N)),c(1,3),sum) # sum over ages
  }
  blist

}


