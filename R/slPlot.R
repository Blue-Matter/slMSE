slplot = function(obj){

  if(class(obj)=="hist") slplot.hist(obj)
  if(class(obj)=="slSimData") slplot.data(obj)

}


dolab=function(i){mtext(paste0("(",letters[i],")"),adj=0.025, line=0.2,cex=0.7); return(i+1)}

slplot.hist = function(obj, scol = c("red","green","blue","black","darkgrey","purple"),labcex = 0.85){

  i = 1
  ylabline = 2.4; xlabline = 2.4
  par(mfrow=c(3,5),mai=c(0.5,0.5,0.2,0.05))
  snames = names(obj@OM@Stock)

  years = Years(obj, Period="Historical")
  # Life history

  len = sapply(obj@OM@Stock, function(x)x@Length@MeanAtAge[1,,1])
  ages = as.numeric(rownames(len))
  matplot(ages, len, type="l", ylim = c(0,max(len)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Length (cm)",2,line=ylabline, cex=labcex)
  legend('bottomright',legend=snames, text.col=scol,bty="n");
  i = dolab(i)


  wt = sapply(obj@OM@Stock, function(x)x@Weight@MeanAtAge[1,,1])
  matplot(ages, wt, type="l", ylim = c(0,max(wt)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Weight (kg)",2,line=ylabline, cex=labcex)
  i = dolab(i)

  mat = sapply(obj@OM@Stock, function(x)x@Maturity@MeanAtAge[1,,1])
  matplot(ages, mat, type="l", ylim = c(0,max(mat)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Spawning Fraction",2,line=ylabline, cex=labcex)
  i = dolab(i)

  #UEN = sapply(obj@Unfished@Equilibrium@Number,function(x) apply(x[,,1,drop=F],2,mean))
  #matplot(ages, UEN, type="l", ylim = c(0,max(UEN)), lty=1, col=scol,xlab="",ylab="");grid();
  #mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Unfished Eqbm. Num.",2,line=ylabline, cex=labcex)

  Mage = sapply(obj@OM@Stock,function(x)apply(x@NaturalMortality@MeanAtAge,2,mean))
  Surv = exp(-Mage*(1:nAge(obj)[[1]]))
  matplot(ages, Surv, type="l", ylim = c(0,max(Surv)), lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Survival (mean, across sims, years).",2,line=ylabline, cex=labcex)
  i = dolab(i)

  # Recruitment

  rec = sapply(obj@Number,function(x)apply(x[,1,,],2,mean))
  matplot(as.numeric(rownames(rec)), rec, type="l", lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Recruitment",2,line=ylabline, cex=labcex)
  i = dolab(i)

  arec = apply(array(rec,c(Seasons(obj),nYear(obj),nStock(obj))),2:3,sum)
  xval = CalcYears(nYear=nYear(obj), pYear = 0, CurrentYear=CurrentYear(obj))
  matplot(xval,arec, type="l", lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Annual Recruitment",2,line=ylabline, cex=labcex)
  i = dolab(i)

  rec_s = apply(array(rec,c(obj@OM@Seasons,nYear(obj),nStock(obj))),c(1,3),mean)
  matplot(rec_s, type="l", lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Season",1,line=xlabline, cex=labcex); mtext("Mean Recruitment",2,line=ylabline, cex=labcex)
  i = dolab(i)

  # Spatial Biomass

  Blist = get_spat_B(obj)
  ylim=c(0,max(unlist(lapply(Blist,max))))
  matplot(years,Blist[[1]],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nArea(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(years, Blist[[ss]],col=scol[ss],type='l',lty=1:nArea(obj),add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Biomass (kg)",2,line=ylabline, cex=labcex)
  legend('topright',legend=paste0("Area",1:nArea(obj)),lty=1:nArea(obj),bty="n")
  i = dolab(i)

  # Biomass

  Blist = get_spat_ann_B(obj)
  ylim=c(0,max(unlist(lapply(Blist,max))))
  xlab =  CalcYears(nYear=nYear(obj), pYear = 0, CurrentYear=CurrentYear(obj))
  matplot(xlab,Blist[[1]],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nArea(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(xlab, Blist[[ss]],col=scol[ss],type='l',lty=1:nArea(obj),add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Mean Annual Biomass (kg)",2,line=ylabline, cex=labcex)
  legend('topright',legend=paste0("Area",1:nArea(obj)),lty=1:nArea(obj),bty="n")
  i = dolab(i)

  # Ontogeny

  Bage = get_spat_age_B(obj)
  ylim=c(0,max(unlist(lapply(Bage,max))))
  matplot(ages,Bage[[1]],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nArea(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(ages, Bage[[ss]],col=scol[ss],type='l',lty=1:nArea(obj),add=T)
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext("Mean Biomass (kg)",2,line=ylabline, cex=labcex)
  legend('topright',legend=paste0("Area",1:nArea(obj)),lty=1:nArea(obj),bty="n")
  i = dolab(i)

  # Exploitation rate

  Fhist = apply(apply(obj@FDead,2:4,mean),2:1,sum)
  matplot(years, Fhist, type="l",lty=1, col=scol,xlab="",ylab="");grid();
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Total Apical Fish. Mort.",2,line=ylabline, cex=labcex)
  i = dolab(i)

  FbFs = apply(obj@FDead,2:4,mean)
  FbFy = apply(array(FbFs,c(nStock(obj),Seasons(obj),nYear(obj),nFleet(obj))),c(1,3,4),sum)
  ylim=c(0,max(FbFy))
  xval = CalcYears(nYear=nYear(obj), pYear = 0, CurrentYear=CurrentYear(obj))
  matplot(xval, FbFy[1,,],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nFleet(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(xval, FbFy[ss,,],col=scol[ss],type='l',lty=1:nFleet(obj),add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext("Annual Total Ap. Fish. Mort.",2,line=ylabline, cex=labcex)
  legend('top',legend=FleetNames(obj),lty=1:nFleet(obj),bty="n")
  i = dolab(i)

  FbF = apply(array(FbFs,c(nStock(obj),Seasons(obj),nYear(obj),nFleet(obj))),c(1,2,4),mean)
  ylim=c(0,max(FbF))
  matplot(FbF[1,,],col=scol[1],type='l',ylim=ylim,xlab="", ylab="", lty=1:nFleet(obj));grid()
  if(nStock(obj)>1)for(ss in 1:nStock(obj))matplot(FbF[ss,,],col=scol[ss],type='l',lty=1:nFleet(obj),add=T)
  mtext("Season",1,line=xlabline, cex=labcex); mtext("Mean Apical Fish. Mort.",2,line=ylabline, cex=labcex)
  legend('top',legend=FleetNames(obj),lty=1:nFleet(obj),bty="n")
  i = dolab(i)

  # Selectivity

  nts = Seasons(obj) * nYear(obj)
  # dim(hist@OM@Fleet$`Small Pacific JFS`$Long_highF@Selectivity@MeanAtAge)
  getfleetsel = function(x,ts) apply(x@Selectivity@MeanAtAge[,,ts,,drop=F],2,mean)
  getsfsel = function(y, ts) sapply(y,getfleetsel,ts=ts)

  sel_s = array(unlist(lapply(obj@OM@Fleet, getsfsel, ts=1)),c(nAge(obj)[[1]],nFleet(obj),nStock(obj)))
  matplot(ages,sel_s[,,1],type="l",col=scol[1],lty=1:nFleet(obj),xlab="", ylab="")
  if(nStock(obj)>1)for(ss in 2:nStock(obj))matplot(ages,sel_s[,,ss],type="l",col=scol[ss],lty=1:nFleet(obj),add=T);grid()
  legend('bottomright',legend=FleetNames(obj),lty=1:nFleet(obj),bty="n")
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext(paste0("Selectivity (",years[1],")"),2,line=ylabline, cex=labcex)
  i = dolab(i)

  sel_f = array(unlist(lapply(obj@OM@Fleet, getsfsel, ts=nts)),c(nAge(obj)[[1]],nFleet(obj),nStock(obj)))
  matplot(ages,sel_f[,,1],type="l",col=scol[1],lty=1:nFleet(obj),xlab="", ylab="")
  if(nStock(obj)>1)for(ss in 2:nStock(obj))matplot(ages,sel_f[,,ss],type="l",col=scol[ss],lty=1:nFleet(obj),add=T);grid()
  legend('bottomright',legend=FleetNames(obj),lty=1:nFleet(obj),bty="n")
  mtext("Age (years)",1,line=xlabline, cex=labcex); mtext(paste0("Selectivity (",round(years[nts],2),")"),2,line=ylabline, cex=labcex)
  i = dolab(i)


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
