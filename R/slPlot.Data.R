
#obj = hist

tsdatplot = function(ts, ylab, ftcol, labcex = 0.85){

  nf = dim(ts)[3]
  nsim = dim(ts)[1]
  yrs = as.numeric(dimnames(ts)[[2]])
  ylim = c(0,max(ts)*1.05)
  plot(range(yrs),ylim,col="white",xlab="",ylab=""); grid()
  for(ff in 1:nf)matplot(yrs,t(ts[,,ff]),type="l",col=ftcol[ff],lty=1:nsim,lwd=1,add=T)
  mtext(ylab,2,line=2.4,cex = labcex)
  mtext("Year",1,line=2.4,cex=labcex)
}

densplot = function(xs,mat,sims, xlab,ftcol,labcex = 0.85){
  matplot(range(xs),c(0,max(mat)*1.05),col="white",xlab="",ylab=""); grid()
  for(sim in sims)matplot(xs,t(mat[sim,,]),type="l",col=ftcol,lty=sim, add=T)
  mtext("Frequency",2,line=2.4,cex = labcex)
  mtext(xlab,1,line=2.4,cex=labcex)
}

slplot.data = function(obj, sims = 1:2, ftcol = c("red","blue","green","black","darkgrey","purple"),Bdenom=1E3){

  i = 1
  par(mfrow=c(3,4),mai=c(0.5,0.5,0.2,0.05))
  Fleetnams = dimnames(obj$Landings)[[3]]

  # Total annual landings
  Ty = array(apply(obj$Landings_y, 1:2,sum),c(dim(obj$Landings_y)[1:2],1))
  dimnames(Ty)[[2]] = dimnames(obj$Landings_y)[[2]]
  tsdatplot(ts=Ty[sims,,,drop=F]/Bdenom,"Total Annual Landings (t)","black")
  legend('top', legend=paste0("Sim ",sims),lty=sims,col ="black", bty='n')
  i = dolab(i)

  # Landings by fleet
  tsdatplot(obj$Landings[sims,,,drop=F]/Bdenom,"Landings by Fleet, Season (t)",ftcol)
  legend('top', legend=Fleetnams, text.col=ftcol, bty='n')
  i = dolab(i)

  # Annual landings by fleet
  tsdatplot(obj$Landings_y[sims,,,drop=F]/Bdenom,"Annual Landings by Fleet (t)",ftcol)
  i = dolab(i)

  # CPUE by fleet
  tsdatplot(obj$CPUE[sims,,,drop=F],"CPUE by Fleet, Season",ftcol)
  i = dolab(i)

  # Annual CPUE by fleet
  tsdatplot(ts = obj$CPUE_y[sims,,,drop=F],"Annual CPUE by Fleet",ftcol)
  i = dolab(i)

  # CPUE by fleet
  tsdatplot(obj$Survey[sims,,,drop=F],"Survey by Fleet, Season",ftcol)
  i = dolab(i)

  # Annual CPUE by fleet
  tsdatplot(ts = obj$Survey_y[sims,,,drop=F],"Annual Survey by Fleet",ftcol)
  i = dolab(i)

  # Catch length composition
  CAL_y = totlen = obj$CAL_y
  cind = MSEtool:::TEG(dim(CAL_y))
  lens = as.numeric(dimnames(CAL_y)[[5]])
  totlen[cind] = CAL_y[cind]*lens[cind[,5]]
  mulen = apply(totlen,c(1,4,3),sum) / apply(CAL_y,c(1,4,3),sum)
  # Annual CPUE by fleet
  tsdatplot(ts = mulen[sims,,,drop=F],"Mean Length in Catch (cm)",ftcol)
  i = dolab(i)

  CAL_s = apply(CAL_y[,,,1,],c(1,3,4),sum)
  CAL_f = apply(CAL_y[,,,dim(CAL_y)[4],],c(1,3,4),sum)
  yrs = as.numeric(dimnames(CAL_y)[[4]])
  densplot(lens, CAL_s, sims, paste0("Length (cm) (",yrs[1],")"), ftcol)
  i = dolab(i)
  densplot(lens, CAL_f, sims, paste0("Length (cm) (",yrs[length(yrs)],")"), ftcol)
  i = dolab(i)

}
