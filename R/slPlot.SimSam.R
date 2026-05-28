
SS_pnt_bias = function(slt,xlab,breaks=20,labcex = 0.8){
  ylabline = 2.2; xlabline = 2.2; labcex = 0.8
  err = ((slt$est-slt$sim)/slt$sim)*100
  xlim = range(err)
  if(xlim[1]>0)xlim[1] = 0
  if(xlim[2]<0)xlim[2] = 0
  hist(err,xlim = xlim, breaks=breaks, xlab="", ylab="",main="");grid()
  abline(v=0,lty=2,lwd=2)
  mtext("Freq. (nsim)",2,line=ylabline, cex=labcex); mtext(xlab,1,line=xlabline, cex=labcex)


}


SS_ts_bias = function(slt, ylab){

  err = ((slt$est-slt$sim)/slt$sim)*100
  proj_plot(err, as.numeric(colnames(err)),ylab =ylab, y0=T)
  abline(h=0,lty=2,lwd=2)

}

# slt = obj$B; xlabs = as.numeric(colnames(obj$U$est)); lab ="bio"; denom = 1E3; doleg=T
SS_ts_sim_bias = function(slt, xlabs, sims, ylab, denom, scol, doleg=F){
  ylabline = 2.2; xlabline = 2.2; labcex = 0.8
  ylim=c(0,max(slt$sim[sims,],slt$est[sims,]))/denom
  matplot(xlabs,t(slt$sim[sims,])/denom, ylim=ylim, type="l", lty=1, col=scol, lwd=2, xlab="", ylab=""); grid()
  matplot(xlabs,t(slt$est[sims,])/denom, ylim=ylim, type="l", lty=2, col=scol, lwd=1, xlab="", ylab="", add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext(ylab,2,line=ylabline, cex=labcex)
  if(doleg) legend('topright',legend=paste0("sim #",sims),text.col=scol,bty='n',cex=0.9)
  if(doleg) legend('topleft',legend=c("Simulated","Estimated"),lty=c(1,2),lwd=c(2,1),bty='n',cex=0.9)

}


slplot.simsam = function(obj, sims = 1:3, scol = c("red","green","blue","black","darkgrey","purple"),labcex = 0.8,Bdenom = 1E3){

  i = 1
  ylabline = 2.2; xlabline = 2.2
  par(mfrow=c(3,3),mai=c(0.5,0.5,0.2,0.05))
  yrs = as.numeric(colnames(obj$U$est))

  # Biomass
  SS_ts_sim_bias(obj$B, yrs, sims, "Biomass (B) (t)", 1E3, scol, T);   i = dolab(i)
  SS_ts_bias(obj$B, "B.Est.Err. (E-S)/S %"); i = dolab(i)

  # Harvest rate
  SS_ts_sim_bias(obj$U, yrs, sims, "Harvest Rate (U)", 1, scol, T);  i = dolab(i)
  SS_ts_bias(obj$U, "U.Est.Err. (E-S)/S %");  i = dolab(i)

  SS_pnt_bias(obj$BMSY,"BMSY");  i = dolab(i)
  SS_pnt_bias(obj$MSY,"MSY");  i = dolab(i)
  SS_pnt_bias(obj$UMSY,"UMSY");  i = dolab(i)


}
