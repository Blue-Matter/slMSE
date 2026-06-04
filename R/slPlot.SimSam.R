
SS_pnt_bias = function(slt,keep, xlab,breaks=20,labcex = 0.8){
  ylabline = 2.2; xlabline = 2.2; labcex = 0.8
  err = ((slt$est[keep]-slt$sim[keep])/slt$sim[keep])*100
  xlim = range(err)
  if(xlim[1]>0)xlim[1] = 0
  if(xlim[2]<0)xlim[2] = 0
  hist(err,xlim = xlim, breaks=breaks, xlab="", ylab="",main="");grid()
  abline(v=0,lty=2,lwd=2)
  mtext("Freq. (nsim)",2,line=ylabline, cex=labcex); mtext(xlab,1,line=xlabline, cex=labcex)


}

# ylab = "test"
SS_ts_bias = function(slt, keep, ylab){

  err = ((slt$est[keep,]-slt$sim[keep,])/slt$sim[keep,])*100
  proj_plot(err, as.numeric(colnames(err)),ylab =ylab, y0=T)
  legend('topleft',legend=paste0("Mean bias = ",round(mean(err),2),"%"),bty='n')
  abline(h=0,lty=2,lwd=2)

}

# slt = obj$B; xlabs = as.numeric(colnames(obj$U$sim)); lab ="bio"; denom = 1E3; doleg=T
SS_ts_sim_bias = function(slt, xlabs, sims, ylab, denom, scol, doleg=F){
  ylabline = 2.2; xlabline = 2.2; labcex = 0.8
  ylim=c(0,max(slt$sim[sims,],slt$est[sims,],na.rm=T))/denom
  matplot(xlabs,t(slt$sim[sims,])/denom, ylim=ylim, type="l", lty=1, col=scol, lwd=2, xlab="", ylab=""); grid()
  matplot(xlabs,t(slt$est[sims,])/denom, ylim=ylim, type="l", lty=2, col=scol, lwd=1, xlab="", ylab="", add=T)
  mtext("Year",1,line=xlabline, cex=labcex); mtext(ylab,2,line=ylabline, cex=labcex)
  if(doleg) legend('topright',legend=paste0("sim #",sims),text.col=scol,bty='n',cex=0.9)
  if(doleg) legend('topleft',legend=c("Simulated","Estimated"),lty=c(1,2),lwd=c(2,1),bty='n',cex=0.9)

}


slplot.simsam = function(obj, nsimplot=3, scol = c("red","green","blue","black","darkgrey","purple"),labcex = 0.8,Bdenom = 1E3){

  nsim = length(obj$BMSY$sim)
  yrs = round(as.numeric(colnames(obj$U$sim)),2)
  conv = is.conv(obj)

  if(all(!conv))stop("! None of the estimation models converged !")
  sims = ((1:nsim)[conv])[1:nsimplot]

  i = 1
  ylabline = 2.2; xlabline = 2.2
  par(mai=c(0.55,0.55,0.2,0.05))
  layout(matrix(c(1,1,1,2,2,2,1,1,1,2,2,2,3,3,3,4,4,4,3,3,3,4,4,4,5,5,6,6,7,7,5,5,6,6,7,7),ncol=6,byrow=T))

  # Biomass
  SS_ts_sim_bias(obj$B, yrs, sims, "Biomass (B) (t)", 1E3, scol, T);   i = dolab(i)
  SS_ts_bias(obj$B, conv, "B.Est.Err. (E-S)/S %"); i = dolab(i)

  # Harvest rate
  SS_ts_sim_bias(obj$U, yrs, sims, "Harvest Rate (U)", 1, scol, T);  i = dolab(i)
  SS_ts_bias(obj$U, conv, "U.Est.Err. (E-S)/S %");  i = dolab(i)

  # MSY quantities
  SS_pnt_bias(obj$BMSY,conv,"BMSY Est.Err. (E-S)/S %");  i = dolab(i)
  SS_pnt_bias(obj$MSY,conv,"MSY Est.Err. (E-S)/S %");  i = dolab(i)
  SS_pnt_bias(obj$UMSY,conv,"UMSY Est.Err. (E-S)/S %");  i = dolab(i)
  if(any(!conv))cat(paste0("Assessment did not converge for simulations: ", paste((1:nsim)[!conv],collapse=", "),". ", sum(!conv)," simulations were removed \n"))

}

is.conv=function(obj){
  est = obj$B$est[,1]
  sim = obj$B$sim[,1]
  test1 = !is.na(est) & ((est-sim)/sim) < 10
  est = obj$BMSY$est
  sim = obj$BMSY$sim
  test2 = (est-sim)/sim < 10
  test1 & test2
}
