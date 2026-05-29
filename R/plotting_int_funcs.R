

# panel labelling plus counter
dolab=function(i){mtext(paste0("(",letters[i],")"),adj=0.025, line=0.2,cex=0.7); return(i+1)}

# generic projection plot
proj_plot = function(mat, lab, sims = 1:2, qsw = c(0.05,0.95), qsn = c(0.25,0.75), y0 = TRUE, ylim = NA, addn=T, addmed=T, addl = F,
                     xlab = "Year", ylab = "Value",fcol= "#10109950", add = F){
  ylabline = 2.2; xlabline = 2.2; labcex=0.85
  yw = apply(mat,2,quantile,p=qsw)
  yn = apply(mat,2,quantile,p=qsn)
  med = apply(mat,2,median)
  if(is.na(ylim[1])){
    ylim = range(yw)
    if(y0){
      if(ylim[1]>0)ylim[1]= 0
      if(ylim[2]<0)ylim[2] = 0
    }
  }
  if(!add)matplot(range(lab),ylim,col="white",xlab = "", ylab = "")
  mtext(xlab,1,line=xlabline, cex=labcex); mtext(ylab,2,line=ylabline, cex=labcex)


  grid()
  polygon(c(lab,rev(lab)),c(yw[1,],rev(yw[2,])),col=fcol, border=fcol)
  if(addn) polygon(c(lab,rev(lab)),c(yn[1,],rev(yn[2,])),col=fcol, border=fcol)
  if(addmed) lines(lab,med,col="white",lwd=2)
  if(addl) matplot(lab,t(mat[sims,]),col="black",lty=1:2,lwd=1,type="l",add=T)

}
