# mahiMSE plotting functions
# MSE = miniMSE

# MSE = Project(Hist, c("mahiMP","mahiMP2"))

# sims = 1:2; qsw = c(0.05,0.95); qsn = c(0.25,0.75); y0 = TRUE; ylim = NA; addn=T; addmed=T; addl = T; xlab = "Projection Year"; ylab = "Value"; fcol= "#10109950"
proj_plot = function(mat, lab, sims = 1:2, qsw = c(0.05,0.95), qsn = c(0.25,0.75), y0 = TRUE, ylim = NA, addn=T, addmed=T, addl = T,
                     xlab = "Projection Year", ylab = "Value",fcol= "#10109950", add = F){
  yw = apply(mat,2,quantile,p=qsw)
  yn = apply(mat,2,quantile,p=qsn)
  med = apply(mat,2,median)
  if(is.na(ylim[1])){
    ylim = range(yw)
    if(y0) ylim[1] = 0
  }
  if(!add)matplot(range(lab),ylim,col="white",xlab = xlab, ylab = ylab)
  grid()
  polygon(c(lab,rev(lab)),c(yw[1,],rev(yw[2,])),col=fcol, border=fcol)
  if(addn) polygon(c(lab,rev(lab)),c(yn[1,],rev(yn[2,])),col=fcol, border=fcol)
  if(addmed) lines(lab,med,col="white",lwd=2)
  if(addl) matplot(lab,t(mat[sims,]),col="black",lty=1:2,lwd=1,type="l",add=T)

}

# varnam = "SBiomass"; varlab = "Spawning Stock Biomass (t)"; MPs = 1:2; denom = 1E3; stock= 1; fcols = c("#10109950","#99000050"); ns = 4; annual = T
MP_comp = function(MSEobj, varnam = "SBiomass", varlab = "Spawning Stock Biomass (t)", MPs = 1:2, denom = 1E3, stock= 1, fcols = c("#0000ff50","#ff000050","#00ff0050"),ns, annual = T, addn=T, addmed=T, addl = T){

  pt = MSEobj@OM@pYear*ns; nt = MSEobj@OM@nYear; lhy = MSEobj@OM@Years[MSEobj@OM@nYear]; py = pt/ns; nsim = nSim(MSEobj); nmp = length(MSEobj@MPs)
  var0 = slot(MSEobj,varnam)/denom
  if(length(dim(var0))==4)var = array(var0[,stock,,,drop=F],dim(var0)[c(1,3,4)])
  if(length(dim(var0))==5)var = apply(var0[,stock,,,,drop=F],c(1,3,5),sum)

  if(annual){
    lab = lhy + 1:py
    newvar = array(var, c(nsim, ns, py, nmp))
    if(varnam == "SBiomass") var = apply(newvar,c(1,3,4),mean) # average SSB over year
    if(varnam == "Landings") var = apply(newvar,c(1,3,4),sum)  # sum of catches over year
  }else{
    lab = rep(lhy + 0:(py-1),each=ns) + seq(1/ns, 1, length.out = ns)
  }
  ylim = c(0,quantile(var,0.9975))
  posMPs = 1:length(MSEobj@MPs)
  MPs = MPs[MPs%in%posMPs]
  for(mp in MPs){
    mat = var[,,mp]
    proj_plot(mat,lab, ylab = varlab,add=mp!=MPs[1],fcol=fcols[match(mp,MPs)],
              ylim = ylim, addn = addn, addmed = addmed, addl = addl)
  }
  for(mp in MPs){
    mat = var[,,mp]
    matplot(lab,apply(mat,2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[match(mp,MPs)])
    matplot(lab,apply(mat,2,median),type="l",lwd=2,lwy=1,add=T,col=fcols[match(mp,MPs)])
  }

}

#' Plot overall yield and biomass projections across MPs
#'
#' @param MSEobj A class of object 'MSE'
#' @param MPs An integer vector of MPs to plot
#' @param annual Boolean, should annual or seasonal results be plotted?
#' @param ns Integer, number of seasons.
#' @param MPcols Character vector of colors for plotting MP projections
#' @param npy Positive integer, number of plot years - No. historical years to plot before projection
#' @param MPnam Character vector, optional, a character vector of MP names
#' @param leg Boolean, add a legend for the MPs?
#' @param addn Boolean, add narrow interquartile range shading?
#' @param addmed Boolean, add a line for the median?
#' @param addl Boolean, add individual simulation lines?
#' @examples
#' slplot.mse(myMSE)
#' @author T. Carruthers
#' @export
slplot.mse = function(MSEobj, MPs = 1:2, annual = T, ns = 4, stock = 1,  MPcols = c("#0000ff50","#ff000050","#00ff0050"),
                    npy=20, MPnams = NA,leg=T, addn=F, addmed=T, addl = F){
  # MSEobj =anMSE; MPs = 1:2; annual = T; ns = 4; stock = 1;  MPcols = c("red","green","blue","orange","grey","purple"); npy=20; MPnams = NA; leg=F; addn=F; addmed=T; addl = F
  par(mfrow=c(1,2),mai=c(0.25,0.85,0.15,0.05),omi=c(0.5,0.01,0.01,0.01))
  MP_comp(MSEobj, ns=ns, stock=stock, MPs = MPs,  addn = addn, addmed = addmed, addl = addl, annual=annual, fcols = MPcols)
  #abline(v=MSElist[[1]]@OM@CurrentYear+0.5,lty=2)
  MP_comp(MSEobj, varnam = "Landings", varlab = "Landings (t)",ns=ns,stock=stock, MPs=MPs, addn = addn, addmed = addmed, addl = addl, annual = annual,fcols=MPcols)
  #abline(v=MSElist[[1]]@OM@CurrentYear+0.5,lty=2)
  if(leg)legend('bottomright',legend = names(MSEobj@MPs)[MPs],text.col=MPcols,bty='n',text.font=2)
  mtext("Projection Year",1,outer=T,line=1.2)
}

