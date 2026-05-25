
#obj = hist

tsdatplot_byft = function(ts, ylab, ftcol, labcex = 0.85){

  nf = dim(ts)[3]
  nsim = dim(ts)[1]
  yrs = as.numeric(dimnames(ts)[[2]])
  ylim = c(0,max(ts)*1.05)
  plot(range(yrs),ylim,col="white",xlab="",ylab=""); grid()
  for(ff in 1:nf)matplot(yrs,t(ts[,,ff]),type="l",col=ftcol[ff],lty=1:nsim,lwd=1,add=T)
  mtext(ylab,2,line=2.4,cex = labcex)
  mtext("Year",1,line=2.4,cex=labcex)
}

get_Landings = function(obj){
  temp = lapply(obj,function(x)x[[1]]@Landings@Value)
  nsim = length(obj)
  nyr = nrow(temp[[1]])
  nft = ncol(temp[[1]])
  temparr = aperm(array(unlist(temp),c(nyr, nft, nsim)),c(3,1,2))
  dimnames(temparr) = list(paste0("sim_",1:nsim),rownames(temp[[1]]), colnames(temp[[1]]))
  temparr
}

convyr = function(ts, Seasons){
  ds = dim(ts)
  temp = apply(array(ts, c(ds[1],Seasons, ds[2]/Seasons, ds[3])),c(1,3,4),sum)
  dimnames(temp) = list(dimnames(ts)[[1]],unique(floor(as.numeric(dimnames(ts)[[2]]))),dimnames(ts)[[3]])
  temp
}


get_CPUE = function(obj){
  temp = lapply(obj,function(x)x[[1]]@CPUE@Value)
  nsim = length(obj)
  nyr = nrow(temp[[1]])
  nft = ncol(temp[[1]])
  temparr = aperm(array(unlist(temp),c(nyr, nft, nsim)),c(3,1,2))
  dimnames(temparr) = list(paste0("sim_",1:nsim),rownames(temp[[1]]), colnames(temp[[1]]))
  temparr

}


process_simdata = function(obj){

    Landings = get_Landings(obj)
  Landings_y = convyr(Landings, Seasons = Seasons)

  CPUE = get_CPUE(obj)
  CPUE_y = convyr(CPUE, Seasons = Seasons)

}

slplot.Data = function(obj, sims = 1:2, ftcol = c("red","blue","green","black","darkgrey","purple")){

  Seasons = obj[[1]][[1]]@Seasons
  par(mfrow=c(3,4),mai=c(0.5,0.5,0.2,0.05))

  # Landings by fleet
  Landings = get_Landings(obj)
  tsdatplot_byft(Landings[sims,,,drop=F],"Landings by Fleet, Season (kg)",ftcol)

  # Annual landings by fleet
  Ly = convyr(Landings, Seasons)
  tsdatplot_byft(ts = Ly[sims,,,drop=F],"Annual Landings by Fleet (kg)",ftcol)

  # Total annual landings
  Ty = array(apply(Ly, 1:2,sum),c(dim(Ly)[1:2],1))
  dimnames(Ty)[1:2] = dimnames(Ly)[1:2]
  tsdatplot_byft(ts = Ty[sims,,,drop=F],"Total Annual Landings (kg)","black")

  # CPUE by fleet
  CPUE = get_CPUE(obj)
  tsdatplot_byft(CPUE[sims,,,drop=F],"CPUE by Fleet, Season",ftcol)

  # Annual CPUE by fleet
  CPUE_y = convyr(CPUE, Seasons)
  tsdatplot_byft(ts = CPUEy[sims,,,drop=F],"Annual CPUE by Fleet",ftcol)



}
