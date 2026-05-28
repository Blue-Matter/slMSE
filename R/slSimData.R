

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


get_Survey = function(obj){
  temp = lapply(obj,function(x)x[[1]]@Survey@Value)
  nsim = length(obj)
  nyr = nrow(temp[[1]])
  nft = ncol(temp[[1]])
  temparr = aperm(array(unlist(temp),c(nyr, nft, nsim)),c(3,1,2))
  dimnames(temparr) = list(paste0("sim_",1:nsim),rownames(temp[[1]]), colnames(temp[[1]]))
  temparr

}

get_CAL = function(obj){
  #obj@LandingsAtSize

}


sampleCAL = function(NAA, sel){


}

invent_CAL = function(obj, ESS = 200){ # hack to create mulitnomial catch at length by fleet and stock

  Seasons = Seasons(obj)
  ns = nStock(obj)
  nf = nFleet(obj)
  CALmids = obj@OM@Stock[[1]]@Length@Classes
  nl = length(CALmids)
  years = Years(obj, Period="Historical")
  ny = length(years)
  nsim = nSim(obj)
  CAL = array(0, c(nsim,ns,nf,ny,nl))
  NAAs = obj@Number

  for(sim in 1:nsim){
    for(ss in 1:ns){
      ALK = obj@OM@Stock[[ss]]@Length@ALK[1,,,1]
      for(ff in 1:nf){
        for(yy in 1:ny){
          Nage = apply(NAAs[[ss]][sim,,yy,,drop=F],2,sum)
          sel = obj@OM@Fleet[[ss]][[ff]]@Selectivity@MeanAtAge[1,,1,1]
          Nlen = (sel*Nage)%*%ALK
          CAL[sim,ss,ff,yy,] = rmultinom(1,ESS,Nlen)
        }
      }
    }
  }

  dimnames(CAL) = list(paste0("sim_",1:nsim), names(obj@OM@Stock), names(obj@OM@Fleet[[1]]),years,CALmids)
  CAL_y = apply(array(CAL,c(nsim,ns,nf,Seasons(obj),nYear(obj),nl)),c(1,2,3,5,6),sum)
  dimnames(CAL_y) = list(paste0("sim_",1:nsim), names(obj@OM@Stock), names(obj@OM@Fleet[[1]]),CalcYears(nYear(obj),0,obj@OM@CurrentYear,1),CALmids)
  list(CAL = CAL, CAL_y = CAL_y)
}

slSimData = function(obj){

  dat = obj@Data

  Landings = get_Landings(dat)
  Landings_y = convyr(Landings, Seasons = Seasons)

  CPUE = get_CPUE(dat)
  CPUE_y = convyr(CPUE, Seasons = Seasons)

  Survey = get_Survey(dat)
  Survey_y = convyr(Survey, Seasons = Seasons)

  # CAL = get_CAL(dat) # coming when obs provides LandingsAtSize
  CALs = invent_CAL(obj)
  CAL = CALs$CAL
  CAL_y = CALs$CAL_y

  # Dimensions
  Seasons = Seasons(obj)
  nStock = nStock(obj)
  nFleet = nFleet(obj)
  CALmids = obj@OM@Stock[[1]]@Length@Classes
  nLen = length(CALmids)
  years = Years(obj, Period="Historical")
  nYear = length(years)
  nSim = nSim(obj)

  # Sim quantities
  B = t(sapply(1:nSim(obj),function(x, obj){get_sim_B(x, obj)}, obj=obj))
  U = apply(Landings_y,1:2,sum)/B

  slsd = list(Landings = Landings, Landings_y = Landings_y,
       CPUE = CPUE, CPUE_y = CPUE_y,
       Survey = Survey, Survey_y = Survey_y,
       CAL = CAL, CAL_y = CAL_y, nSeason = Seasons, nStock = nStock, nFleet = nFleet,
       nLen = nLen, nYear = nYear, nSim = nSim, CALmids = CALmids, years = years,
       B = B, U = U)

  class(slsd) = "slSimData"
  slsd

}



get_sim_B = function(sim, hist){

  blist = list()
  for(ss in 1:nStock(hist)){
    N = hist@Number[[ss]][sim,,,]
    W = hist@OM@Stock[[ss]]@Weight@MeanAtAge[1,,] # doesn't vary with sim
    Btemp = apply(array(N*array(W,dim(N)),dim(N)),2,sum) # sum over ages and areas
    blist[[ss]] = apply(array(Btemp,c(Seasons(hist),nYear(hist))),2,mean)
  }
  allB = array(unlist(blist),c(nYear(hist),nStock(hist)))
  sumB = apply(allB,1,sum) # sum over all stocks
  names(sumB) = CalcYears(nYear(hist),0,hist@OM@CurrentYear)
  sumB
}

