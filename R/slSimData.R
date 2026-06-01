

get_Landings = function(obj){
  temp = lapply(obj,function(x)x[[1]]@Landings@Value)
  nsim = length(obj)
  nyr = nrow(temp[[1]])
  nft = ncol(temp[[1]])
  temparr = aperm(array(unlist(temp),c(nyr, nft, nsim)),c(3,1,2))
  dimnames(temparr) = list(paste0("sim_",1:nsim),rownames(temp[[1]]), colnames(temp[[1]]))
  temparr
}


convq = function(ts, Seasons){
  if(Seasons ==4){
    temp = ts
  }else{
    tslab = floor(as.numeric(colnames(ts)))
    yrs = unique(tslab)
    n2agg = Seasons/4
    ds = dim(ts)
    temp = apply(array(ts, c(ds[1],n2agg, ds[2]/n2agg, ds[3])), c(1,3,4), sum)
    dimnames(temp) = list(dimnames(ts)[[1]],rep(yrs,each=4)+rep(c(0,0.25,0.5,0.75),length(yrs)),dimnames(ts)[[3]])
  }

  temp
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
  if(Seasons == 4){
    CAL_q = CAL
  }else{
    n2agg = Seasons/4
    tslab = Years(obj,"Hist")
    CAL_q = apply(array(CAL, c(nsim,ns,nf, n2agg, length(tslab)/n2agg,nl),c(1,2,3,5,6)),sum)
    yrs = CalcYears(nYear,0,hist@OM@CurrentYear)
    tslab = rep(yrs,each=2)+rep((0:3)/4,length(yrs))
    dimnames(CAL_q) = list(paste0("sim_",1:nsim), names(obj@OM@Stock), names(obj@OM@Fleet[[1]]),tslab,CALmids)
  }
  list(CAL = CAL, CAL_q = CAL_q, CAL_y = CAL_y)
}

slSimData = function(obj){

  dat = obj@Data

  Landings = get_Landings(dat)
  Landings_q = convq(Landings, Seasons = Seasons)
  Landings_y = convyr(Landings, Seasons = Seasons)

  CPUE = get_CPUE(dat)
  CPUE_q = convq(CPUE, Seasons = Seasons)
  CPUE_y = convyr(CPUE, Seasons = Seasons)

  Survey = get_Survey(dat)
  Survey_q = convq(Survey, Seasons = Seasons)
  Survey_y = convyr(Survey, Seasons = Seasons)

  # CAL = get_CAL(dat) # coming when obs provides LandingsAtSize
  CALs = invent_CAL(obj)
  CAL = CALs$CAL
  CAL_q = CALs$CAL_q
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
  B_q = t(sapply(1:nSim(obj),function(x, obj){get_sim_Bq(x, obj)}, obj=obj))
  U_q = apply(Landings_q,1:2,sum)/B_q
  B_y = t(sapply(1:nSim(obj),function(x, obj){get_sim_By(x, obj)}, obj=obj))
  U_y = apply(Landings_y,1:2,sum)/B_y


  # Spec quantities
  getmuage=function(x,slt)apply(slot(x,slt)@MeanAtAge, 2, mean)
  len_age_guess = apply(sapply(obj@OM@Stock, getmuage, slt="Length"), 1, mean)
  wt_age_guess = apply(sapply(obj@OM@Stock, getmuage, slt="Weight"), 1, mean)
  mat_age_guess = apply(sapply(obj@OM@Stock, getmuage, slt="Maturity"), 1, mean)
  mat_age_guess = mat_age_guess/max(mat_age_guess)
  sel_age_guess = apply(sapply(obj@OM@Fleet[[1]], getmuage, slt="Selectivity"), 1, mean)


  slsd = list(Landings = Landings, Landings_q = Landings_q, Landings_y = Landings_y,
       CPUE = CPUE, CPUE_q = CPUE_q, CPUE_y = CPUE_y,
       Survey = Survey, Survey_q = Survey_q, Survey_y = Survey_y,
       CAL = CAL, CAL_q = CAL_q, CAL_y = CAL_y,
       nSeason = Seasons, nStock = nStock, nFleet = nFleet,
       nLen = nLen, nYear = nYear, nSim = nSim, CALmids = CALmids, years = years,
       B_y = B_y, B_q = B_q, U_y = U_y, U_q = U_q,
       BMSY = rep(1,nSim), MSY = rep(1,nSim), UMSY = rep(1,nSim), BMSY_B0 = rep(1,nSim),
       len_age_guess = len_age_guess, wt_age_guess = wt_age_guess,
       mat_age_guess = mat_age_guess, sel_age_guess = sel_age_guess)

  class(slsd) = "slSimData"
  slsd

}



get_sim_By = function(sim, hist){

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




get_sim_Bq = function(sim, hist){

  blist = list()
  Seasons = Seasons(hist)
  if(Seasons == 4){
    for(ss in 1:nStock(hist)){
      N = hist@Number[[ss]][sim,,,]
      W = hist@OM@Stock[[ss]]@Weight@MeanAtAge[1,,] # doesn't vary with sim
      blist[[ss]] = apply(array(N*array(W,dim(N)),dim(N)),2,sum) # sum over ages and areas
    }
  }else{
    nt = length(Years(hist,"Hist")) #length(CalcYears(nYear(hist),0,hist@OM@CurrentYear))
    n2agg = Seasons/4
    for(ss in 1:nStock(hist)){
      N = hist@Number[[ss]][sim,,,]
      W = hist@OM@Stock[[ss]]@Weight@MeanAtAge[1,,] # doesn't vary with sim
      Btemp = apply(array(N*array(W,dim(N)),dim(N)),2,sum) # sum over ages and areas
      blist[[ss]] = apply(array(Btemp,c(n2agg,nt/n2agg)),2,mean)
    }
  }
  allB = array(unlist(blist),c(nYear(hist)*4,nStock(hist)))
  sumB = apply(allB,1,sum) # sum over all stocks
  yrs = CalcYears(nYear(hist),0,hist@OM@CurrentYear)
  names(sumB) = rep(yrs,each=4) + rep((0:3)/4,length(yrs))
  sumB
}
