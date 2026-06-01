# RCM Rapid Conditioning Model of Huynh 2026 - https://samtool.openmse.com/reference/RCM.html

#  c_oe = 0.05; i_oe = 0.2; C_eq_fac = 1; C_eq_nyrs = 3; nsubyr = 4; Name = "An RCM OM"; R0init = 1E7; M=0.4
# Len_age = NA; Wt_age = NA; Mat_age = NA; Sel_age = NA
#  Steepness = 0.8; SRrel = "Ricker"; oe = 0.2; pe = 0.6; proyears = 15; CurrentYear = 2026

RCM_data = function(sim, simdata, Name = "An RCM OM w data", c_oe = 0.05, i_oe = 0.2,
                    ESS = 50, C_eq_fac = 1, C_eq_nyrs = 5, nsubyr = 4,
                    R0init, M, Len_age = NA, Wt_age = NA, Mat_age = NA, Sel_age = NA,
                    Steepness, SRrel = 2, pe = 0.6,  # SRrel = 1 is B-H, SRrel = 2 is Ricker
                    proyears = 15, CurrentYear = 2026){


  nsim = 1
  dat = new('RCMdata')
  nfleet = dim(simdata$Landings_q)[3]
  nyear = dim(simdata$Landings_q)[2]

  if(is.na(Wt_age[1])){Wt_age = c(0,simdata$wt_age_guess); cat("Weight at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}
  if(is.na(Len_age[1])){Len_age = c(0,simdata$len_age_guess); cat("Length at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}
  if(is.na(Mat_age[1])){Mat_age = c(0,simdata$mat_age_guess); cat("Maturity at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}
  if(is.na(Sel_age[1])){Sel_age = c(0,simdata$sel_age_guess); cat("Selectivity at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}

  # Landings
  dat@Chist = simdata$Landings_q[sim,,]
  dat@C_sd = array(rep(c_oe,each=nyear),dim(dat@Chist))

  # Equilibrium Landings
  dat@C_eq = apply(dat@Chist[1:(nsubyr*C_eq_nyrs),],2,mean)
  cbys = array(dat@Chist[1:(nsubyr*C_eq_nyrs),],c(nsubyr,C_eq_nyrs,nfleet))
  smu = apply(cbys,2:3,mean)
  scv = apply(smu,2,sd)/apply(smu,2,mean)
  dat@C_eq_sd = scv

  # Surveys
  dat@Index = cbind(simdata$Survey_q[sim,,1],simdata$CPUE_q[sim,,])
  dat@I_sd = array(rep(i_oe,each=nyear),dim(dat@Index))

  # Fishery-dependent indices
  dat@length_bin = as.numeric(dimnames(simdata$CAL_q)[[5]])
  dat@CAL = apply(simdata$CAL_q[sim,,,,],c(3,4,2),sum)
  dat@CAL_ESS = array(rep(ESS,nyear),c(nyear, nfleet))

  allyears = nyear + proyears
  nage = length(Wt_age)
  maxage = nage-1 # age zero formulation

  OM = new('OM')
  OM@Name = Name
  OM@R0 = rep(R0init, 2)
  OM@M = rep(M, 2)
  OM@h = rep(Steepness,2)
  OM@SRrel = SRrel
  OM@proyears = proyears
  OM@nyears = nyear
  OM@maxage = maxage
  OM@nsim = 2
  OM@Perr = rep(pe, 2)

  cpars = list()



  cpars$Wt_age  = array(rep(Wt_age,  each=nsim), c(nsim, nage, allyears))
  cpars$Mat_age = array(rep(Mat_age, each=nsim), c(nsim, nage, allyears))
  cpars$Len_age = array(rep(Len_age, each=nsim), c(nsim, nage, allyears))
  cpars$V_real  = array(rep(Sel_age, each=nsim), c(nsim, nage, allyears))

  OM@L5 = rep(max(Len_age)*0.4,2)
  OM@LFS = rep(max(Len_age)*0.7,2)
  OM@Vmaxlen = rep(0.99,2)

  # Nuisance slots
  OM@D = rep(0.5,2)

  cpars$Data = dat

  OM@cpars = cpars

  OM
}

# helpers:
# input = readRDS("C:/GitHub/csrf_hh_data/OMs/Geoduck/Reference_Case/Objects/RCMinput_106.rda")


#  plot = F; sim = 1; c_oe = 0.05; i_oe = 0.2; C_eq_fac = 1; C_eq_nyrs = 3; nsubyr = 4; Name = "An RCM OM"; R0init = 1E7; M=0.6; Len_age = NA; Wt_age = NA; Mat_age = NA; Sel_age = NA; Steepness = 0.8; SRrel = 2; oe = 0.2; pe = 5.0;  max_F = 3.0; proyears = 15; CurrentYear = 2026

do_RCM=function(sim, simdata, mode = "ASPM",Name = "An RCM OM w data", c_oe = 0.05, i_oe = 0.2,
                ESS = 50, C_eq_fac = 1, C_eq_nyrs = 5, nsubyr = 4,
                R0init, M, Len_age = NA, Wt_age = NA, Mat_age = NA, Sel_age = NA,
                Steepness, SRrel = 2, pe = 5.0,  max_F = 3.0, # SRrel = 1 is B-H, SRrel = 2 is Ricker
                proyears = 15, CurrentYear = 2026, plot=F){

  Sdata = RCM_data(sim, simdata, Name, c_oe, i_oe, ESS, C_eq_fac, C_eq_nyrs, nsubyr, R0init, M, Len_age, Wt_age, Mat_age, Sel_age, Steepness, SRrel, pe, proyears, CurrentYear)

 #  Sinput = SPiCT_config(Sdata, r.pr, bk.pr, shape.pr, oe, pe, fdevs, ce,  q.pr, timing, dteuler)
  #ck = check.inp(Sinput)
  #fit = fit.spict(Sinput)

  if(mode == 'ASPM'){
    Sdata@cpars$Data@CAL = array()
    Sdata@cpars$V_real = Sdata@cpars$Mat_age
  }

  RCMfit = RCM(Sdata, Sdata@cpars$Data, s_selectivity=c("B",1,2),
               max_F = max_F, mean_fit = F, condition = "catch", cores = 1,
               comp_like="multinom", drop_nonconv=T, drop_highF=T, resample = F, silent=F,
               pbc_recdev = rep(0, Sdata@nyears), pbc_earlyrecdev = rep(0, Sdata@maxage))


  # Year = Year, s_name = fitobj@OM@Misc$Slabs, f_name = fitobj@OM@Misc$Flabs
  if(plot) plot(RCMfit)

  # You got here!!
  RCM_output(fit)

 # s_selectivity = c(rep("dome_age",nsub-1),"logistic_age",rep(1,isCPUE)


}

RCM_output = function(fit){
  predts = fit$inp$time
  BMSY = getsval(fit, "Bmsy",ts=F)
  MSY = getsval(fit, "MSY",ts=F)
  BMSY_B0 = exp(getsval(fit, "logBpK",ts=F))
  Brel = exp(getsval(fit, nam = "logBBmsy"))
  B = BMSY * Brel
  C = exp(getsval(fit, nam = "logCpred"))
  U = C/B
  list(B = B, C = C, U = U, BMSY = BMSY, MSY = MSY, BMSY_B0 = BMSY_B0, predts = predts, fit=fit)
}

SimSam_RCM = function(simdata, timestep = "year", parallel =T,
                        r.pr = c(0.5,0.2,1),
                        bk.pr = c(0.5,0.3,1),
                        shape.pr = c(2, 0.001, 1),
                        oe = c(0.2, 0.5, 1),
                        pe = c(0.2,0.5,1),
                        fdevs = c(4, 0.5, 1),
                        ce = c(0.05, 0.001, 1),
                        q.pr = NULL,
                        timing = 0.625,
                        dteuler = 0.25){

  if(parallel){
    setup()
    sfLibrary(slMSE)
    sfLibrary(spict)
    Est = sfLapply(1:nSim,do_spict, simdata = simdata, timestep = timestep, r.pr = r.pr, bk.pr = bk.pr,
                   shape.pr = shape.pr, oe = oe, pe = pe, fdevs = fdevs, ce = ce,
                   q.pr = q.pr, timing = timing, dteuler = dteuler)
  }else{
    Est = lapply(1:nSim,do_spict, simdata = simdata, timestep = timestep, r.pr = r.pr, bk.pr = bk.pr,
                 shape.pr = shape.pr, oe = oe, pe = pe, fdevs = fdevs, ce = ce,
                 q.pr = q.pr, timing = timing, dteuler = dteuler)
  }


  SimSam = list()
  estU = t(sapply(Est,function(x)x$U))
  estB = t(sapply(Est,function(x)x$B))

  estBMSY = sapply(Est,function(x)x$BMSY)
  estMSY = sapply(Est,function(x)x$MSY)
  estUMSY = estMSY / estBMSY

  if(timestep == "year"){
    Usim = simdata$U_y
    Bsim = simdata$B_y
  }else if(timestep =="quarter"){
    Usim = simdata$U_q
    Bsim = simdata$B_q
  }
  SimSam$U = list(sim = Usim, est = estU)
  SimSam$B = list(sim = Bsim, est = estB)
  SimSam$BMSY = list(sim = simdata$BMSY, est = estBMSY)
  SimSam$MSY = list(sim = simdata$MSY, est = estMSY)
  SimSam$UMSY = list(sim = simdata$UMSY, est = estUMSY)
  class(SimSam) = "slSimSam"
  SimSam
}




