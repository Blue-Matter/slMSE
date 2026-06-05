# RCM Rapid Conditioning Model of Huynh 2026 - https://samtool.openmse.com/reference/RCM.html

#  c_oe = 0.05; i_oe = 0.2; C_eq_fac = 1; C_eq_nyrs = 3; nsubyr = 4; Name = "An RCM OM"; R0init = 1E7; M=0.4
# Len_age = NA; Wt_age = NA; Mat_age = NA; Sel_age = NA
#  Steepness = 0.8; SRrel = 2; oe = 0.2; pe = 0.6; proyears = 15; CurrentYear = 2026

RCM_data = function(sim, simdata, Name = "An RCM OM w data", c_oe = 0.05, i_oe = 0.2,
                    ESS = 50, C_eq_fac = 1, C_eq_nyrs = 5, nsubyr = 4,
                    R0init, M, Len_age = NA, Wt_age = NA, Mat_age = NA, Sel_age = NA,
                    Steepness, SRrel = 2, pe = 0.6,  # SRrel = 1 is B-H, SRrel = 2 is Ricker
                    proyears = 15, CurrentYear = 2026, nSim = 2){


  nsim = 2
  dat = new('RCMdata')
  nfleet = dim(simdata$Landings_q)[3]
  nyear = dim(simdata$Landings_q)[2]

  if(is.na(Wt_age[1])){  Wt_age  = c(0, simdata$wt_age_guess);  cat("Weight at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}
  if(is.na(Len_age[1])){ Len_age = c(0, simdata$len_age_guess); cat("Length at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}
  if(is.na(Mat_age[1])){ Mat_age = c(0, simdata$mat_age_guess); cat("Maturity at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}
  if(is.na(Sel_age[1])){ Sel_age = c(0, simdata$sel_age_guess); cat("Selectivity at age not specified, taking a bad guess from simulations (averaged over stocks and simulations) \n")}

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
  cpars$Data = dat
  OM@cpars = cpars

  OM@L5 = rep(max(Len_age)*0.4,2)
  OM@LFS = rep(max(Len_age)*0.7,2)
  OM@Vmaxlen = rep(0.99,2)

  # Nuisance slots
  OM = fill_OM(OM)
  OM@nsim = nSim

  OM
}

# helpers:
# input = readRDS("C:/GitHub/csrf_hh_data/OMs/Geoduck/Reference_Case/Objects/RCMinput_106.rda")


#   sim = 1; mode = "ASPM"; Name = "An RCM OM w data"; c_oe = 0.05; i_oe = 0.2; ESS = 50; C_eq_fac = 1; C_eq_nyrs = 5; nsubyr = 4
#   R0init = 1E7; M = 0.5; Len_age = NA; Wt_age = NA; Mat_age = NA; Sel_age = NA; Steepness = 0.8; SRrel = 2; pe = 5.0;  max_F = 3.0
#   proyears = 15; CurrentYear = 2026; silent = T; plot = T; resample = F

do_RCM=function(sim, simdata, mode = "ASPM",Name = "An RCM OM w data", c_oe = 0.05, i_oe = 0.2,
                ESS = 50, C_eq_fac = 1, C_eq_nyrs = 5, nsubyr = 4,
                R0init = 1E7, M = 0.5, Len_age = NA, Wt_age = NA, Mat_age = NA, Sel_age = NA,
                Steepness = 0.8, SRrel = 2, pe = 5.0,  max_F = 3.0, # SRrel = 1 is B-H, SRrel = 2 is Ricker
                proyears = 15, CurrentYear = 2026, silent = T, plot = F, resample = F){

  tslab = round(as.numeric(dimnames(simdata$Landings_q)[[2]]),2)
  fnam = dimnames(simdata$Landings_q)[[3]]

  Sdata = RCM_data(sim, simdata, Name, c_oe, i_oe, ESS, C_eq_fac, C_eq_nyrs, nsubyr, R0init, M, Len_age, Wt_age, Mat_age, Sel_age, Steepness, SRrel, pe, proyears, CurrentYear)

  if(mode == 'ASPM'){
    Sdata@cpars$Data@CAL = array()
    Sdata@cpars$V_real = Sdata@cpars$Mat_age
  }

  fit = RCM(Sdata, Sdata@cpars$Data, s_selectivity=c("B",1,2),
               max_F = max_F, mean_fit = F, condition = "catch", cores = 1,
               comp_like="multinom", drop_nonconv=T, drop_highF=T, resample = resample, silent=silent,
               pbc_recdev = rep(0, Sdata@nyears), pbc_earlyrecdev = rep(0, Sdata@maxage))


  if(plot) plot(fit, Year = tslab, f_nam = fnam, s_name = c("Survey",paste0("CR_index ",fnam)))

  RCM_output(fit, simdata, Year = tslab, f_nam = fnam, s_name = c("Survey",paste0("CR_index ",fnam)))

}

RCM_output = function(fit, simdata, Year, f_nam, s_name){

  if(class(fit)!="RCModel")stop("Not an object of class RCModel - an RCM fit")
  fnam = dimnames(simdata$Landings_q)[[3]]
  predts = round(as.numeric(dimnames(simdata$Landings_q)[[2]]),2)
  nts = length(predts)
  histRCM = Simulate(fit@OM)
  BMSY = histRCM@Ref$ReferencePoints$BMSY[1]
  MSY = histRCM@Ref$ReferencePoints$MSY[1]
  UMSY = MSY / BMSY
  BMSY_B0 = histRCM@Ref$ReferencePoints$BMSY[1]/histRCM@Ref$ReferencePoints$B0
  B = fit@report[[1]]$B[1:nts]
  C = apply(fit@data@Chist,1,sum) # observed catches (as in sim test) vs predicted catches: apply(fit@report[[1]]$Cpred,1,sum)
  U = C/B
  list(B = B, C = C, U = U, BMSY = BMSY, MSY = MSY, UMSY = UMSY, BMSY_B0 = BMSY_B0, predts = predts, fit=fit, Year = Year, f_nam = f_nam, s_name= s_name)
}


# mode = "ASPM" ; Name = "An RCM OM"; c_oe = 0.05; i_oe = 0.2; ESS = 50; C_eq_fac = 1; C_eq_nyrs = 3; nsubyr = 4;  R0init = 1E7; M=0.6; Len_age = NA; Wt_age = NA; Mat_age = NA; Sel_age = NA; Steepness = 0.8; SRrel = 2; pe = 5.0;  max_F = 3.0; proyears = 15; CurrentYear = 2026

SimSam_RCM = function(simdata, mode = "ASPM", parallel = T, Name = "An RCM OM w data", c_oe = 0.05, i_oe = 0.2,
                      ESS = 50, C_eq_fac = 1, C_eq_nyrs = 5, nsubyr = 4,
                      R0init = 1E7, M = 0.5, Len_age = NA, Wt_age = NA, Mat_age = NA, Sel_age = NA,
                      Steepness, SRrel = 2, pe = 5.0,  max_F = 3.0, # SRrel = 1 is B-H, SRrel = 2 is Ricker
                      proyears = 15, CurrentYear = 2026){

  if(parallel){
    setup()
    sfLibrary(slMSE)
    sfLibrary(SAMtool)
    Est = sfLapply(1:nSim, do_RCM, simdata = simdata, mode = mode, Name = Name, c_oe = c_oe, i_oe = i_oe,
                   ESS = ESS, C_eq_fac = C_eq_fac, C_eq_nyrs = C_eq_nyrs, nsubyr = nsubyr,
                   R0init = R0init, M = M, Len_age = Len_age, Wt_age = Wt_age, Mat_age = Mat_age, Sel_age = Sel_age,
                   Steepness = Steepness, SRrel = SRrel, pe = pe,  max_F = max_F, # SRrel = 1 is B-H, SRrel = 2 is Ricker
                   proyears = proyears, CurrentYear = CurrentYear)
  }else{
    Est = lapply(1:nSim,do_RCM, simdata = simdata, mode = mode, Name = Name, c_oe = c_oe, i_oe = i_oe,
                 ESS = ESS, C_eq_fac = C_eq_fac, C_eq_nyrs = C_eq_nyrs, nsubyr = nsubyr,
                 R0init = R0init, M = M, Len_age = Len_age, Wt_age = Wt_age, Mat_age = Mat_age, Sel_age = Sel_age,
                 Steepness = Steepness, SRrel = SRrel, pe = pe,  max_F = max_F, # SRrel = 1 is B-H, SRrel = 2 is Ricker
                 proyears = proyears, CurrentYear = CurrentYear)
  }


  SimSam = list()
  estU = t(sapply(Est,function(x)x$U))
  estB = t(sapply(Est,function(x)x$B))

  estBMSY = sapply(Est,function(x)x$BMSY)
  estMSY = sapply(Est,function(x)x$MSY)
  estUMSY = estMSY / estBMSY

  Usim = simdata$U_q
  Bsim = simdata$B_q

  SimSam$U = list(sim = Usim, est = estU)
  SimSam$B = list(sim = Bsim, est = estB)
  SimSam$BMSY = list(sim = simdata$BMSY, est = estBMSY)
  SimSam$MSY = list(sim = simdata$MSY, est = estMSY)
  SimSam$UMSY = list(sim = simdata$UMSY, est = estUMSY)
  class(SimSam) = "slSimSam"
  SimSam
}


# RCM helpers:

fill_OM = function(OM){

  # Length requirements
  OM@L50 = rep(getL50(OM)$L50,2)
  OM@L50_95 = rep(getL50(OM)$L50_95,2)

  # Depletion requirements
  OM@D = rep(0.5,2)

  #Observation
  OM@Cobs = rep(0.025, 2)     # Catch observation error
  OM@Cbiascv = rep(0, 2)      # unbiased
  OM@CAA_nsamp = rep(200, 2)  # no CAA data are actually sampled
  OM@CAA_ESS = rep(100, 2)    # assume effective sample size of 100
  OM@CAL_nsamp = rep(2000, 2) # between 1300-3500 2018-2022
  OM@CAL_ESS = rep(100, 2)    # assume effective sample size is around 100
  OM@Iobs = rep(0.15, 2)      # hypothetical biomass indices observed with a 15% CV
  OM@Btobs = rep(0.15, 2)     # hypothetical absolute biomass estimates observed with a 15% CV
  OM@Dobs = rep(0.2, 2)       # hypothetical depletion observations have CV of 20%
  OM@Eobs = rep(0.1, 2)       # hypothetical effort observations have a CV of 10%

  # No bases in simulated data:
  OM@Btbiascv = OM@LenMbiascv = OM@Mbiascv = OM@Kbiascv =
    OM@t0biascv = OM@Linfbiascv = OM@LFCbiascv = OM@LFSbiascv =
    OM@FMSY_Mbiascv = OM@BMSY_B0biascv = OM@Irefbiascv = OM@Brefbiascv =
    OM@Crefbiascv = OM@Dbiascv = OM@hbiascv = OM@Recbiascv = OM@sigmaRbiascv =
    OM@Ebiascv = rep(0.001, 2)

  OM@beta = c(1, 1)          # assuming relative abundance index is linearly related

  #Implementation
  OM@TACFrac = rep(1, 2)     # TAC implemented without consistent overages or underages
  OM@TACSD = rep(0, 2)       # no variability in TAC implementation (annual variation around recommendation)
  OM@TAEFrac = rep(1, 2)     # hypothetical TAE implemented without consistent overages or underages
  OM@TAESD = rep(0, 2)       # no variability in hypothetical TAE implementation (annual variation around recommendation)
  OM@SizeLimFrac = rep(1, 2) # hypothetical Size limit implemented without consistent overages or underages
  OM@SizeLimSD = rep(0, 2)   # no variability in hypothetical Size limit implementation (annual variation around recommendation)

  OM
}

getL50 = function(OM){

  len = OM@cpars$Len_age[1,,1]
  mat = OM@cpars$Mat_age[1,,1]
  suppressWarnings({
    L50 = approx(mat,len,0.5)$y
    L50_95 = approx(mat,len,0.95)$y - L50
  })
  list(L50 = L50, L50_95 = L50_95)
}

doMSYcalc=function(x, MSYin){
  out = MSEtool:::MSYCalcs(
    logF = MSYin$logF_s[x],
    M_at_Age = MSYin$M_at_Age_s[x,],
    Wt_at_Age = MSYin$Wt_at_Age_s[x,],
    Mat_at_Age = MSYin$Mat_at_Age_s[x,],
    Fec_at_Age = MSYin$Fec_at_Age_s[x,],
    V_at_Age = MSYin$V_at_Age_s[x,],
    Wt_at_Age_C = MSYin$Wt_at_Age_s[x,],
    maxage = MSYin$maxage,
    relRfun = function(){},
    SRRpars=list(),
    R0x = MSYin$R0x_s[x],
    SRrelx = 1,
    hx = MSYin$hx_s[x],
    SSBpR = MSYin$SSBpR_s[x],
    opt=2L)
  out

}


