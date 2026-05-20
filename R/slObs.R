
slObs = function(Landings_CV = 0.05, Effort_CV = 0.1, Survey_CV = 0.15, CPUE_CV = 0.25, CAL_ESS = 400){

  obs = Obs()
  obs@Landings@CV = Landings_CV
  obs@Effort@CV = Effort_CV
  obs@Survey@CV = Survey_CV
  obs@CPUE@CV = CPUE_CV
  obs@LandingsAtSize@ESS = CAL_ESS
  obs@LandingsAtSize@SampleSize = CAL_ESS
  obslist=list()
  obslist[[1]]=list()
  obslist[[1]][[1]] = obs
  obslist

}
