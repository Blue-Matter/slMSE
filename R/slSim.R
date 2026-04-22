# Short-lived Sim
# https://ices-library.figshare.com/WGCEPH
# https://ices-library.figshare.com/articles/report/WGCEPH_2025_Data_call_for_2024_landings_discards_biological_effort_and_survey_data/28495505
# https://sprfmo.int/assets/Meetings/02-SC/13th-SC-2025/Working-Papers/SC13-WP17-SQUIDSIM-Papers-presented-by-Chile.pdf

#' Make a Short-Lived Species Stock Object
#'
#' A function that builds the stock operating model component for a short-lived species from
#' specifications for growth, stock-recruitment, natural mortality rate and spatial distribution. This
#' is a wrapper function and produces a fully customizable MSEtool object of class 'stock'. You can edit
#' that object to include very complex stock dynamics if desired. Note that you can also create lists of these
#' stock objects if you want multiple species or growth-type-groups of the same species.
#'
#' @param Name Character string. The name of the stock object. E.g., "North Pacific Neon Flying Squid"
#' @param Species Character string. The genus and species. E.g., "Ommastrephes bartramii"
#' @param CommonName Character string. The common name. E.g., "Neon Flying Squid"
#' @param nYear Positive Integer. The number of historical years. E.g., 10, (2015 - 2024)
#' @param pYear  Positive Integer. The number of projection years. E.g., 10 (2025 - 2034)
#' @param Seasons Positive Integer. The number of sub-year time steps. E.g., 12 for a monthly model
#' @param CurrentYear Integer. The Calendar year of the last historical year. E.g. 2024.
#' @param nSim Positive Integer. The number of independent simulations. E.g., 4 for testing, 24 for preliminary results, 48 for representative results, 192+ for informing management.
#' @param rec_age Positive Integer. The age the fish is recruited (in seasonal time steps). E.g., given a monthly model (Seasons = 12) a value of 4 means that the stock recruitment function calculates recruiment into the time step 4 months later.
#' @param nages Positive Integer. The number of seasonal time steps that population dynamics will be calculated for. E.g., 24 would be two years in a monthly model. Note that this only necessary if creatures live that long or dynamics are not suitably approximated with a plus group.
#' @param PlusGroup Boolean. Should population numbers be aggregated in the nages age class or just be assumed to die off at nages+1?
#' @param spawndist Numeric vector Seasons long. The seasonal distribution of spawning. E.g., c(0, 0, 0.1, 0.5, 0.3, 0.2, 0, 0, 0, 0, 0, 0), is a monthly spawning pattern where 50 per cent of spawning occurs in April.
#' @param Linf Positive real number or vector of positive real numbers nSim long. The asymptotic length.
#' @param K Positive real number or vector of positive real numbers nSim long. The somatic growth rate (von B.) per season.
#' @param t0 Negative real number or vector of negative real numbers nSim long. The theoretical length at age (seasonal) zero.
#' @param Len_CV Positive real number. The coefficient of variation of length at age. Typically between 0.05 and 0.25.
#' @param a Positive real number. Weight (W) at length (L) parameter a. W = aL^b. Converts units of length to weight.
#' @param b Positive real number. Weight (W) at length (L) parameter b. W = aL^b. Typically approximately cubic (3-ish)
#' @param M Positive real number or vector of positive real numbers nSim long. The instantaneous rate of natural mortality per season.
#' @param amat50 Positive real number or vector of positive real numbers nSim long. The age (in seasons) that 50% of individuals are mature, in a logistic maturity model. E.g., 6 in a monthly model would be 50% mature after 6 months.
#' @param amatSLP Positive real number or vector of positive real numbers nSim long. The slope of the logistic maturity model. E.g., 2 in a monthly model.
#' @param h Positive real number or vector of positive real numbers nSim long. Steepness of the stock-recruitment relationship (Bev-Holt)
#' @param sigmaR Positive real number or vector of positive real numbers nSim long. The lognormal standard deviation in recruitment deviations.
#' @param trunc_sigmaR Positive real number. The number of standard deviations to truncate recruitment deviations. E.g., 2 resamples recruitment deviations if over two standard deviations from the mean (upper tail q of 1.96).
#' @param R_AC Positive real number or vector of positive real numbers nSim long. The lag-1 autocorrelation in recruitment deviations.
#' @param R0 Positive real number or vector of positive real numbers nSim long. Unfished recruitment (number of individuals).
#' @param nareas Positive real number. The number of spatial areas. Should match the spec of Frac_area which has nareas-1 rows.
#' @param Frac_area Matrix of fractions nareas-1 x nages. The fraction of stock numbers found in each area by age. Frac_area[1,3] = 0.9 means that 90 percent of individuals of age 3 are found in area 1.
#' @param prob_stay Fraction. The tendency for individuals to remain in the same area among seasonal time steps. Using gravity equations, the function optimizes for movement that obtains Frac_area while staying as close to prob_stay as possible.
#' @return An object of MSEtool class stock
#' @author T. Carruthers
#' @examples
#' A_short_lived_stock <- slStock("Default stock values for demo")
#' class(A_short_lived_stock)
#' @seealso \link{slFleet} for making a short-lived fleet object and \link{slOM} for specifying the entire operating model from fleet and stock objects.
#' @export
slStock = function(Name = "A short-lived creature", Species = "Shortus liveus", CommonName = "Short-lived creature",
                         nYear = 10, pYear = 10, Seasons = 12, CurrentYear = 2026, nSim = 4,
                         rec_age = 1, nages = 24, PlusGroup = F,
                         spawndist = c(0, 0, 0.1, 0.5, 0.3, 0.2, 0, 0, 0, 0, 0, 0),
                         Linf = 1, K = 0.2, t0 = 0, Len_CV = 0.2, a = 1E-5, b = 3,
                         M = 0.2, amat50 = 6, amatSLP = 2,
                         h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                         nareas = 2, Frac_area = matrix(c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.5,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.01,0.01),nrow=1),
                         prob_stay = 0.9){

  # Default args for testing:
  # Name = "A short-lived creature"; Species = "Shortus liveus"; CommonName = "Short-lived creature"
  # nYear = 10; pYear = 10; Seasons = 12; CurrentYear = 2026; nSim = 4; rec_age = 1; nages = 24; PlusGroup = F
  # Linf = 1; K = 0.2; t0 = 0; Len_CV = 0.2; a = 1; b = 3; M = 0.2; amat50 = 6; amatSLP = 2;
  # h = 0.9; R0 = 1E6, sigmaR = 1.0; trunc_sigmaR = 2.0; R_AC = 0.5; nareas = 2;  prob_stay = 0.9
  # prob_stay = 0.9; Frac_area = matrix(c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.5,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.01,0.01),nrow=1)



  # pre calcs
  na = nages - rec_age + 1

  stock = new('stock')

  # Names and dimensions
  stock@Name = Name
  stock@Species = Species
  stock@CommonName = CommonName
  stock@nYear = nYear
  stock@pYear = pYear
  stock@CurrentYear = CurrentYear
  stock@Years = CurrentYear+ (-(nYear-1):pYear)
  stock@Seasons = Seasons
  stock@nSim = nSim

  # Ages
  stock@Ages = Ages(MaxAge = nages, MinAge = rec_age, Units = "month", PlusGroup = F)

  # Length
  Length = new('length')
  #Length@Model = "vonBert"
  #Length@Pars = list(Linf = rep(Linf,2), K = rep(K, 2), t0 = rep(t0, 2))
  Length@MeanAtAge = Linf * 1-exp(-K * ((rec_age:nages)-t0))
  Length@Units = "mm"
  Length@CVatAge = rep(Len_CV, na)
  Length@Dist = "normal"
  Length@TruncSD = 2
  stock@Length = Length

  # Weight
  Weight = new('weight')
  #Weight@Model = ExampleOM@Stock@Weight@Model #   NULL #MSEtool::WeightatMeanLength()
  #Weight@Pars = list(alpha = a, beta = b)
  Weight@MeanAtAge = a * Length@MeanAtAge ^ b
  Weight@Units = "g"
  Weight@Dist = "lognormal"
  Weight@TruncSD = 2
  stock@Weight = Weight

  # Natural Mortality
  stock@NaturalMortality = NaturalMortality(MeanAtAge = rep(M,na))

  # Maturity
  avec = ((rec_age : nages)-amat50) * amatSLP
  atrans = exp(avec)/(1+exp(avec))
  amat = array(rep(atrans,each=nSim),c(nSim,na,nYear*Seasons+pYear*Seasons))

  stock@Maturity = Maturity(MeanAtAge = amat)
  #stock@Maturity = Maturity(Pars = list(L50 = amat50, L50_95 = amat50/10))


  # Stock-Recruitment
  stock@SRR = SRR(Model = "BevertonHolt", Pars = list(h = h), R0 = R0, SD = sigmaR, AC = R_AC, TruncSD = trunc_sigmaR)

  # Spatial
  # nSim, nArea, nArea, nAge, and nTS,
  Movement = array(0,c(1,nareas,nareas,na,1))
  dist1 = c(Frac_area[,1],1-sum(Frac_area[,1]))
  Movement[1,,,1,1] = matrix(dist1,ncol=nareas,nrow=nareas,byrow=T) # Initial movement is fully mixed
  for(aa in 2:na){
    fracsin = c(Frac_area[,aa-1],1-sum(Frac_area[,aa-1]))
    fracsout = c(Frac_area[,aa],1-sum(Frac_area[,aa]))
    movout = get_mov_D(fracsin, fracsout, prob = rep(prob_stay,nareas))
    Movement[1,,,aa,1] = movout$mov
  }
  stock@Spatial = Spatial(Movement=Movement)

  # Optional slots:
  # Fecundity - assumes proportional to SBiomass if left empty
  # Spatial - leave empty if no spatial structure
  # Depletion - leave empty if you don't need to optimize q for specific depeltion
  # Fleet - Retention: don't need it if all retained

  stock
}


slStock_check = function(stock){
  # stock spec checks

}

#' Make a Short-Lived Species Fleet Object
#'
#' A function that builds the fleet component of an operating model for a short-lived species. Requires from
#' specifications for growth, stock-recruitment, natural mortality rate and spatial distribution. This
#' is a wrapper function and produces a fully customizable MSEtool object of class 'stock'. You can edit
#' that object to include very complex stock dynamics if desired.
#'
#' @param Name Character string. The name of the stock object. E.g., "North Pacific Neon Flying Squid"
#' @param Species Character string. The genus and species. E.g., "Ommastrephes bartramii"
#' @param CommonName Character string. The common name. E.g., "Neon Flying Squid"
#' @param nYear Positive Integer. The number of historical years. E.g., 10, (2015 - 2024)
#' @param pYear  Positive Integer. The number of projection years. E.g., 10 (2025 - 2034)
#' @param Seasons Positive Integer. The number of sub-year time steps. E.g., 12 for a monthly model
#' @param CurrentYear Integer. The Calendar year of the last historical year. E.g. 2024.
#' @param nSim Positive Integer. The number of independent simulations. E.g., 4 for testing, 24 for preliminary results, 48 for representative results, 192+ for informing management.
#' @param rec_age Positive Integer. The age the fish is recruited (in seasonal time steps). E.g., given a monthly model (Seasons = 12) a value of 4 means that the stock recruitment function calculates recruiment into the time step 4 months later.
#' @param nages Positive Integer. The number of seasonal time steps that population dynamics will be calculated for. E.g., 24 would be two years in a monthly model. Note that this only necessary if creatures live that long or dynamics are not suitably approximated with a plus group.
#' @param Effort A matrix of positive real numbers nSim x time steps (nYear x Seasons). Given a default q value of 1 this is the apical fishing mortality rate. Defaul is 'NA' and in this case effort pattern is simulated with a seasonal and temporal trend.
#' @param sel50 Positive real number. The age (in seasons) that 50% of individuals are selected, in a logistic selectivity model. E.g., 6 in a monthly model would be 50% selected after 6 months.
#' @param selSLP Positive real number. The slope of the logistic selectivity model. E.g., 2 in a monthly model.
#' @return An object of MSEtool class fleet
#' @author T. Carruthers
#' @examples
#' A_short_lived_fleet<- slFleet("Default fleet values for demo")
#' class(A_short_lived_fleet)
#' @seealso \link{slFleet} for making a short-lived stock object and \link{slOM} for specifying the entire operating model from fleet and stock objects.
#' @export
slFleet = function(Name = "A fleet", nYear = 10, pYear = 10, Seasons = 12, CurrentYear = 2026,
                         nSim = 4, rec_age = 1, nages = 24, Effort = NA, sel50 = 6, selSLP = 2){

  # Default args for testing:
  # Name = "A fleet"; nYear = 10; pYear = 10; Seasons = 12; nArea = CurrentYear = 2026; nSim = 4; rec_age = 4; nages = 24; Effort = NA;
  # sel50 = 6; selSLP = 2

  fleet = Fleet()
  fleet@Name = Name
  fleet@nYear = nYear
  fleet@pYear = pYear
  fleet@CurrentYear = CurrentYear
  fleet@Years = CurrentYear+ (-(nYear-1):pYear)
  fleet@Seasons = Seasons
  fleet@nSim = nSim

  if(class(Effort)!="array"){
    Effort = Effort_sim(nSim, nYear, Seasons, nArea, ymin = 0.25, yfac = 0.5, ECV = 0.15, maxF = 0.1, plot=T)
  }
  fleet@Effort = Effort(Effort=Effort)

  fleet@Catchability = Catchability(Efficiency = 1)

  svec = ((rec_age : nages)-sel50)*selSLP
  fleet@Selectivity = Selectivity(MeanAtAge = exp(svec)/(1+exp(svec)))

  # Optional slots:
  # Distribution - only needed if you don't want openMSE to solve the distribution according to vulnerable biomass

  fleet

}

SL_fleet_check = function(fleet){
  # fleet spec checks

}

#' Construct Complete Short-Lived Species Operating Model
#'
#' A function that combines stock and fleet object in a complete operating model for a short-lived species. Can
#' include lists of fleets and stocks for more complex operating models
#'
#' @param Name Character string. The name of the operating model. E.g., "NP NFS Ref OM #3: high M, high h, low rec"
#' @param Agency Character string. The relevant fishery agency "North Pacific Fishery Council (NPFC)"
#' @param Author Character string. Name of the author(s) of the operating model E.g., "A. Person"
#' @param Email Character string. Email address of the corresponding author
#' @param Region Character string. The location of the stock / management area.
#' @param Latitude Real number. Degrees N. Negative is southern hemisphere
#' @param Longitude Real number. Degrees E. Negative is western hemisphere
#' @param Sponsor Character string. The name of the person(s) or organization(s) that paid for the research.
#' @param nSim Positive Integer. The number of independent simulations. E.g., 4 for testing, 24 for preliminary results, 48 for representative results, 192+ for informing management.
#' @param nYear Positive Integer. The number of historical years. E.g., 10, (2015 - 2024)
#' @param pYear  Positive Integer. The number of projection years. E.g., 10 (2025 - 2034)
#' @param Seasons Positive Integer. The number of sub-year time steps. E.g., 12 for a monthly model
#' @param CurrentYear Integer. The Calendar year of the last historical year. E.g. 2024.
#' @param Interval Positive integer. The management interval (in seasons) - how frequently is new advice provided. In a monthly model, a value of 6 would mean every 6 months after June and December.
#' @param Seed Real number. The seed for sampling of random variables.
#' @param stock A stock object or list of stock objects
#' @param fleet A fleet object or list of fleet objects
#' @return An object of MSEtool class om
#' @author T. Carruthers
#' @examples
#' A_short_lived_om <- slOM("Default fleet values for demo")
#' class(A_short_lived_om)
#' @seealso \link{slStock} for making a short-lived stock object and \link{slFleet} for making a short-lived fleet object.
#' @export
slOM = function(Name = "Short-lived simulation", Agency = "A fishery agency", Author = "A fishery analyst",
                 Email = "a.person@email.com", Region = "A fishery management area", Latitude = NA, Longitude = NA,
                 Sponsor = "A generous funder", nSim = 4, nYear = 10, pYear = 10, Seasons = 12, CurrentYear = 2026,
                 Interval = 6, Seed = 1, stock = NA, fleet = NA){

  #  Name = "Short-lived simulation"; Agency = "A fishery agency"; Author = "A fishery analyst"; Email = "a.person@email.com"; Region = "A fishery management area"; Latitude = NA; Longitude = NA
  #  Sponsor = "A generous funder"; nSim = 4; nYear = 10; pYear = 10; Seasons = 12; CurrentYear = 2026; Interval = 6; Seed = 1

  om = new('om')     # Operating model object

  # Dimensions and labels
  om@Name = Name
  om@Agency = Agency
  om@Author = Author
  om@Email = "An email address"
  om@Region = "A region"
  om@Latitude = 0
  om@Longitude = 0
  om@Sponsor = "A sponsor"
  om@nSim = nSim
  om@nYear = nYear
  om@pYear = pYear
  om@CurrentYear = CurrentYear
  om@Seasons = Seasons
  om@Interval = Interval
  om@Seed = Seed

  # Observation model --------------------------------------------------------------

  #obs = Obs()
  #obs@Landings@CV = 0.025
  #obs@Survey@CV = 0.1
  #obs@CAL@ESS = 200
  #om@Obs = obs

  # Implementation model -----------------------------------------------------------

  #imp = Imp()
  #imp@TAC@Mean=1
  #imp@TAC@SD = 0.001
  #om@Imp = imp

  # Populate stock and fleet slots -------------------------------------------------

  if(!(class(stock)%in%c("list","stock")))stock = slStock()
  if(!(class(fleet)%in%c("list","fleet")))fleet = slFleet()

  om@Stock = stock
  om@Fleet = fleet
  om = Populate(om)


  list("Astock" = stock)
  class(om@Stock) = "StockList"

  sfl = list(); class(sfl) = "StockFleetList"
  sfl[['Astock']] = list()
  class(sfl[['Astock']]) = "FleetList"
  sfl$Astock$Fleet1 = fleet
  om@Fleet = sfl



  # hist = Simulate(om)

  om

}

slOM_check=function(){
  # om spec checks

}

# Invent Effort

Effort_sim = function(nSim, nYear, Seasons, nArea, ymin = 0.25, yfac = 0.5, ECV = 0.15, maxF = 0.1, plot=F){
  nt =  nYear*Seasons
  Effort = array(0,c(nSim,nt))
  ye = 1 + ymin + (sin(((9+1:(nYear))/3))*1)
  se = dnorm(1:Seasons,Seasons/2,Seasons/5)
  seind = t(array(1:Seasons,c(nt,nSim)))
  yind = t(array(rep(1:nYear,each=Seasons),c(nt,nSim)))
  Effort[] = ye[yind] * se[seind] * trlnorm(nt*nSim,1,ECV)
  Effort = Effort / apply(Effort,1,max) * maxF
  if(plot)matplot(t(Effort),type="l",lty=1)
  Effort
}


dograv2 = function(log_visc,log_grav,fracsin){
  log_grav_1 = c(0,log_grav)
  nareas = length(fracsin)
  lmov = matrix(log_grav_1,byrow=T,nrow=nareas,ncol=nareas)
  diag(lmov) = diag(lmov) + log_visc
  emov = exp(lmov)
  mov = emov/apply(emov,1,sum)
  grav_out=list()
  grav_out$predfracsout = fracsin %*% mov
  grav_out$mov = mov
  grav_out$psum = mean(diag(mov))
  grav_out
}

opt_mov_D <- function(x, fracsin, fracsout, prob, probCV = 0.35, distCV = 0.01) {

  grav_out <- dograv2(log_visc = x[1:length(prob)], log_grav = x[(length(prob)+1):length(x)], fracsin=fracsin)

  nll_prior = dnorm(x,0,5,TRUE)
  nll_dist <- dnorm(log(grav_out$predfracsout), log(fracsout), distCV, TRUE)
  if(length(prob) == 1) {
    nll_stay <- dnorm(log(grav_out$psum), log(prob), probCV, TRUE)
  } else {
    nll_stay <- dnorm(log(diag(grav_out$mov)), log(prob), probCV, TRUE)
  }
  nll <- sum(nll_dist, nll_stay, nll_prior)
  return(-nll)
}

# fracsin = c(0.1, 0.2, 0.3, 0.4); fracsout = c(0.4,0.3,0.2,0.1); prob = c(0.5, 0.8, 0.9, 0.95)
get_mov_D <- function(fracsin = c(0.1, 0.2, 0.3, 0.4), fracsout = c(0.4,0.3,0.2,0.1), prob = c(0.5, 0.8, 0.9, 0.95)) {

  nareas <- length(fracsin)
  nprob <- length(prob)

  opt <- stats::nlminb(rep(0, nprob + nareas - 1), opt_mov_D,
                       fracsin = fracsin, fracsout=fracsout, prob = prob,
                       control = list(iter.max = 5e3, eval.max = 5e3))

  hess = numDeriv::hessian(opt_mov_D,opt$par, prob=prob, fracsin=fracsin,fracsout=fracsout)
  vcv = solve(hess)

  mov <- dograv2(log_visc = opt$par[1:length(prob)],
                 log_grav = opt$par[(length(prob)+1):length(opt$par)],
                 fracsin = fracsin)$mov

  out=list()
  out$opt = opt
  out$vcv = vcv
  out$par = opt$par
  out$mov = mov
  out$predfracsout = fracsin %*% mov
  out$predprobs = diag(mov)
  out
}


