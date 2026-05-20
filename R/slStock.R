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
#' @param Len_age Vector of positive real numbers, nages long. The Length at age.
#' @param Len_CV Positive real number. The coefficient of variation of length at age. Typically between 0.05 and 0.25.
#' @param a Positive real number. Weight (W) at length (L) parameter a. W = aL^b. Converts units of length to weight.
#' @param b Positive real number. Weight (W) at length (L) parameter b. W = aL^b. Typically approximately cubic (3-ish)
#' @param Wt_age Vector of positive real numbers, nages long. The weight at age.
#' @param M Positive real number or vector of positive real numbers nSim long. The instantaneous rate of natural mortality per season.
#' @param amat50 Positive real number or vector of positive real numbers nSim long. The age (in seasons) that 50% of individuals are mature, in a logistic maturity model. E.g., 6 in a monthly model would be 50% mature after 6 months.
#' @param amatSLP Positive real number or vector of positive real numbers nSim long. The slope of the logistic maturity model. E.g., 2 in a monthly model.
#' @param Mat_age Vector of positive real numbers, nages long. The spawning fraction ('maturity') at age.
#' @param SR_type Character string. The type of stock recruitment relationship e.g., BevertonHolt', 'Ricker'
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
                   Linf = 1, K = 0.2, t0 = 0, Len_CV = 0.2, Len_age = NA,
                   a = 1E-5, b = 3, Wt_age = NA,
                   M = 0.2, amat50 = 6, amatSLP = 2,
                   Mat_age = NA,
                   SR_type = "BevertonHolt", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                   nareas = 2, Frac_area = matrix(c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.5,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.01,0.01,
                                                    0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.5,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99,0.99,0.99),byrow=T,nrow=2),
                   prob_stay = 0.9){

  # Default args for testing:
  # Name = "A short-lived creature"; Species = "Shortus liveus"; CommonName = "Short-lived creature"
  # nYear = 10; pYear = 10; Seasons = 12; CurrentYear = 2026; nSim = 4; rec_age = 1; nages = 24; PlusGroup = F
  # Linf = 1; K = 0.2; t0 = 0; Len_CV = 0.2; a = 1; b = 3; M = 0.2; amat50 = 6; amatSLP = 2;
  # h = 0.9; R0 = 1E6; sigmaR = 1.0; trunc_sigmaR = 2.0; R_AC = 0.5; nareas = 2;  prob_stay = 0.9
  # Frac_area = matrix(c(0.95,0.9,0.85,0.8,0.75,0.7,0.65,0.6,0.55,0.5,0.5,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05,0.01,0.01,0.01, 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.5,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99,0.99,0.99),nrow=2,byrow=T)
  # Len_age = NA; Wt_age = NA; Mat_age = NA; SR_type = "BevertonHolt"
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
  if(!is.na(Len_age[1]))  Length@MeanAtAge = Len_age
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
  if(!is.na(Wt_age[1]))  Length@MeanAtAge = Wt_age
  Weight@Units = "g"
  Weight@Dist = "lognormal"
  Weight@TruncSD = 2
  stock@Weight = Weight

  # Natural Mortality
  stock@NaturalMortality = NaturalMortality(MeanAtAge = rep(M,na))

  # Maturity
  avec = ((rec_age : nages)-amat50) * amatSLP
  atrans = exp(avec)/(1+exp(avec))
  if(!is.na(Mat_age[1])) atrans = Mat_age
  amat = array(rep(atrans,each=nSim),c(nSim,na,nYear*Seasons+pYear*Seasons))
  stock@Maturity = Maturity(MeanAtAge = amat)
  #stock@Maturity = Maturity(Pars = list(L50 = amat50, L50_95 = amat50/10))

  # Fecundity slot for spawndist


  # Stock-Recruitment
  stock@SRR = SRR(Model = SR_type, Pars = list(h = h), R0 = R0, SD = sigmaR, AC = R_AC, TruncSD = trunc_sigmaR)

  # Spatial
  # nSim, nArea, nArea, nAge, and nTS,
  Movement = array(0,c(1,nareas,nareas,na,1))
  dist1 = Frac_area[,1]
  Movement[1,,,1,1] = matrix(dist1,ncol=nareas,nrow=nareas,byrow=T) # Initial movement is fully mixed
  for(aa in 2:na){
    fracsin = Frac_area[,aa-1]
    fracsout = Frac_area[,aa]
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


# === Stock Internal ================================================================================


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






