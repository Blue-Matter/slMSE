# Two-cohort monthly operating model for North Pacific neon flying squid
# built with slMSE (slStock / slFleet / slCombineFleets / slOM) on MSEtool 4.
#
# Stocks : WS  = winter-spring cohort (hatch Jan-May, tail to Aug)
#          AUT = autumn cohort        (hatch Sep-Feb)
# Areas  : 1 = West of 170E, 2 = East of 170E
# Fleets : West (all members' catch west of 170E, incl. Japanese / Russian NW waters)
#          East (all members' catch east of 170E)
# Time   : monthly; ages 1-13 months; no plus group (semelparous, ~1 yr life span)

# --- small helpers -------------------------------------------------------------

wnorm <- function(mu, sigma, n = 12) MSEtool::WrappedNormal(n, mu, sigma)

logistic <- function(x, x50, slp) 1 / (1 + exp(-(x - x50) * slp))

# fleet selectivity at age (months): logistic recruitment to the jig grounds,
# declining for spawners (ages 11-13) that have left the feeding grounds
sel_at_age <- function(ages, sel50, selslp, dome50, domeslp) {
  s <- logistic(ages, sel50, selslp) * (1 - logistic(ages, dome50, domeslp))
  s / max(s)
}

# --- default settings ---------------------------------------------------------

NFS_settings <- function() {
  list(
    nages = 13,
    # recruitment (age 1 month) timing = hatching + 1 month
    WS  = list(Name = "Winter-spring cohort", rec_mu = 4,  rec_sd = 1.5,  # hatch peak Mar
               Linf = 450, K = 0.15, t0 = 0,          # ML (mm) at age (months)
               West = 0.9999,                         # fraction west of 170E (170E rule: all)
               M    = c(0.1, 0.1),                    # monthly M (1.2 / yr), fixed
               R0   = c(0.05e9, 4e9)),             # annual unfished recruits, log-U(lower, upper)
    AUT = list(Name = "Autumn cohort",        rec_mu = 12, rec_sd = 1.5,  # hatch peak Nov
               Linf = 550, K = 0.15, t0 = 0,
               West = 0.0001,                         # (170E rule: none)
               M    = c(0.1, 0.1),                    # monthly M (1.2 / yr), fixed
               R0   = c(0.1e9, 2e9)),
    a = 2.60e-5, b = 3.021,                           # W(g) = a ML(mm)^b
    amat50 = 9, amatSLP = 1.5,                        # maturity at age (months)
    h = 0.8, sigmaR = 0.6, R_AC = 0.8,                # Bev-Holt, monthly rec devs
    prob_stay = 0.95,
    # fleets: seasonal effort profile (wrapped normal over calendar months)
    West = list(mu = 9,   sd = 1.3),                  # Jul-Nov, peak Sep (WP04-WP06)
    East = list(mu = 7.5, sd = 1.3),                  # May-Sep, peak Aug (WP03)
    sel  = c(sel50 = 4, selslp = 1.5, dome50 = 11.5, domeslp = 2),
    maxF = 2                                          # monthly F cap
  )
}

# --- stock --------------------------------------------------------------------

make_cohort <- function(cs, set, nSim, nYear, pYear, CurrentYear) {
  nages <- set$nages
  # numbers-at-age area split (constant with age for now)
  Frac <- rbind(rep(cs$West, nages), rep(1 - cs$West, nages))
  st <- slStock(Name = cs$Name, Species = "Ommastrephes bartramii",
                CommonName = "Neon flying squid",
                nYear = nYear, pYear = pYear, Seasons = 12,
                CurrentYear = CurrentYear, nSim = nSim,
                rec_age = 1, nages = nages, PlusGroup = FALSE,
                spawndist = wnorm(cs$rec_mu, cs$rec_sd),
                Linf = cs$Linf, K = cs$K, t0 = cs$t0, Len_CV = 0.1,
                lenclasses = seq(10, 700, 10),
                a = set$a, b = set$b,
                M = mean(cs$M),                      # replaced below by per-sim M
                amat50 = set$amat50, amatSLP = set$amatSLP,
                SR_type = "BevertonHolt", h = set$h, sigmaR = set$sigmaR,
                trunc_sigmaR = 2, R_AC = set$R_AC, R0 = 1,  # R0 scaled below
                nareas = 2, Frac_area = Frac, prob_stay = set$prob_stay)
  st
}

# --- fleet --------------------------------------------------------------------

make_fleet <- function(name, area, set, nSim, nYear, pYear, CurrentYear) {
  fl <- slFleet(Name = name, nYear = nYear, pYear = pYear, Seasons = 12,
                CurrentYear = CurrentYear, nSim = nSim, rec_age = 1,
                nAges = set$nages, maxF = 0.1)       # placeholder effort, replaced later
  ages <- 1:set$nages
  s <- set$sel
  fl@Selectivity <- Selectivity(MeanAtAge = sel_at_age(ages, s["sel50"], s["selslp"],
                                                       s["dome50"], s["domeslp"]))
  fl
}

# --- operating model ----------------------------------------------------------

#' @param catch data.frame Year, West, East (tonnes), from NFS_catch_by_area()
#' @param burnin number of years before the first catch year, fished at the
#'   first-year catch, so that 2005 does not start from unfished conditions
make_NFS_OM <- function(catch, nSim = 24, pYear = 10, burnin = 2,
                        set = NFS_settings(), seed = 1) {

  set.seed(seed)
  catch       <- catch[order(catch$Year), ]
  firstY      <- min(catch$Year) - burnin
  CurrentYear <- max(catch$Year)
  nYear       <- CurrentYear - firstY + 1

  stocks <- list(WS  = make_cohort(set$WS,  set, nSim, nYear, pYear, CurrentYear),
                 AUT = make_cohort(set$AUT, set, nSim, nYear, pYear, CurrentYear))

  fleets <- list(West = make_fleet("West", 1, set, nSim, nYear, pYear, CurrentYear),
                 East = make_fleet("East", 2, set, nSim, nYear, pYear, CurrentYear))
  fleet  <- slCombineFleets(fleets, stocks)

  om <- slOM(Name = "NPFC neon flying squid: two-cohort monthly OM",
             Agency = "NPFC", Author = "Blue Matter Science", Region = "North Pacific",
             ComplexName = "Neon flying squid",
             nSim = nSim, nYear = nYear, pYear = pYear, Seasons = 12,
             CurrentYear = CurrentYear, Interval = 1, Seed = seed,
             stock = stocks, fleet = fleet)
  om@maxF <- set$maxF

  # --- per-simulation parameters ----------------------------------------------
  # M: one uniform quantile per simulation, mapped onto each cohort's own range
  # (so the cohorts' M are perfectly correlated across simulations)
  uM <- runif(nSim)
  Mv <- list(WS  = set$WS$M[1]  + uM * diff(set$WS$M),
             AUT = set$AUT$M[1] + uM * diff(set$AUT$M))
  rlu <- function(r) exp(runif(nSim, log(r[1]), log(r[2])))   # log-uniform
  R0 <- list(WS = rlu(set$WS$R0), AUT = rlu(set$AUT$R0))

  nT <- (nYear + pYear) * 12
  for (i in seq_along(stocks)) {
    s    <- names(stocks)[i]
    st   <- om@Stock[[i]]
    na   <- set$nages
    Marr <- array(Mv[[s]], c(nSim, na, 1),
                  dimnames = list(Sim = 1:nSim,
                                  Age = dimnames(st@NaturalMortality@MeanAtAge)$Age,
                                  Year = firstY))
    st@NaturalMortality@MeanAtAge <- Marr

    pi_m <- wnorm(set[[s]]$rec_mu, set[[s]]$rec_sd)
    R0a  <- outer(R0[[s]], rep(pi_m, nYear + pYear))
    dimnames(R0a) <- dimnames(st@SRR@R0)
    if (is.null(dimnames(R0a))) dimnames(R0a) <- list(Sim = 1:nSim, Year = NULL)
    st@SRR@R0 <- R0a
    om@Stock[[i]] <- st
  }

  # --- fleets: area, seasonal effort, placeholder annual effort ----------------
  tlab  <- dimnames(om@Fleet[[1]][[1]]@Effort@Effort)$Year
  month <- rep(1:12, nYear)
  yr    <- rep(firstY:CurrentYear, each = 12)
  for (f in c("West", "East")) {
    prof <- wnorm(set[[f]]$mu, set[[f]]$sd)
    prof <- prof / max(prof)
    Eff  <- matrix(rep(prof[month], each = nSim), nSim, nYear * 12,
                   dimnames = list(Sim = 1:nSim, Year = tlab)) * 0.05
    area <- if (f == "West") c(1, 0) else c(0, 1)
    Dist <- array(area, c(1, 1, 2),
                  dimnames = list(Sim = 1, Year = tlab[1], Area = 1:2))
    for (i in seq_along(stocks)) {
      om@Fleet[[i]][[f]]@Effort@Effort       <- Eff
      om@Fleet[[i]][[f]]@Effort@Distribution <- Dist
      om@Fleet[[i]][[f]]@Effort@Mode         <- "Biomass"
    }
  }

  # catch targets incl. burn-in years (fished at first-year catch)
  tgt <- rbind(data.frame(Year = firstY:(min(catch$Year) - 1),
                          West = catch$West[1], East = catch$East[1])[seq_len(burnin), ],
               catch[, c("Year", "West", "East")])

  list(om = om, set = set, M = Mv, R0 = R0, catch = tgt, firstY = firstY,
       month = month, yr = yr)
}
