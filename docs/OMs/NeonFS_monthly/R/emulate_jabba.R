# Emulate JABBA cohort assessments in the OM through recruitment deviations.
#
# For each simulation i, a posterior draw of the JABBA cohort assessment gives
# K_i and a biomass trajectory P_i[y] = B_y / K (2005-2025). The OM is made to
# follow it:
#   1. Scale:   R0 is set so that the OM's unfished mean annual biomass equals K_i.
#   2. Depletion: an annual recruitment multiplier d[y] is solved so that
#      SSB_y / SSB0 equals P_i[y], where SSB_y is spawning production summed
#      over a spawning year that contains a single cohort's spawning, relative
#      to unfished equilibrium over the same months. d[y] applies to every
#      month of that cohort's recruitment:
#        winter-spring: recruits of calendar year y - 1; SSB over Jul y-1 to Jun y
#        autumn:        recruits of Jul y - 1 to Jun y;  SSB over calendar year y
#   3. Catch:  effort is solved at the same time so that annual catches are kept.
#
# Each cohort's catch and spawning come from the same fish: winter-spring
# recruits of year y - 1 are caught by the west fleet in Jul-Nov y - 1 and spawn
# in Jul y-1 to Jun y; autumn recruits of Jul y-1 to Jun y are caught by the east
# fleet in May-Sep y and spawn in year y. The two controls are therefore updated
# as a pair, which is stable even when catch is a large share of the cohort:
#   recruitment:  d   <- d * (C* + S*) / (C + S)        (cohort total)
#   effort:       E   <- E * log(1 - phi*) / log(1 - phi), phi = C / (C + S)
# where C is catch, S = kappa * SSB_y/SSB0 the escapement in catch units
# (kappa = unfished spawning biomass of the cohort), and * marks targets.
# Catch years not linked to a target cohort use E <- E * C*/C.
#
# The multipliers are written to SRR@RecDevHist, which MSEtool uses as given
# (multiplicative scale; no truncation is applied to supplied deviations).

# --- run JABBA for one cohort (WP07 base-case settings) ------------------------
run_jabba_cohort <- function(cohort, years, catchv, index, ni = 30000, nt = 5,
                             nb = 10000, nc = 3) {
  inp <- JABBA::build_jabba(
    catch = data.frame(Year = years, Catch = catchv),
    cpue  = data.frame(Year = years, GAM = index), se = NULL,
    assessment = "NFS", scenario = cohort, model.type = "Schaefer",
    r.dist = "range", r.prior = c(0.6, 1.5),            # WP07: r in 0.6-1.5
    K.dist = "range", K.prior = c(1, 10) * max(catchv), # WP07: 1-10 x max catch
    igamma = c(0.001, 0.001),                           # WP07: diffuse process variance
    verbose = FALSE)
  JABBA::fit_jabba(inp, ni = ni, nt = nt, nb = nb, nc = nc, save.all = TRUE,
                   verbose = FALSE)
}

# --- posterior draws, one per simulation -----------------------------------------
jabba_draws <- function(fit, nSim, years, seed = 1) {
  set.seed(seed)
  po  <- fit$posteriors
  idx <- sample(nrow(po$P), nSim)
  P   <- po$P[idx, seq_along(years), drop = FALSE]
  dimnames(P) <- list(NULL, years)
  list(P = P, K = po$K[idx], r = po$r[idx], ProcDev = po$Proc.Dev[idx, , drop = FALSE],
       draw = idx)
}

# --- unfished mean annual biomass (t) per simulation and stock -----------------------
unfished_B <- function(hist, year) {
  B0 <- hist@Unfished@Equilibrium@Biomass                      # Sim x Stock x Time (g)
  yr <- floor(as.numeric(dimnames(B0)[[3]]) + 1e-6)
  apply(B0[, , yr == year, drop = FALSE], 1:2, mean) / 1e6
}

# --- the solver ---------------------------------------------------------------------
emulate_jabba <- function(x, draws, maxit = 40, tol_catch = 0.01, tol_dep = 0.02,
                          damp = 0.5, verbose = TRUE) {
  om     <- x$om
  nSim   <- om@nSim
  tgt    <- x$catch
  years  <- tgt$Year
  fleets <- c("West", "East")
  tlab   <- as.numeric(dimnames(om@Fleet[[1]][[1]]@Effort@Effort)$Year)
  tyr    <- floor(tlab + 1e-6)
  tmon   <- rep(1:12, length.out = length(tlab))
  py     <- as.integer(colnames(draws[[1]]$P))               # target years

  # recruitment window (time-step indices) feeding spawning in year y
  window <- list(WS  = function(y) which(tyr == y - 1),
                 AUT = function(y) which((tyr == y - 1 & tmon >= 7) | (tyr == y & tmon <= 6)))
  # spawning-year window (time-step indices) over which SSB_y is summed
  spwin  <- list(WS  = function(y) which((tyr == y - 1 & tmon >= 7) | (tyr == y & tmon <= 6)),
                 AUT = function(y) which(tyr == y))
  rel_win <- function(h) {
    sp  <- h@SProduction[, , seq_along(tlab), drop = FALSE]
    sp0 <- h@Unfished@Equilibrium@SProduction[, , seq_along(tlab), drop = FALSE]
    out <- array(NA, c(nSim, 2, length(py)))
    for (i in 1:2) for (k in seq_along(py)) {
      w <- spwin[[names(x$R0)[i]]](py[k])
      out[, i, k] <- apply(sp[, i, w, drop = FALSE], 1, sum) /
                     apply(sp0[, i, w, drop = FALSE], 1, sum)
    }
    out
  }

  # 1. scale R0 to K (unfished mean annual biomass in the last historical year)
  h  <- Simulate(om, silent = TRUE, control = fast_control())
  B0 <- unfished_B(h, max(years))
  for (i in 1:2) {
    s <- names(x$R0)[i]
    x$om <- om
    x <- set_R0(x, i, x$R0[[s]] * draws[[s]]$K / B0[, i])
    om <- x$om
  }

  # 2-3. joint fixed-point iteration on effort and recruitment multipliers
  dev <- lapply(1:2, function(i) matrix(1, nSim, length(tlab),
                                         dimnames = list(Sim = 1:nSim, Year = tlab)))
  names(dev) <- names(x$R0)
  fleet_of <- c(WS = "West", AUT = "East")
  catch_lag <- c(WS = 1, AUT = 0)                 # catch year = target year - lag
  allyrs <- c(min(tyr):max(years))
  trace <- NULL
  for (it in seq_len(maxit)) {
    for (i in 1:2) om@Stock[[i]]@SRR@RecDevHist <- dev[[i]]
    h  <- Simulate(om, silent = TRUE, control = fast_control())
    al <- annual_landings(h)
    rel <- rel_win(h)
    # kappa: unfished spawning biomass (t) of a cohort, mean over its spawning window
    sp0 <- h@Unfished@Equilibrium@SProduction[, , seq_along(tlab), drop = FALSE] / 1e6

    cdev <- matrix(0, nSim, 2, dimnames = list(NULL, fleets))
    ddev <- matrix(0, nSim, 2, dimnames = list(NULL, names(dev)))
    for (i in 1:2) {
      s  <- names(dev)[i]; f <- fleet_of[[s]]
      E  <- om@Fleet[[1]][[f]]@Effort@Effort
      Cp <- matrix(0, nSim, length(allyrs)); sub <- al[al$Fleet == f, ]
      Cp[cbind(sub$Sim, match(sub$Y, allyrs))] <- sub$Value
      Co <- matrix(c(rep(tgt[[f]][1], min(years) - min(allyrs)), tgt[[f]]),
                   nSim, length(allyrs), byrow = TRUE)
      cdev[, f] <- apply(abs(log((Cp + 1e-3) / (Co + 1e-3)))[, allyrs >= min(years), drop = FALSE], 1, max)
      ddev[, i] <- apply(abs(log(draws[[s]]$P / rel[, i, ])), 1, max)
      linked <- rep(FALSE, length(allyrs))
      for (k in seq_along(py)) {
        cy <- py[k] - catch_lag[[s]]; j <- match(cy, allyrs)
        if (is.na(j)) next
        linked[j] <- TRUE
        w     <- spwin[[s]](py[k])
        kappa <- apply(sp0[, i, w, drop = FALSE], 1, mean)
        S  <- kappa * rel[, i, k];  St <- kappa * draws[[s]]$P[, k]
        C  <- pmax(Cp[, j], 1e-6);  Ct <- Co[, j]
        rR <- (Ct + St) / (C + S)
        phi <- pmin(C / (C + S), 0.999); phit <- pmin(Ct / (Ct + St), 0.999)
        rE <- ifelse(Ct <= 0, 0, log(1 - phit) / log(1 - pmax(phi, 1e-9)))
        rw <- window[[s]](py[k])
        dev[[i]][, rw] <- dev[[i]][, rw] * rR^damp
        em <- which(tyr == cy)
        E[, em] <- pmin(E[, em] * pmin(pmax(rE, 0.05), 20)^damp, 10 * om@maxF)
      }
      for (j in which(!linked)) {                  # catch years without a target cohort
        em <- which(tyr == allyrs[j])
        rt <- ifelse(Co[, j] <= 0, 0, Co[, j] / pmax(Cp[, j], 1e-9))
        E[, em] <- pmin(E[, em] * pmin(rt, 20), 10 * om@maxF)
      }
      for (j2 in seq_along(om@Stock)) om@Fleet[[j2]][[f]]@Effort@Effort <- E
    }
    trace <- rbind(trace, data.frame(Iteration = it,
                                     Catch = max(cdev), WS = max(ddev[, 1]), AUT = max(ddev[, 2]),
                                     CatchMedian = median(apply(cdev, 1, max)),
                                     WSMedian = median(ddev[, 1]), AUTMedian = median(ddev[, 2])))
    if (verbose) message(sprintf("iter %2d  max |log| catch %.3f  WS %.3f  AUT %.3f",
                                 it, max(cdev), max(ddev[, 1]), max(ddev[, 2])))
    if (max(cdev) < tol_catch && max(ddev) < tol_dep) break
  }
  for (i in 1:2) om@Stock[[i]]@SRR@RecDevHist <- dev[[i]]
  x$om <- om
  x$rel <- rel                                   # last evaluated SSB_y / SSB0 (Sim x Stock x Year)
  x$recdev <- dev
  x$draws <- draws
  x$emu_trace <- trace
  x$catch_ok <- apply(cdev, 1, max) < tol_catch
  x$dep_ok   <- ddev < tol_dep
  x$converged <- x$catch_ok & x$dep_ok[, 1] & x$dep_ok[, 2]
  x
}
