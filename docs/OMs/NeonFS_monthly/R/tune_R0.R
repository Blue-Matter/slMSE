# R0 helpers and per-simulation R0 tuning with the JABBA-emulated recruitment
# deviations held fixed.

# set annual R0 (vector, one value per simulation) for stock i, keeping the
# seasonal recruitment pattern
set_R0 <- function(x, i, R0) {
  s    <- names(x$R0)[i]
  pi_m <- wnorm(x$set[[s]]$rec_mu, x$set[[s]]$rec_sd)
  R0a  <- x$om@Stock[[i]]@SRR@R0
  R0a[] <- outer(R0, rep(pi_m, ncol(R0a) / 12))
  x$om@Stock[[i]]@SRR@R0 <- R0a
  x$R0[[s]] <- R0
  x
}

# SSB_y / SSB0 over spawning-year windows that contain a single cohort's
# spawning (as in emulate_jabba()):
#   winter-spring: Jul y-1 to Jun y;  autumn: calendar year y.
# Returns Sim x Stock x Year.
spawn_rel <- function(h, tlab, years = 2005:2025, stocks = c("WS", "AUT")) {
  tyr  <- floor(tlab + 1e-6)
  tmon <- rep(1:12, length.out = length(tlab))
  win  <- list(WS  = function(y) which((tyr == y - 1 & tmon >= 7) | (tyr == y & tmon <= 6)),
               AUT = function(y) which(tyr == y))
  sp   <- h@SProduction[, , seq_along(tlab), drop = FALSE]
  sp0  <- h@Unfished@Equilibrium@SProduction[, , seq_along(tlab), drop = FALSE]
  out  <- array(NA, c(dim(sp)[1], length(stocks), length(years)),
                dimnames = list(Sim = NULL, Stock = stocks, Year = years))
  for (i in seq_along(stocks)) for (k in seq_along(years)) {
    w <- win[[stocks[i]]](years[k])
    out[, i, k] <- apply(sp[, i, w, drop = FALSE], 1, sum) /
                   apply(sp0[, i, w, drop = FALSE], 1, sum)
  }
  out
}

# Tune R0 per simulation and cohort so that SSB2025/SSB0 equals that
# simulation's JABBA B2025/K draw, holding the emulated recruitment deviations
# (SRR@RecDevHist) fixed and re-solving effort at every step to keep the
# catches. Bisection on log(R0) within a multiplicative range around the
# K-based R0 of the emulation. Lower R0 -> lower depletion; an R0 at which a
# cohort's own fleet cannot take its catches is treated as too low.
tune_R0_jabba <- function(xj, year = 2025, range = c(0.7, 1.4), nstep = 8,
                          verbose = TRUE) {
  x      <- xj
  nSim   <- x$om@nSim
  tlab   <- as.numeric(dimnames(x$om@Fleet[[1]][[1]]@Effort@Effort)$Year)
  fleet_of <- c(WS = "West", AUT = "East")
  E0 <- lapply(c(West = "West", East = "East"),
               function(f) x$om@Fleet[[1]][[f]]@Effort@Effort)
  reset_effort <- function(x) {             # restart from the emulation's solved effort
    for (f in names(E0)) for (i in seq_along(x$om@Stock))
      x$om@Fleet[[i]][[f]]@Effort@Effort <- E0[[f]]
    x
  }
  x$R0_K <- x$R0                            # K-based R0 from the emulation
  for (s in c("WS", "AUT")) {
    i   <- match(s, names(x$R0))
    tgt <- x$draws[[s]]$P[, as.character(year)]
    lo  <- log(x$R0_K[[s]] * range[1]); hi <- log(x$R0_K[[s]] * range[2])
    for (k in seq_len(nstep)) {
      x <- set_R0(x, i, exp((lo + hi) / 2))
      x <- solve_effort(reset_effort(x), maxit = 12, verbose = FALSE)
      h <- Simulate(x$om, silent = TRUE, control = fast_control())
      d <- spawn_rel(h, tlab, year)[, i, 1]
      too_high <- x$fleet_ok[, fleet_of[[s]]] & d > tgt
      hi[too_high]  <- log(x$R0[[s]][too_high])
      lo[!too_high] <- log(x$R0[[s]][!too_high])
      if (verbose) message(sprintf("%s step %d  feasible %d/%d  median |log(dep/target)| %.3f",
                                   s, k, sum(x$fleet_ok[, fleet_of[[s]]]), nSim,
                                   median(abs(log(d / tgt)))))
    }
    x <- set_R0(x, i, exp((lo + hi) / 2))
  }
  x <- solve_effort(reset_effort(x), maxit = 15, verbose = FALSE)
  h <- Simulate(x$om, silent = TRUE, control = fast_control())
  x$rel <- spawn_rel(h, tlab, as.integer(colnames(x$draws$WS$P)))
  fallback_to_K(x, xj, year = year)
}

# Simulations whose tuned SSB_year/SSB0 misses the target by more than `tol`
# (log scale) keep the K-based R0 and the solved effort of the emulation, which
# match the JABBA draw by construction. This happens where the catch takes most
# of a cohort: depletion is then very sensitive to R0, and the effort solver
# converges slowly enough at lower R0 for the bisection to misjudge feasibility.
fallback_to_K <- function(x, xj, year = 2025, tol = 0.05) {
  tlab <- as.numeric(dimnames(x$om@Fleet[[1]][[1]]@Effort@Effort)$Year)
  bad <- lapply(c(WS = "WS", AUT = "AUT"), function(s)
    which(abs(log(x$rel[, s, as.character(year)] / x$draws[[s]]$P[, as.character(year)])) > tol))
  x$fallback <- bad
  sims <- sort(unique(unlist(bad)))
  if (length(sims) == 0) return(x)
  for (i in 1:2) {
    s <- names(x$R0)[i]
    R0 <- x$R0[[s]]; R0[bad[[s]]] <- x$R0_K[[s]][bad[[s]]]
    x <- set_R0(x, i, R0)
  }
  for (f in c("West", "East")) {
    E <- x$om@Fleet[[1]][[f]]@Effort@Effort
    E[sims, ] <- xj$om@Fleet[[1]][[f]]@Effort@Effort[sims, ]
    for (j in seq_along(x$om@Stock)) x$om@Fleet[[j]][[f]]@Effort@Effort <- E
  }
  x <- solve_effort(x, maxit = 15, verbose = FALSE)
  h <- Simulate(x$om, silent = TRUE, control = fast_control())
  x$rel <- spawn_rel(h, tlab, as.integer(colnames(x$draws$WS$P)))
  x
}
