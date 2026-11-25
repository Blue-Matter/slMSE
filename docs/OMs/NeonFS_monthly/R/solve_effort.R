# Solve monthly effort so that each simulation reproduces the observed annual
# catch by fleet (West / East of 170E). The within-year seasonal effort profile
# of each fleet is kept fixed; only its annual level is scaled.
#
# Fixed-point iteration: E[f, y] <- E[f, y] * (Cobs[f, y] / Cpred[f, y])^damp.
# Each fleet's catch depends mostly on its own effort (the fleets fish
# different areas), so this converges in a handful of iterations unless the
# catch is infeasible (F hits OM@maxF).

fast_control <- function()
  SimControl(DynamicUnfished = FALSE, ConditionObs = FALSE, GenerateData = FALSE,
             MSYRefs = FALSE, RefPoints = FALSE, MGT = FALSE, CalcCatchAtSize = FALSE)

# annual landings (tonnes) by Sim x Year x Fleet, summed over stocks
annual_landings <- function(hist) {
  L <- Landings(hist)
  L$Y <- floor(L$Year + 1e-6)
  L <- subset(L, Period == "Historical")
  a <- aggregate(Value ~ Sim + Y + Fleet, data = L, sum)
  a$Value <- a$Value / 1e6                       # g -> t
  a
}

# record = TRUE keeps the annual catch predicted at each iteration and a
# convergence trace (x$trace_landings, x$trace), for diagnostics.
solve_effort <- function(x, maxit = 15, tol = 0.01, damp = 1, verbose = TRUE,
                         record = FALSE) {
  om    <- x$om
  tgt   <- x$catch
  nSim  <- om@nSim
  years <- tgt$Year
  yr    <- x$yr
  fleets <- c("West", "East")

  simdev <- rep(Inf, nSim)
  trace  <- NULL
  tl     <- NULL
  for (it in seq_len(maxit)) {
    h  <- Simulate(om, silent = TRUE, control = fast_control())
    al <- annual_landings(h)
    if (record) tl <- rbind(tl, data.frame(Iteration = it, al))
    maxdev <- 0
    simdev <- rep(0, nSim)
    fltdev <- matrix(0, nSim, 2, dimnames = list(NULL, fleets))
    for (f in fleets) {
      E    <- om@Fleet[[1]][[f]]@Effort@Effort
      pred <- matrix(0, nSim, length(years))
      sub  <- al[al$Fleet == f, ]
      pred[cbind(sub$Sim, match(sub$Y, years))] <- sub$Value
      obs  <- matrix(tgt[[f]], nSim, length(years), byrow = TRUE)
      rat  <- ifelse(obs == 0, 0, (obs + 1e-9) / pmax(pred, 1e-9))
      rat  <- pmin(rat, 20)^damp
      dev  <- abs(log((pred + 1e-3) / (obs + 1e-3))) * (obs > 0)
      fltdev[, f] <- apply(dev, 1, max)
      simdev <- pmax(simdev, fltdev[, f])
      maxdev <- max(maxdev, dev)
      E <- pmin(E * rat[, match(yr, years)], 10 * om@maxF)   # cap runaway effort
      for (i in seq_along(om@Stock)) om@Fleet[[i]][[f]]@Effort@Effort <- E
    }
    trace <- rbind(trace, data.frame(Iteration = it, MaxAbsLogDev = maxdev,
                                     MedianSimMaxDev = median(simdev),
                                     SimsWithinTol = sum(simdev < tol)))
    if (verbose) message(sprintf("iter %2d  max |log(pred/obs)| = %.4f  (%d/%d sims within tol)",
                                 it, maxdev, sum(simdev < tol), nSim))
    if (maxdev < tol) break
  }
  x$om <- om
  x$converged <- simdev < tol          # per simulation: catches matched
  x$fleet_ok  <- fltdev < tol          # per simulation and fleet
  x$simdev    <- simdev
  x$maxdev    <- maxdev
  x$trace     <- trace
  if (record) x$trace_landings <- tl
  x
}
