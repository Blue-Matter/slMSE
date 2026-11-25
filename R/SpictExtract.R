#' Extract biomass and catch estimates from a SPiCT fit by year or quarter
#'
#' Point estimates only (no Hessian / sdreport calculations).
#'
#' Biomass:
#'   bio_type = "start": B at the start of each interval (e.g. B(2010.0), B(2010.25)),
#'                       from the logB random effect (as spict::get.par("logB", fit, exp=TRUE)).
#'   bio_type = "mean" : arithmetic mean of B over the Euler steps in the interval.
#' Catch:
#'   Predicted catch over [t_start, t_end), computed as in spict.cpp:
#'   sum over Euler steps i of F_i * B_i * dt_i, with F = exp(logFs) (incl. seasonality).
#'
#' @param fit        A fitted 'spictcls' object (from spict::fit.spict).
#' @param resolution "year" or "quarter".
#' @param years      Integer vector of years. Default: all years with catch observations.
#'                   Only whole intervals inside the catch observation period are returned.
#' @param bio_type   "start" or "mean".
#' @return data.frame, one row per interval.
spict_extract <- function(fit, resolution = c("year", "quarter"), years = NULL,
                          bio_type = c("start", "mean")) {
  stopifnot(inherits(fit, "spictcls"))
  resolution <- match.arg(resolution)
  bio_type   <- match.arg(bio_type)

  inp  <- fit$inp
  tg   <- inp$time                 # Euler time grid (state times)
  dte  <- inp$dteuler
  ns   <- length(tg)
  step <- if (resolution == "year") 1 else 0.25
  tol  <- dte / 10

  if (abs(round(step / dte) * dte - step) > tol)
    stop("Interval length ", step, " is not a multiple of inp$dteuler (", dte, ").")

  # ---- intervals ---------------------------------------------------------
  if (is.null(years)) years <- seq(floor(min(inp$timeC) + tol), floor(max(inp$timeC) + tol))
  t0 <- if (resolution == "year") as.numeric(years) else
    as.vector(t(outer(as.numeric(years), c(0, 0.25, 0.5, 0.75), "+")))
  t1 <- t0 + step
  # restrict to the period covered by the catch observations (inp$obsC):
  # from the first catch time to the end of the last catch interval
  iobs    <- seq_len(inp$nobsC)
  obs_rng <- c(min(inp$timeC), max(inp$timeC + inp$dtcp[iobs]))
  keep <- t0 >= max(min(tg), obs_rng[1]) - tol & t1 <= min(max(tg), obs_rng[2]) + tol
  if (any(!keep))
    message("Dropping ", sum(!keep), " interval(s) outside the catch observation period [",
            round(obs_rng[1], 3), ", ", round(obs_rng[2], 3), "].")
  t0 <- t0[keep]; t1 <- t1[keep]
  if (!length(t0)) stop("No whole intervals inside the catch observation period.")

  match_grid <- function(t) {
    i <- vapply(t, function(x) { d <- abs(tg - x); j <- which.min(d)
    if (d[j] > tol) NA_integer_ else j }, 1L)
    if (anyNA(i)) stop("Interval boundaries not on the Euler grid: ",
                       paste(t[is.na(i)], collapse = ", "))
    i
  }
  ic  <- match_grid(t0)                          # index of first state in interval
  nc  <- as.integer(round((t1 - t0) / dte))      # number of Euler steps in interval
  idx <- lapply(seq_along(ic), function(k) ic[k] + seq_len(nc[k]) - 1L)

  # ---- estimates from the fitted object ------------------------------------
  logB <- fit$par.random[names(fit$par.random) == "logB"]
  logFs <- fit$value[names(fit$value) == "logFs"]
  stopifnot(length(logB) == ns, length(logFs) == ns)
  B <- exp(logB)

  # catch per Euler step, as in spict.cpp
  Csub <- exp(logFs[-ns] + logB[-ns]) * inp$dt[seq_len(ns - 1)]
  sum_steps <- function(ic, nc) vapply(seq_along(ic), function(k)
    sum(Csub[ic[k] + seq_len(nc[k]) - 1L]), 1)

  # sanity check against the fitted catch predictions
  orig <- fit$value[names(fit$value) == "logCpred"]
  chk  <- max(abs(log(sum_steps(inp$ic, inp$nc)) - orig))
  if (chk > 1e-6) warning("Catch reconstruction differs from fitted logCpred (max |diff| on log scale ",
                          signif(chk, 3), ").")

  B_est <- if (bio_type == "start") B[ic] else vapply(idx, function(j) mean(B[j]), 1)
  C_est <- sum_steps(ic, nc)

  # ---- observed catch summed over interval (only where fully covered) ------
  tC <- inp$timeC; dC <- inp$dtcp[seq_len(inp$nobsC)]; oC <- inp$obsC
  C_obs <- vapply(seq_along(t0), function(k) {
    inside <- tC >= t0[k] - tol & tC < t1[k] - tol
    if (!any(inside)) return(NA_real_)
    if (abs(sum(dC[inside]) - (t1[k] - t0[k])) > 2 * dte) return(NA_real_)  # partial coverage
    sum(oC[inside])
  }, 1)

  year <- floor(t0 + tol)
  out <- data.frame(time = t0, year = year)
  if (resolution == "quarter") out$quarter <- as.integer(round((t0 - year) / 0.25)) + 1L
  out$t_start <- t0
  out$t_end   <- t1
  out$B       <- B_est
  out$C       <- C_est
  out$C_obs   <- C_obs
  attr(out, "bio_type") <- bio_type
  out
}
