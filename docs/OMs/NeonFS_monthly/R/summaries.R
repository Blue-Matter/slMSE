# Annual summaries from a monthly hist object ---------------------------------

# Sim x Stock x Year array of annual sums of a monthly Sim x Stock x Time array
annual_sum <- function(arr) {
  yr <- floor(as.numeric(dimnames(arr)$Year) + 1e-6)
  out <- apply(arr, 1:2, function(v) tapply(v, yr, sum))   # Year x Sim x Stock
  aperm(out, c(2, 3, 1))
}

# Annual spawning production relative to unfished equilibrium (by cohort)
# type = "Dynamic" uses the unfished trajectory with the same recruitment
# deviations (isolates the effect of fishing)
SP_rel <- function(hist, type = c("Dynamic", "Equilibrium")) {
  type <- match.arg(type)
  sp  <- annual_sum(hist@SProduction)
  sp0 <- annual_sum(slot(hist@Unfished, type)@SProduction)
  yrs <- intersect(dimnames(sp)[[3]], dimnames(sp0)[[3]])
  sp[, , yrs, drop = FALSE] / sp0[, , yrs, drop = FALSE]
}
