# Figures for the two-cohort monthly OM ---------------------------------------
library(ggplot2)

COH_COL  <- c("Winter-spring cohort" = "#2a78d6", "Autumn cohort" = "#eb6834")
INK      <- "#3d3d3a"
theme_nfs <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e6e6e3", linewidth = 0.3),
          strip.text = element_text(face = "bold", hjust = 0),
          legend.position = "top", legend.title = element_blank(),
          axis.title = element_text(colour = INK))
}

q_df <- function(df, by) {               # median and 90% interval
  do.call(data.frame, aggregate(df$Value, df[by], function(v)
    c(lo = unname(quantile(v, 0.05)), med = median(v), hi = unname(quantile(v, 0.95)))))
}

# 1. Observed catch by area -----------------------------------------------------
plot_catch_area <- function(ca) {
  d <- rbind(data.frame(Year = ca$Year, Area = "West of 170°E", Catch = ca$West),
             data.frame(Year = ca$Year, Area = "East of 170°E", Catch = ca$East))
  d$Area <- factor(d$Area, c("West of 170°E", "East of 170°E"))
  ggplot(d, aes(Year, Catch / 1000)) +
    geom_col(fill = INK, width = 0.75) +
    facet_wrap(~Area, ncol = 2) +
    labs(x = NULL, y = "Catch (kt)") + theme_nfs()
}

# 2. Schedules ------------------------------------------------------------------
plot_schedules <- function(x) {
  set <- x$set; ages <- 1:set$nages
  m <- month.abb
  rec <- rbind(data.frame(Month = 1:12, Cohort = "Winter-spring cohort",
                          v = wnorm(set$WS$rec_mu, set$WS$rec_sd)),
               data.frame(Month = 1:12, Cohort = "Autumn cohort",
                          v = wnorm(set$AUT$rec_mu, set$AUT$rec_sd)))
  eff <- rbind(data.frame(Month = 1:12, Fleet = "West fleet", v = wnorm(set$West$mu, set$West$sd)),
               data.frame(Month = 1:12, Fleet = "East fleet", v = wnorm(set$East$mu, set$East$sd)))
  len <- rbind(data.frame(Age = ages, Cohort = "Winter-spring cohort",
                          v = set$WS$Linf * (1 - exp(-set$WS$K * (ages - set$WS$t0)))),
               data.frame(Age = ages, Cohort = "Autumn cohort",
                          v = set$AUT$Linf * (1 - exp(-set$AUT$K * (ages - set$AUT$t0)))))
  s <- set$sel
  sm <- rbind(data.frame(Age = ages, What = "Selectivity",
                         v = sel_at_age(ages, s["sel50"], s["selslp"], s["dome50"], s["domeslp"])),
              data.frame(Age = ages, What = "Maturity",
                         v = logistic(ages, set$amat50, set$amatSLP)))
  p1 <- ggplot(rec, aes(Month, v, colour = Cohort)) + geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) + scale_colour_manual(values = COH_COL) +
    scale_x_continuous(breaks = 1:12, labels = substr(m, 1, 1)) +
    labs(title = "Recruitment (age 1 month) by month", x = NULL, y = "Fraction of annual R0") +
    theme_nfs()
  p2 <- ggplot(eff, aes(Month, v)) + geom_col(fill = INK, width = 0.7) +
    facet_wrap(~Fleet) + scale_x_continuous(breaks = 1:12, labels = substr(m, 1, 1)) +
    labs(title = "Seasonal effort profile", x = NULL, y = "Fraction of annual effort") + theme_nfs()
  p3 <- ggplot(len, aes(Age, v, colour = Cohort)) + geom_line(linewidth = 0.8) +
    scale_colour_manual(values = COH_COL) + scale_x_continuous(breaks = ages) +
    labs(title = "Mantle length at age", x = "Age (months)", y = "ML (mm)") +
    theme_nfs() + theme(legend.position = "none")
  p4 <- ggplot(sm, aes(Age, v, linetype = What)) + geom_line(linewidth = 0.8, colour = INK) +
    scale_x_continuous(breaks = ages) +
    labs(title = "Selectivity (both fleets) and maturity", x = "Age (months)", y = NULL) + theme_nfs()
  patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2)
}

# 3. Monthly biomass by cohort ---------------------------------------------------
plot_biomass <- function(h, keep, from = 2005) {
  B <- as.data.frame(Biomass(h)); B <- subset(B, Sim %in% keep & Year >= from & Period == "Historical")
  B$Value <- B$Value / 1e9                           # g -> kt
  q <- q_df(B, c("Stock", "Year"))
  ggplot(q, aes(Year)) +
    geom_ribbon(aes(ymin = x.lo, ymax = x.hi, fill = Stock), alpha = 0.2) +
    geom_line(aes(y = x.med, colour = Stock), linewidth = 0.5) +
    scale_colour_manual(values = COH_COL) + scale_fill_manual(values = COH_COL) +
    facet_wrap(~Stock, ncol = 1, scales = "free_y") +
    labs(x = NULL, y = "Biomass (kt), median and 90% interval") +
    theme_nfs() + theme(legend.position = "none")
}

# 5. Spawning production relative to unfished --------------------------------------
# type = "Equilibrium": SSB / SSB0 (comparable to JABBA B/K)
# type = "Dynamic"    : relative to the unfished trajectory with the same rec devs
plot_depletion <- function(h, keep, from = 2005, to = 2025, type = "Equilibrium",
                           target = NULL) {
  r <- SP_rel(h, type)[keep, , as.character(from:to), drop = FALSE]
  d <- as.data.frame.table(r, responseName = "Value")
  names(d)[1:3] <- c("Sim", "Stock", "Year")
  d$Year <- as.numeric(as.character(d$Year))
  q <- q_df(d, c("Stock", "Year"))
  ggplot(q, aes(Year)) +
    geom_hline(yintercept = 1, colour = "#b0b0ab", linewidth = 0.3) +
    geom_ribbon(aes(ymin = x.lo, ymax = x.hi, fill = Stock), alpha = 0.2) +
    geom_line(aes(y = x.med, colour = Stock), linewidth = 0.8) +
    scale_colour_manual(values = COH_COL) + scale_fill_manual(values = COH_COL) +
    { if (!is.null(target))
        geom_point(data = data.frame(Year = to, Stock = names(COH_COL), v = unlist(target)),
                   aes(Year, v, colour = Stock), shape = 4, size = 4, stroke = 1.4) } +
    coord_cartesian(ylim = c(0, max(1.1, q$x.hi))) +
    labs(x = NULL, y = if (type == "Equilibrium") "SSB / SSB0" else
           "SSB / dynamic unfished SSB") + theme_nfs()
}

# 6. Monthly apical F by fleet (median) ----------------------------------------------
apical_F <- function(h, keep) {
  Fa <- h@FDeadArea[[1]]                          # Sim x Age x Time x Fleet x Area
  ap <- apply(Fa[keep, , , , , drop = FALSE], c(1, 3, 4), max)  # max over age & area
  med <- apply(ap, 2:3, median)
  tt <- as.numeric(dimnames(Fa)[[3]])
  d <- data.frame(Year = floor(tt + 1e-6), Month = rep(1:12, length.out = length(tt)),
                  West = med[, 1], East = med[, 2])
  rbind(data.frame(d[, 1:2], Fleet = "West fleet", F = d$West),
        data.frame(d[, 1:2], Fleet = "East fleet", F = d$East))
}
plot_F <- function(h, keep, from = 2005) {
  d <- subset(apical_F(h, keep), Year >= from)
  d$Fleet <- factor(d$Fleet, c("West fleet", "East fleet"))
  ggplot(d, aes(Month, Year, fill = F)) + geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "#f1f5fb", high = "#0d3b7a", name = "Monthly apical F",
                        guide = guide_colourbar(barwidth = unit(6, "cm"), barheight = unit(0.3, "cm"))) +
    scale_x_continuous(breaks = 1:12, labels = substr(month.abb, 1, 1), expand = c(0, 0)) +
    scale_y_reverse(expand = c(0, 0)) + facet_wrap(~Fleet) +
    labs(x = NULL, y = NULL) + theme_nfs() +
    theme(legend.title = element_text(size = 9), panel.grid = element_blank())
}


# 7. Catch fit -------------------------------------------------------------------------
plot_catch_fit <- function(h, keep, catch) {
  al <- subset(annual_landings(h), Sim %in% keep & Y >= min(catch$Year))
  al$Fleet <- factor(paste(al$Fleet, "fleet"), c("West fleet", "East fleet"))
  ob <- rbind(data.frame(Y = catch$Year, Fleet = "West fleet", Value = catch$West),
              data.frame(Y = catch$Year, Fleet = "East fleet", Value = catch$East))
  ob$Fleet <- factor(ob$Fleet, levels(al$Fleet))
  ggplot(al, aes(Y, Value / 1000)) +
    geom_line(aes(group = Sim), colour = "#b0b0ab", linewidth = 0.3) +
    geom_point(data = ob, colour = INK, size = 1.8) +
    facet_wrap(~Fleet, scales = "free_y") +
    labs(x = NULL, y = "Catch (kt)") + theme_nfs()
}

# 10. Effort solver: catch before and after solving --------------------------------
plot_effort_demo <- function(demo, catch) {
  tl <- demo$trace_landings
  last <- max(tl$Iteration)
  d <- subset(tl, Iteration %in% c(1, last) & Y >= min(catch$Year))
  d$Stage <- factor(ifelse(d$Iteration == 1, "Initial effort (iteration 1)",
                           paste0("Solved effort (iteration ", last, ")")))
  d$Stage <- factor(d$Stage, levels = c("Initial effort (iteration 1)",
                                        paste0("Solved effort (iteration ", last, ")")))
  d$Fleet <- factor(paste(d$Fleet, "fleet"), c("West fleet", "East fleet"))
  ob <- rbind(data.frame(Y = catch$Year, Fleet = "West fleet", Value = catch$West),
              data.frame(Y = catch$Year, Fleet = "East fleet", Value = catch$East))
  ob$Fleet <- factor(ob$Fleet, levels(d$Fleet))
  ob <- merge(ob, data.frame(Stage = levels(d$Stage)))
  ggplot(d, aes(Y, Value / 1000)) +
    geom_line(aes(group = Sim), colour = "#b0b0ab", linewidth = 0.3) +
    geom_point(data = ob, colour = INK, size = 1.4) +
    facet_grid(Fleet ~ Stage, scales = "free_y") +
    labs(x = NULL, y = "Catch (kt)") + theme_nfs()
}

# 11. Effort solver: convergence trace ----------------------------------------------
plot_effort_trace <- function(demo, tol = 0.01) {
  tr <- demo$trace
  d <- rbind(data.frame(Iteration = tr$Iteration, Measure = "Worst simulation", v = tr$MaxAbsLogDev),
             data.frame(Iteration = tr$Iteration, Measure = "Median simulation", v = tr$MedianSimMaxDev))
  ggplot(d, aes(Iteration, v, linetype = Measure)) +
    geom_hline(yintercept = tol, linetype = 3, colour = INK, linewidth = 0.4) +
    geom_line(colour = INK, linewidth = 0.6) + geom_point(colour = INK, size = 2) +
    scale_y_log10() + scale_x_continuous(breaks = tr$Iteration) +
    labs(x = "Iteration", y = "Max |log(predicted / observed catch)|") + theme_nfs()
}

# Catch-fit summary of an OM: relative error by fleet over years and simulations
catch_fit_table <- function(h, keep, catch) {
  al <- subset(annual_landings(h), Sim %in% keep & Y >= min(catch$Year))
  al$Obs <- ifelse(al$Fleet == "West", catch$West[match(al$Y, catch$Year)],
                   catch$East[match(al$Y, catch$Year)])
  al$RelErr <- 100 * (al$Value / al$Obs - 1)
  do.call(rbind, lapply(c("West", "East"), function(f) {
    a <- al[al$Fleet == f, ]
    data.frame(Fleet = f, `Observed catch 2005-2025 (kt)` = sum(catch[[f]]) / 1000,
               `Median relative error (%)` = median(a$RelErr),
               `Max absolute relative error (%)` = max(abs(a$RelErr)),
               check.names = FALSE)
  }))
}

# JABBA emulation ---------------------------------------------------------------------
COH_KEY <- c(WS = "Winter-spring cohort", AUT = "Autumn cohort")

# JABBA process-error deviations (median and 95% CI) by cohort
plot_jabba_procdev <- function(fits) {
  d <- do.call(rbind, lapply(names(fits), function(s) {
    ts <- fits[[s]]$timeseries[, , "procB"]
    data.frame(Stock = COH_KEY[[s]], Year = as.numeric(rownames(ts)),
               med = ts[, "mu"], lo = ts[, "lci"], hi = ts[, "uci"])
  }))
  d$Stock <- factor(d$Stock, names(COH_COL))
  ggplot(d, aes(Year)) +
    geom_hline(yintercept = 0, colour = "#b0b0ab", linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Stock), alpha = 0.2) +
    geom_line(aes(y = med, colour = Stock), linewidth = 0.8) +
    scale_colour_manual(values = COH_COL) + scale_fill_manual(values = COH_COL) +
    facet_wrap(~Stock) + labs(x = NULL, y = "JABBA process deviation (log B)") +
    theme_nfs() + theme(legend.position = "none")
}

# OM SSB/SSB0 against the JABBA B/K draws it was solved to (median and 90% interval)
plot_emu_traj <- function(x) {
  yrs <- as.numeric(colnames(x$draws$WS$P))
  q <- function(m) t(apply(m, 2, quantile, c(0.05, 0.5, 0.95)))
  d <- do.call(rbind, lapply(1:2, function(i) {
    s <- names(x$draws)[i]
    rbind(data.frame(Stock = COH_KEY[[s]], Source = "JABBA B/K (posterior draws)", Year = yrs,
                     q(x$draws[[s]]$P)),
          data.frame(Stock = COH_KEY[[s]], Source = "OM SSB/SSB0", Year = yrs,
                     q(x$rel[, i, ])))
  }))
  names(d)[4:6] <- c("lo", "med", "hi")
  d$Stock <- factor(d$Stock, names(COH_COL))
  jb <- subset(d, Source == "JABBA B/K (posterior draws)")
  om <- subset(d, Source == "OM SSB/SSB0")
  # the OM matches each draw to within 1%, so it is drawn as points and bars on
  # top of the JABBA line and band (a line would sit exactly under the JABBA median)
  ggplot(jb, aes(Year, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Stock), alpha = 0.2) +
    geom_line(aes(colour = Stock), linewidth = 0.8) +
    geom_errorbar(data = om, aes(ymin = lo, ymax = hi), colour = INK, width = 0.3,
                  linewidth = 0.3) +
    geom_point(data = om, aes(shape = Source), colour = INK, fill = "white", size = 1.8,
               stroke = 0.6) +
    scale_shape_manual(values = 21, labels = "OM SSB/SSB0 (median and 90% interval)") +
    scale_colour_manual(values = COH_COL, guide = "none") +
    scale_fill_manual(values = COH_COL, guide = "none") +
    facet_wrap(~Stock, scales = "free_y") +
    labs(x = NULL, y = "Biomass relative to unfished") + theme_nfs()
}

# log recruitment multiplier of the cohort that spawns in each target year
# (winter-spring: recruits of calendar year y - 1; autumn: recruits of Jul y-1 to Jun y)
emu_cohort_logdev <- function(x, s) {
  d   <- x$recdev[[s]]
  tt  <- as.numeric(colnames(d)); tyr <- floor(tt + 1e-6); tmon <- rep(1:12, length.out = length(tt))
  py  <- as.numeric(colnames(x$draws[[s]]$P))
  out <- sapply(py, function(y) {
    w <- if (s == "WS") which(tyr == y - 1) else which((tyr == y - 1 & tmon >= 7) | (tyr == y & tmon <= 6))
    log(apply(d[, w, drop = FALSE], 1, mean))
  })
  colnames(out) <- py
  out
}
plot_emu_recdev <- function(x, fits) {
  d <- do.call(rbind, lapply(c("WS", "AUT"), function(s) {
    ld <- emu_cohort_logdev(x, s)
    om <- data.frame(Stock = COH_KEY[[s]], Source = "OM recruitment multiplier (log)",
                     Year = as.numeric(colnames(ld)), med = apply(ld, 2, median),
                     lo = apply(ld, 2, quantile, 0.05), hi = apply(ld, 2, quantile, 0.95))
    ts <- fits[[s]]$timeseries[, , "procB"]
    jb <- data.frame(Stock = COH_KEY[[s]], Source = "JABBA process deviation",
                     Year = as.numeric(rownames(ts)), med = ts[, "mu"], lo = NA, hi = NA)
    rbind(om, jb)
  }))
  d$Stock <- factor(d$Stock, names(COH_COL))
  ggplot(d, aes(Year, med, linetype = Source)) +
    geom_hline(yintercept = 0, colour = "#b0b0ab", linewidth = 0.3) +
    geom_ribbon(data = subset(d, !is.na(lo)), aes(ymin = lo, ymax = hi, fill = Stock),
                alpha = 0.2, linetype = 0) +
    geom_line(aes(colour = Stock), linewidth = 0.8) +
    scale_colour_manual(values = COH_COL, guide = "none") +
    scale_fill_manual(values = COH_COL, guide = "none") +
    scale_linetype_manual(values = c(2, 1)) + facet_wrap(~Stock) +
    labs(x = "Spawning year", y = "Log deviation") + theme_nfs()
}

plot_emu_trace <- function(x, tol = 0.01) {
  tr <- x$emu_trace
  d <- rbind(data.frame(Iteration = tr$Iteration, What = "Catch", v = tr$Catch),
             data.frame(Iteration = tr$Iteration, What = "Winter-spring B/K", v = tr$WS),
             data.frame(Iteration = tr$Iteration, What = "Autumn B/K", v = tr$AUT))
  ggplot(d, aes(Iteration, v, linetype = What)) +
    geom_hline(yintercept = tol, linetype = 3, colour = INK, linewidth = 0.4) +
    geom_line(colour = INK, linewidth = 0.6) + geom_point(colour = INK, size = 1.6) +
    scale_y_log10() + labs(x = "Iteration", y = "Worst |log(OM / target)|") + theme_nfs()
}

# JABBA fits to the CPUE indices: observed (points) and predicted (median, 95% CI)
plot_jabba_fits <- function(fits) {
  d <- do.call(rbind, lapply(names(fits), function(s) {
    g <- fits[[s]]$diags
    data.frame(Stock = COH_KEY[[s]], Year = g$year, obs = g$obs,
               hat = g$hat, lo = g$hat.lci, hi = g$hat.uci)
  }))
  d$Stock <- factor(d$Stock, names(COH_COL))
  ggplot(d, aes(Year)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Stock), alpha = 0.2) +
    geom_line(aes(y = hat, colour = Stock), linewidth = 0.8) +
    geom_point(aes(y = obs), colour = INK, size = 1.8) +
    scale_colour_manual(values = COH_COL, guide = "none") +
    scale_fill_manual(values = COH_COL, guide = "none") +
    facet_wrap(~Stock, scales = "free") +
    labs(x = NULL, y = "GAM-standardized CPUE") + theme_nfs()
}

# summary of the index fits (log scale)
jabba_fit_table <- function(fits) {
  do.call(rbind, lapply(names(fits), function(s) {
    g <- fits[[s]]$diags
    data.frame(Cohort = c(WS = "Winter-spring", AUT = "Autumn")[[s]],
               Years = paste(range(g$year), collapse = "-"),
               `RMSE (log)` = sqrt(mean(g$residual^2)),
               `Correlation, log observed vs predicted` = cor(log(g$obs), log(g$hat)),
               `Observation SD (posterior median)` = sqrt(median(fits[[s]]$posteriors$tau2)),
               `Process SD (posterior median)` = sqrt(median(fits[[s]]$posteriors$sigma2)),
               check.names = FALSE)
  }))
}
