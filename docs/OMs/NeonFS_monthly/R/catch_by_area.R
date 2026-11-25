# Annual neon flying squid catch west and east of 170E, by member and summed.
# Source: NPFC-2026-AR Annual Summary Footprint - Squids (Rev3), sheet 'Catch'.
#
# Redaction rules ('x' = fewer than 3 vessels):
#  - if one areal component is 'x' and the total and other components are known,
#    the component is inferred as total minus the others;
#  - otherwise an 'x' component contributes 0 (i.e. a redacted total is replaced
#    by the sum of the reported areal components);
#  - blanks and 'n/a' are 0.
# Japan 'NW' (national waters) and Russia NW gears are counted as west of 170E.

NFS_catch_by_area <- function(file, years = 2005:2025) {

  raw <- suppressMessages(readxl::read_excel(file, sheet = "Catch", col_names = FALSE,
                                             col_types = "text"))
  raw <- as.data.frame(raw)
  yr  <- suppressWarnings(as.integer(raw[[1]]))
  raw <- raw[!is.na(yr) & yr %in% years, ]
  yr  <- as.integer(raw[[1]])

  num <- function(v) {                         # NA marks a redaction
    v <- trimws(v)
    out <- suppressWarnings(as.numeric(v))
    out[is.na(v) | v %in% c("", "n/a")] <- 0
    out
  }

  # 1-based column indices: Total, West components, East components
  members <- list(
    China          = list(T = 2,  W = 3,          E = 4),
    Japan          = list(T = 8,  W = c(9, 11),   E = 10),
    Korea          = list(T = 14, W = 15,         E = 16),
    Russia         = list(T = 17, W = c(18, 20),  E = 19),
    Chinese_Taipei = list(T = 30, W = 31,         E = 32),
    Vanuatu        = list(T = 33, W = 34,         E = 35)
  )

  out <- list()
  for (m in names(members)) {
    cl  <- members[[m]]
    Tt  <- num(raw[[cl$T]])
    W   <- sapply(cl$W, function(j) num(raw[[j]]))
    E   <- sapply(cl$E, function(j) num(raw[[j]]))
    W   <- matrix(W, nrow = length(yr)); E <- matrix(E, nrow = length(yr))
    comps <- cbind(W, E)
    for (i in seq_along(yr)) {
      miss <- which(is.na(comps[i, ]))
      if (length(miss) == 1 && !is.na(Tt[i]))
        comps[i, miss] <- max(0, Tt[i] - sum(comps[i, -miss]))
      comps[i, is.na(comps[i, ])] <- 0
    }
    nW <- ncol(W)
    out[[m]] <- data.frame(Year = yr, Member = m,
                           West = rowSums(comps[, 1:nW, drop = FALSE]),
                           East = rowSums(comps[, -(1:nW), drop = FALSE]))
  }
  bymember <- do.call(rbind, out)
  total <- aggregate(cbind(West, East) ~ Year, data = bymember, sum)
  total <- total[order(total$Year), ]
  list(by_member = bymember, total = total)
}
