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
