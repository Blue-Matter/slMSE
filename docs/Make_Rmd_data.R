# R code extracted from index.Rmd
# Simulation Methods for the Assessment and MSE for Short-Lived Species
# Tom Carruthers <tom@bluematterscience.com>

runRmd <- function(rmdfile, ignorechunk = list("Installation"), envir = globalenv(),
                   run_eval_false = TRUE, use_rmd_dir = TRUE, echo = FALSE) {

  if (!file.exists(rmdfile)) stop("File not found: ", rmdfile)
  ignorechunk <- unlist(ignorechunk)
  lines  <- readLines(rmdfile, warn = FALSE)
  starts <- grep("^\\s*```\\{r", lines)
  ends   <- grep("^\\s*```\\s*$", lines)

  code <- character(0); skipped <- character(0)
  for (i in seq_along(starts)) {
    s <- starts[i]; e <- ends[ends > s][1]
    header <- sub("^\\s*```\\{r\\s*", "", sub("\\}\\s*$", "", lines[s]))
    opts   <- trimws(strsplit(header, ",")[[1]])
    # label: first unnamed option, or label = "..."
    lab <- if (length(opts) && opts[1] != "" && !grepl("=", opts[1])) opts[1] else {
      l <- grep("^label\\s*=", opts, value = TRUE)
      if (length(l)) gsub("^label\\s*=\\s*|[\"']", "", l[1]) else paste0("unnamed-chunk-", i)
    }
    evalF <- any(grepl("^eval\\s*=\\s*(F|FALSE)$", opts))

    if (lab %in% ignorechunk || as.character(i) %in% ignorechunk ||
        (evalF && !run_eval_false)) {
      skipped <- c(skipped, lab); next
    }
    if (e > s + 1) code <- c(code, paste("# ---- chunk:", lab), lines[(s + 1):(e - 1)])
  }

  bad <- setdiff(ignorechunk, c(skipped, as.character(seq_along(starts))))
  if (length(bad)) warning("ignorechunk labels not found: ", paste(bad, collapse = ", "))

  rfile <- tempfile(fileext = ".R")
  writeLines(code, rfile)
  if (use_rmd_dir) {
    oldwd <- setwd(dirname(normalizePath(rmdfile)))
    on.exit(setwd(oldwd), add = TRUE)
  }
  source(rfile, local = envir, echo = echo, max.deparse.length = Inf)
  invisible(list(script = rfile, skipped = skipped))
}

savefiles <- function(folder, ...) {
  objs  <- list(...)
  nams  <- sapply(as.list(substitute(list(...)))[-1], deparse)
  given <- names(objs)
  if (!is.null(given)) nams[given != ""] <- given[given != ""]   # allow name = obj overrides

  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  paths <- file.path(folder, paste0(nams, ".rds"))
  for (i in seq_along(objs)) saveRDS(objs[[i]], paths[i])
  invisible(paths)
}


runRmd("C:/GitHub/slMSE/docs/index.Rmd",ignorechunk = list("Installation"))

savefiles("C:/Github/slMSE/docs/Rmd_data", Sout, SS_spict, Sout_q, Rout, SS_spict_q, SS_RCM_ASPM,
          SS_RCM_SCAL, hist, anMSE, NFS_Spict_y, NFS_spict_q, NFS_ASPM_q)




