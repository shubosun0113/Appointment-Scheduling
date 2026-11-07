# Run from a terminal: Rscript run_study.R STUDY_NUMBER [CORES]
# STUDY_NUMBER is 1, 2, or 3. A partially completed study resumes from its RDS file.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L || !args[[1L]] %in% c("1", "2", "3")) {
    stop("Usage: Rscript run_study.R STUDY_NUMBER [CORES]")
  }
  study_number <- args[[1L]]
  cores <- if (length(args) >= 2L) as.integer(args[[2L]]) else 100L
  if (is.na(cores) || cores < 1L) stop("CORES must be a positive integer")

  script_arg <- commandArgs(trailingOnly = FALSE)
  script_arg <- script_arg[grepl("^--file=", script_arg)]
  if (length(script_arg) == 1L) {
    setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
  }
  source("erlang_studies.R")
  dir.create("results", showWarnings = FALSE)
  result_file <- file.path("results", paste0("study_", study_number,
                                            "_result_all_list.rds"))

  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doSNOW::registerDoSNOW(cl)

  runner <- switch(study_number,
                   `1` = run_erlang_study_1,
                   `2` = run_erlang_study_2,
                   `3` = run_erlang_study_3)
  result_all_list <- runner(cl = cl, checkpoint_file = result_file)
  saveRDS(result_all_list, result_file)
  message("Complete: ", normalizePath(result_file))
  invisible(result_all_list)
}

if (sys.nframe() == 0L) main()
