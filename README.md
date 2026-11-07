# Erlang service-time simulation studies

This folder contains the revised R code and its original helper files. The three studies use `model = "erlang_real"`, `B0 = 1000`, Erlang shape `K = 3`, and the other fixed settings from the supplied experiment. Each parameter case runs all three search methods with five starts.

For Codex Cloud, add these files to a GitHub repository, create a Cloud environment for that repository, and use `bash cloud_install.sh` as the install script. Publish the environment after setup succeeds, then start a cloud task in that environment. The setup script installs R and the packages used by these studies.

## Run

From this folder, run one study at a time:

```text
Rscript run_study.R 1 16
Rscript run_study.R 2 16
Rscript run_study.R 3 16
```

The second argument is the number of worker processes; omit it to use the original setting of 100. The script writes `results/study_1_result_all_list.rds`, etc. It saves after every completed parameter case, so rerunning the same command resumes an interrupted study. These are large simulations; no study is started merely by sourcing `erlang_studies.R`.

The required R packages are `foreach`, `doParallel`, `doSNOW`, `RhpcBLASctl`, `Rcpp`, `RcppArmadillo`, `distributions3`, `extraDistr`, `magrittr`, `purrr`, `rlist`, `assertthat`, and `fields`. `Main.R` compiles `workload_discrete.cpp` through `Rcpp::sourceCpp()`.

## Result structure

- Study 1: `result_all_list$late_0_8$cw_0.5$mu_10`, with the same structure under `late_0_15`.
- Study 2: `result_all_list$late_0_15$mu_12$cw_1.0`.
- Study 3: `result_all_list$early_minus8_0$mu_10`; the other interval keys are `late_0_8` and `symmetric_minus8_8`.

Each parameter entry has `settings`, `result_local_search_plus_dpp`, `result_local_search_only`, and `result_punctual_local_search`. Each method result is a list of five starts. The punctual method is a local search optimized for punctual arrivals and evaluated under the selected non-punctual scenario; it is not the Zacharias and Yunes algorithm.

For example, after `result_all_list <- readRDS("results/study_1_result_all_list.rds")`:

```r
result_all_list$late_0_8$cw_0.5$mu_10$
  result_local_search_plus_dpp[[1]]$ci95$V
```

The five `ci95` entries (`N`, `OT`, `WT`, `IT`, `V`) replace the original `raw` vectors. Each is a one-row data frame with `n`, `mean`, `lower`, and `upper`, giving a 95% Student-t interval for the mean across final outer simulation replications. `WT` is total scheduled-patient waiting time, as in the original `raw$WT`. The existing `perf0`, `perfnew`, schedules, loss paths, and timing fields remain.

## Arrival handling

`uniform` draws offsets on `[0, omega]`; `uniform_early` uses `[-omega, 0]`; `uniform_symmetric` uses `[-omega, omega]`. As in the supplied arrival generator, actual arrival times before 0 are set to 0, and those after `T` are excluded. Thus the effective arrival-time distribution changes at the session boundaries.

No result files are included. The studies require an R installation and substantial computation.
