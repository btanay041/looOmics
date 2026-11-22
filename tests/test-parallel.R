# Multisession workers are fresh R processes and can only load looOmics
# when it is installed in a library. devtools::test() loads the package
# from source via load_all(), so workers started there cannot attach it;
# skip in that setting. The comparison is exercised against the installed
# package under R CMD check / BiocCheck.
loo_installed_in_library <- function() {
  any(vapply(.libPaths(), function(lib)
    dir.exists(file.path(lib, "looOmics")), logical(1)))
}

test_that("sequential and multisession plans give identical results", {
  skip_if(!loo_installed_in_library(),
          "looOmics is not installed in a library (running via load_all())")
  d <- make_sim_data()
  res_seq <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                             parallel = FALSE, seed = 99)
  old <- future::plan()
  on.exit(future::plan(future::sequential), add = TRUE)
  future::plan(future::multisession, workers = 2)
  res_par <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                             parallel = TRUE, seed = 99)
  future::plan(future::sequential)
  
  expect_identical(res_seq@loo_betas, res_par@loo_betas)
  expect_identical(res_seq@influence, res_par@influence)
})

test_that("fixed seed gives reproducible results", {
  d <- make_sim_data()
  r1 <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                        parallel = FALSE, seed = 7)
  r2 <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                        parallel = FALSE, seed = 7)
  expect_identical(r1@shifts, r2@shifts)
})

test_that("rank-deficient LOO refit is recorded, not fatal", {
  # Full design is fine, but dropping S07 makes batch perfectly confounded
  # with condition -> that single LOO refit is rank deficient and must be
  # recorded in @failures instead of aborting the run.
  d <- make_sim_data()
  d$col_data$batch <- factor(c("kit1", "kit1", "kit1", "kit1", "kit1", "kit1",
                               "kit1", "kit2", "kit2", "kit2", "kit2", "kit2"))
  # all condition-A samples and S07 (B) are kit1; other B samples are kit2
  res <- loo_sensitivity(d$assay, ~ condition + batch, col_data = d$col_data,
                         parallel = FALSE)
  expect_s4_class(res, "LOOResult")
  expect_equal(nrow(res@failures), 1)
  expect_equal(res@failures$sample, "S07")
  expect_true(grepl("rank", res@failures$error, ignore.case = TRUE))
})

test_that("fail-soft: worker errors land in failures slot", {
  d <- make_sim_data()
  # engine that throws for left-out S05 only: force via a col covariate
  d$col_data$trigger <- 0
  d$col_data$trigger[5] <- 1
  bad_engine <- function(input, design, contrast = NULL, block = NULL, ...) {
    # throws only for the LOO fit that dropped the trigger sample S05
    if (!any(input@col_data$trigger == 1))
      stop("boom")
    list(betas = matrix(0, nrow(input@assay), ncol(design)),
         se = matrix(1, nrow(input@assay), ncol(design)),
         p = matrix(0.5, nrow(input@assay), ncol(design)))
  }
  register_engine("boomer", bad_engine)
  res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                         engine = "boomer", parallel = FALSE)
  expect_equal(nrow(res@failures), 1)
  expect_equal(res@failures$sample, "S05")
  expect_true(res@influence$error[res@influence$sample == "S05"])
})