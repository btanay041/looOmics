test_that("OmicsInput validity catches dimension mismatches", {
  d <- make_sim_data()
  inp <- as_omics_input(d$assay, ~ condition, col_data = d$col_data)
  bad <- inp
  bad@col_data <- d$col_data[1:6, ]
  expect_error(validObject(bad), "ncol")
})

test_that("LOOResult validity catches betas/se mismatch", {
  res <- make_sim_res()
  bad <- res
  bad@loo_se <- bad@loo_se[, 1:6, drop = FALSE]
  expect_error(validObject(bad), "identical dimensions")
})

test_that("show and summary methods run", {
  res <- make_sim_res()
  expect_invisible(show(res))
  s <- summary(res)
  expect_s3_class(s, "summary.LOOResult")
  expect_true(all(c("term", "abs_dfbetas") %in% names(s)))
  expect_equal(colnames(s$abs_dfbetas), c("median", "max"))
  expect_invisible(print(s))
})