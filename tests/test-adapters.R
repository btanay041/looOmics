test_that("matrix adapter builds a valid OmicsInput", {
  d <- make_sim_data()
  inp <- as_omics_input(d$assay, ~ condition, col_data = d$col_data)
  expect_s4_class(inp, "OmicsInput")
  expect_identical(dim(inp@assay), dim(d$assay))
  expect_identical(colnames(inp@assay), rownames(d$col_data))
  expect_identical(inp@engine, "lm")
})

test_that("adapter rejects mismatched col_data", {
  d <- make_sim_data()
  cd <- d$col_data
  rownames(cd) <- paste0("X", 1:12)
  expect_error(as_omics_input(d$assay, ~ condition, col_data = cd),
               "rownames")
})

test_that("adapter rejects missing values in the assay", {
  d <- make_sim_data()
  a <- d$assay; a[1, 1] <- NA
  expect_error(as_omics_input(a, ~ condition, col_data = d$col_data),
               "finite")
})

test_that("adapter rejects duplicate sample names", {
  d <- make_sim_data()
  inp <- as_omics_input(d$assay, ~ condition, col_data = d$col_data)
  bad <- inp
  colnames(bad@assay) <- rep("S01", 12)
  expect_error(validObject(bad), "unique")
})

test_that("transforms work", {
  d <- make_sim_data()
  inp <- as_omics_input(abs(d$assay), ~ condition, col_data = d$col_data,
                        transform = "log2")
  expect_equal(inp@assay[1, 1], log2(abs(d$assay[1, 1]) + 1))
  expect_error(as_omics_input(d$assay, ~ condition, col_data = d$col_data,
                              transform = "bogus"), "Unknown transform")
})

test_that("CLR centers each sample across features and keeps orientation", {
  a <- matrix(c(10, 100, 1000, 20, 200, 2000), nrow = 3,
              dimnames = list(c("f1", "f2", "f3"), c("s1", "s2")))
  cd <- data.frame(x = 1:2, row.names = c("s1", "s2"))
  inp <- as_omics_input(a, ~ x, col_data = cd, transform = "clr")
  # orientation preserved: features x samples
  expect_identical(dim(inp@assay), dim(a))
  # each column (sample) is centered: column means of CLR values are 0
  expect_equal(colMeans(inp@assay), c(s1 = 0, s2 = 0))
  # reference value computed by hand for s1/f1
  logs <- log(a + 1)
  expect_equal(inp@assay["f1", "s1"], logs["f1", "s1"] - mean(logs[, "s1"]))
})

test_that("pseudocount is applied exactly once", {
  a <- matrix(c(0, 4, 0, 4), nrow = 2, dimnames = list(c("f1", "f2"), c("s1", "s2")))
  cd <- data.frame(x = 1:2, row.names = c("s1", "s2"))
  inp <- as_omics_input(a, ~ x, col_data = cd, transform = "clr", pseudocount = 1)
  # f1 has all-zero row: log(0 + 1) = 0 for both samples
  expect_equal(inp@assay["f1", "s1"], 0 - mean(c(0, log(5))))
})

test_that("data.frame adapter passes through matrix adapter", {
  d <- make_sim_data()
  inp <- as_omics_input(as.data.frame(d$assay), ~ condition,
                        col_data = d$col_data)
  expect_s4_class(inp, "OmicsInput")
})


test_that("model formula is stored as an evaluated formula string", {
  d <- make_sim_data()
  input <- as_omics_input(d$assay, formula = ~ condition,
                          col_data = d$col_data, engine = "lm")
  expect_match(input@formula, "condition")
  expect_false(identical(input@formula, "formula"))
})


test_that("SummarizedExperiment adapter preserves assay and sample metadata", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  d <- make_sim_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = d$assay),
    colData = S4Vectors::DataFrame(d$col_data)
  )
  inp <- as_omics_input(se, ~ condition, engine = "lm",
                        assay_name = "logcounts")
  expect_s4_class(inp, "OmicsInput")
  expect_identical(dim(inp@assay), dim(d$assay))
  expect_identical(rownames(inp@col_data), colnames(d$assay))
  expect_match(inp@formula, "condition")
})


test_that("matrix adapter rejects infinite assay values", {
  d <- make_sim_data()
  a <- d$assay
  a[1, 1] <- Inf
  expect_error(as_omics_input(a, ~ condition, col_data = d$col_data),
               "finite values")
})