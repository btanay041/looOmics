test_that("LOO betas match a manual lm loop exactly", {
  d <- make_sim_data()
  res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                         parallel = FALSE)
  # Manual oracle for feature g1, term conditionB, left-out S01
  i <- 1; g <- "g1"
  dat <- data.frame(y = d$assay[g, -i],
                    cond = droplevels(d$col_data$condition[-i]))
  manual <- coef(lm(y ~ cond, data = dat))
  got <- res@loo_betas[paste(g, "conditionB", sep = "::"), "S01"]
  expect_equal(unname(manual[["condB"]]), unname(got), tolerance = 1e-10)
  
  # Across all samples for one feature
  manual_all <- vapply(seq_len(12), function(i) {
    dat <- data.frame(y = d$assay[g, -i],
                      cond = droplevels(d$col_data$condition[-i]))
    unname(coef(lm(y ~ cond, data = dat))["condB"])
  }, numeric(1))
  expect_equal(unname(manual_all), unname(res@loo_betas[paste(g, "conditionB", sep = "::"), ]),
               tolerance = 1e-10)
})

test_that("full-sample coefficients match lm directly", {
  d <- make_sim_data()
  res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                         parallel = FALSE)
  dat <- data.frame(y = d$assay["g1", ], cond = d$col_data$condition)
  fit <- lm(y ~ cond, data = dat)
  cf <- coef(summary(fit))
  row <- res@coefficients[res@coefficients$feature == "g1" &
                            res@coefficients$term == "conditionB", ]
  expect_equal(row$beta, unname(cf["condB", "Estimate"]), tolerance = 1e-10)
})

test_that("planted influential sample S07 tops the influence score", {
  res <- make_sim_res()
  top <- res@influence$sample[which.max(res@influence$max_dfbetas)]
  expect_equal(top, "S07")
})

test_that("unknown engine errors informatively", {
  d <- make_sim_data()
  expect_error(loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                               engine = "nope"),
               "Unknown engine")
})

test_that("glm engine matches manual LOO glm fits", {
  d <- make_sim_data()
  res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                         engine = "glm", parallel = FALSE)
  i <- 3; g <- "g2"
  dat <- data.frame(y = d$assay[g, -i],
                    cond = droplevels(d$col_data$condition[-i]))
  manual <- coef(glm(y ~ cond, data = dat))
  got <- res@loo_betas[paste(g, "conditionB", sep = "::"), "S03"]
  expect_equal(unname(manual[["condB"]]), unname(got), tolerance = 1e-8)
})

test_that("custom engine registration works", {
  my_engine <- function(input, design, contrast = NULL, block = NULL, ...) {
    # trivial identity engine: beta = rowMeans difference proxy
    g <- nrow(input@assay)
    p <- ncol(design)
    m <- matrix(rnorm(g * p), g, p, dimnames = list(rownames(input@assay),
                                                    colnames(design)))
    list(betas = m, se = abs(m) + 1, p = matrix(0.1, g, p))
  }
  register_engine("toy", my_engine)
  d <- make_sim_data()
  res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                         engine = "toy", parallel = FALSE)
  expect_equal(res@engine, "toy")
  expect_equal(nrow(res@coefficients), 30 * 2)
})

test_that("glm engine is family-aware: gaussian uses t-based inference", {
  d <- make_sim_data()
  res_glm <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                             engine = "glm", parallel = FALSE)
  res_lm <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                            parallel = FALSE)
  # Gaussian glm p-values must match the lm engine (t-based, estimated dispersion)
  expect_equal(res_glm@coefficients$p, res_lm@coefficients$p, tolerance = 1e-8)
  expect_equal(res_glm@coefficients$se, res_lm@coefficients$se, tolerance = 1e-8)
})

test_that("glm engine: binomial uses fixed dispersion and z-based inference", {
  set.seed(7)
  # The binomial family models binary 0/1 outcomes (or counts with
  # weights) - not fractional values, which would trigger
  # "non-integer #successes" warnings. Simulate true binary responses
  # with a group-dependent success probability.
  grp <- factor(rep(c("A", "B"), each = 4))
  cd <- data.frame(condition = grp, row.names = sprintf("S%02d", 1:8))
  p <- ifelse(grp == "B", 0.85, 0.25)
  a <- matrix(rbinom(5 * 8, size = 1, prob = p), 5, 8,
              dimnames = list(paste0("b", 1:5), rownames(cd)))
  res <- loo_sensitivity(a, ~ condition, col_data = cd, engine = "glm",
                         family = binomial(), parallel = FALSE)
  # compare full-sample fit against stats::glm directly
  dat <- data.frame(y = a["b1", ], cond = cd$condition)
  m <- summary(glm(y ~ cond, data = dat, family = binomial()))
  row <- res@coefficients[res@coefficients$feature == "b1" &
                            res@coefficients$term == "conditionB", ]
  expect_equal(row$beta, unname(m$coefficients["condB", "Estimate"]),
               tolerance = 1e-6)
  expect_equal(row$p, unname(m$coefficients["condB", "Pr(>|z|)"]),
               tolerance = 1e-4)
})