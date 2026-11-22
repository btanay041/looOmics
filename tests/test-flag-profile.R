test_that("threshold method flags samples with max DFBETAS > 1", {
    res <- flag_influential(make_sim_res(), method = "threshold")
    expect_equal(res@flags$flagged, res@flags$score > 1)
})

test_that("mad method flags the planted sample", {
    res <- flag_influential(make_sim_res(), method = "mad", threshold = 4)
    expect_true("S07" %in% res@flags$sample[res@flags$flagged])
})

test_that("topk returns exactly k rows sorted by score", {
    res <- flag_influential(make_sim_res(), method = "topk", k = 3)
    expect_equal(sum(res@flags$flagged), 3)
    top3 <- head(res@flags$sample[order(-res@flags$score)], 3)
    expect_setequal(res@flags$sample[res@flags$flagged], top3)
})

test_that("invalid inputs raise informative errors", {
    expect_error(flag_influential(make_sim_res(), method = "bogus"),
                 "'arg' should be")
    expect_error(flag_influential(list()), "LOOResult")
})

test_that("profile detects the planted batch association", {
    d <- make_sim_data()
    # Plant a metadata association: flagged sample S07 gets kit "bad"
    d$col_data$extraction <- factor(rep(c("good", "good", "bad", "good"), 3))
    d$col_data$extraction[7] <- "bad"
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags(model_covariates = "condition")

    pr <- profiles(res)
    expect_s3_class(pr, "FlagProfile")
    uni <- pr$univariate
    expect_true("extraction" %in% uni$covariate)
    # Model covariates remain available for biological interpretation and
    # are explicitly labeled in the metadata association results.
    expect_true(any(uni$in_model[uni$covariate == "condition"]))
    expect_output(print(pr), "FlagProfile")
})

test_that("profile requires flags first", {
    res <- make_sim_res()
    expect_error(profile_flags(res), "flags")
})

test_that("continuous covariates use wilcox and match wilcox.test", {
    d <- make_sim_data()
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags(covariates = "depth")
    uni <- profiles(res)$univariate
    expect_equal(uni$test, "wilcox")
    fl <- res@flags$flagged; names(fl) <- res@flags$sample
    expect_equal(uni$p_value,
                 suppressWarnings(wilcox.test(d$col_data$depth ~ fl))$p.value)
})

test_that("rare levels are collapsed and single-level covariates dropped", {
    d <- make_sim_data()
    d$col_data$sparse <- factor(rep("lvl", 12))
    d$col_data$two_lvl <- factor(c("rare", "rare", rep(c("a", "b"), 5)))
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags()
    uni <- profiles(res)$univariate
    expect_false("sparse" %in% uni$covariate)
    expect_true("two_lvl" %in% uni$covariate)
})

test_that("n_features_shifted counts unique features, not feature-term pairs", {
    res <- make_sim_res()
    s <- res@shifts
    # A feature with 2 shifted coefficients must count once
    sub <- s[s$sample == "S07", ]
    expected_unique <- length(unique(sub$feature[abs(sub$dfbetas) > 1]))
    expect_equal(res@influence$n_features_shifted[res@influence$sample == "S07"],
                 expected_unique)
})

test_that("per-term influence scores are reported separately", {
    res <- make_sim_res()
    bt <- influence_by_term(res)
    expect_true(all(c("term", "sample", "max_abs_shift", "n_sig_flips") %in%
                    colnames(bt)))
    expect_setequal(unique(bt$term), c("(Intercept)", "conditionB"))
    # pooled max over terms equals the overall max_dfbetas
    pooled <- tapply(bt$max_abs_shift, bt$sample, max, na.rm = TRUE)
    overall <- setNames(res@influence$max_dfbetas, res@influence$sample)
    expect_equal(as.numeric(pooled[names(overall)]), as.numeric(overall))
})

test_that("profile handles NA covariates and reports dropped covariates", {
    d <- make_sim_data()
    d$col_data$depth[c(1, 3, 5)] <- NA          # partial NA: still testable
    d$col_data$allna <- NA                       # fully NA: dropped
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags(covariates = c("depth", "allna"))
    pr <- profiles(res)
    expect_true("depth" %in% pr$univariate$covariate)
    expect_equal(pr$univariate$n_used[pr$univariate$covariate == "depth"], 9)
    expect_true("allna" %in% pr$dropped$covariate)
})

test_that("profile works when every covariate is dropped", {
    d <- make_sim_data()
    d$col_data$only_level <- "x"
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags(covariates = "only_level")
    pr <- profiles(res)
    expect_equal(nrow(pr$univariate), 0)
    expect_true(nrow(pr$dropped) >= 1)
    expect_true(!is.null(pr$mv_note))
    expect_output(print(pr), "not tested")
})

test_that("Fisher is selected using expected cell counts", {
    d <- make_sim_data()
    # 3 "y" vs 9 "n" (survives min_cell_size = 3): expected cell count for
    # the "y" row is 3 * 6 / 12 = 1.5 < 5, so Fisher must be selected
    d$col_data$rare <- factor(c("y", "y", "y", rep("n", 9)))
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags(covariates = "rare")
    uni <- profiles(res)$univariate
    expect_equal(uni$test, "fisher_simulated")
})


test_that("model covariates are profiled and labeled by default", {
    d <- make_sim_data()
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold") |>
        profile_flags(multivariable = FALSE)
    uni <- profiles(res)$univariate
    expect_true("condition" %in% uni$covariate)
    expect_true(all(uni$in_model[uni$covariate == "condition"]))
})

test_that("simulated Fisher tests are reproducible and preserve caller RNG", {
    d <- make_sim_data()
    d$col_data$rare <- factor(c("a", "a", "a", "b", "b", "b",
                                "c", "c", "c", "d", "d", "d"))
    res <- loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                           parallel = FALSE) |>
        flag_influential(method = "threshold")
    set.seed(999)
    before <- .Random.seed
    p1 <- profiles(profile_flags(res, covariates = "rare",
                                 multivariable = FALSE, seed = 123))$univariate$p_value
    expect_identical(.Random.seed, before)
    p2 <- profiles(profile_flags(res, covariates = "rare",
                                 multivariable = FALSE, seed = 123))$univariate$p_value
    expect_equal(p1, p2)
})
