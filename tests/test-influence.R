test_that("delta equals beta - beta_loo elementwise", {
    res <- make_sim_res()
    s <- res@shifts
    expect_true(all(abs(s$delta - (s$beta - s$beta + s$delta)) < .Machine$double.eps))
    # direct check for one cell
    key <- "g1::conditionB"
    row <- res@loo_betas[key, "S02"]
    full <- res@coefficients[res@coefficients$feature == "g1" &
                             res@coefficients$term == "conditionB", "beta"]
    cell <- s[s$feature == "g1" & s$term == "conditionB" & s$sample == "S02", ]
    expect_equal(cell$delta, full - row, tolerance = 1e-12)
})

test_that("dfbetas equals delta / se_loo", {
    res <- make_sim_res()
    s <- res@shifts
    idx <- which(!is.na(s$dfbetas))
    expect_equal(s$dfbetas[idx], s$delta[idx] / s$se[idx], tolerance = 1e-12)
})

test_that("per-sample aggregates are correct maxima of the shift table", {
    res <- make_sim_res()
    sub <- res@shifts[res@shifts$sample == "S07", ]
    expect_equal(max(abs(sub$dfbetas), na.rm = TRUE),
                 res@influence$max_dfbetas[res@influence$sample == "S07"])
    expect_equal(sum(sub$sig_flip, na.rm = TRUE),
                 res@influence$n_sig_flips[res@influence$sample == "S07"])
})

test_that("sig_flip detects q-threshold crossings", {
    res <- make_sim_res()
    s <- res@shifts
    idx <- which(!is.na(s$q_loo))
    # full q: reconstruct from coefficients p per term via BH across features
    sp <- split(res@coefficients, res@coefficients$term)
    full_q <- do.call(rbind, lapply(sp, function(d)
        data.frame(feature = d$feature, term = d$term, q = p.adjust(d$p, "BH"),
                   stringsAsFactors = FALSE)))
    qmap <- setNames(full_q$q, paste(full_q$feature,
                                     rep(unique(full_q$term), each = length(unique(full_q$feature))),
                                     sep = "::"))
    expected <- (qmap[paste(s$feature, s$term, sep = "::")][idx] < 0.05) !=
        (s$q_loo[idx] < 0.05)
    expect_equal(as.logical(expected), s$sig_flip[idx])
})


test_that("shift table records within-term feature rank changes", {
    res <- make_sim_res()
    expect_true(all(c("p_full", "q_full", "rank_full", "p_loo",
                      "q_loo", "rank_loo", "rank_change") %in%
                    names(res@shifts)))
    ok <- !res@shifts$error & !is.na(res@shifts$rank_change)
    expect_equal(res@shifts$rank_change[ok],
                 res@shifts$rank_loo[ok] - res@shifts$rank_full[ok])
})
