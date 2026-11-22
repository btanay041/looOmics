test_that("coef, tidy, as.data.frame accessors", {
    res <- make_sim_res()
    cf <- coef(res)
    expect_equal(nrow(cf), 30 * 2)
    expect_equal(colnames(cf), c("feature", "term", "beta", "se", "p"))

    tb <- tidy(res, what = "shifts")
    expect_s3_class(tb, "tbl_df")
    expect_equal(nrow(tb), nrow(res@shifts))

    df <- as.data.frame(res)
    expect_equal(nrow(df), nrow(res@shifts))
})

test_that("autoplot returns ggplot for each type", {
    res <- flag_influential(make_sim_res())
    p1 <- autoplot(res, type = "shift")
    p2 <- autoplot(res, type = "scree")
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_error(autoplot(res, type = "bogus"), "'arg' should be")
    expect_s3_class(plot(res), "ggplot")
})
