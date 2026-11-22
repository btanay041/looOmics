test_that("DESeq2 requires raw integer counts and stores the model formula", {
    skip_if_not_installed("DESeq2")
    set.seed(12)
    n <- 12L
    cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6),
                                        levels = c("A", "B")))
    rownames(cd) <- paste0("S", seq_len(n))
    counts <- matrix(rnbinom(80L * n, mu = 80, size = 10), 80L, n,
                     dimnames = list(paste0("g", seq_len(80L)), rownames(cd)))
    input <- as_omics_input(counts, formula = ~ condition, col_data = cd,
                            engine = "deseq2", transform = "none")
    expect_match(input@formula, "condition")
    expect_error(
        loo_sensitivity(counts + 0.5, ~ condition, engine = "deseq2",
                        col_data = cd, parallel = FALSE),
        "non-negative integer raw counts")
    res <- loo_sensitivity(input, ~ condition, parallel = FALSE)
    expect_s4_class(res, "LOOResult")
    expect_true(all(c("conditionB") %in% res@diagnostics$terms))
})

test_that("DESeq2 coefficient names map to non-reference factor levels", {
    cd <- data.frame(group = factor(c("A", "B", "C"),
                                    levels = c("A", "B", "C")))
    expect_equal(looOmics:::.deseq_result_name(
        "groupB", cd, "group_B_vs_A"), "group_B_vs_A")
    expect_equal(looOmics:::.deseq_contrast("groupC", cd),
                 c("group", "C", "A"))
})

test_that("limma engine returns coefficient matrices when available", {
    skip_if_not_installed("limma")
    set.seed(13)
    cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)))
    rownames(cd) <- paste0("S", seq_len(nrow(cd)))
    x <- matrix(rnorm(60L * 12L), 60L, 12L,
                dimnames = list(paste0("g", seq_len(60L)), rownames(cd)))
    x[1:10, cd$condition == "B"] <-
        x[1:10, cd$condition == "B"] + 1
    res <- loo_sensitivity(x, ~ condition, col_data = cd,
                           engine = "limma", parallel = FALSE)
    expect_s4_class(res, "LOOResult")
    expect_true(all(c("feature", "term", "beta", "p") %in%
                    names(res@coefficients)))
})
