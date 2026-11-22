# Shared simulated test data (deterministic, with a planted influential sample)

make_sim_data <- function() {
    set.seed(42)
    n_feat <- 30
    n_samp <- 12
    samples <- sprintf("S%02d", 1:n_samp)
    condition <- factor(rep(c("A", "B"), each = n_samp / 2))
    batch <- factor(rep(c("kit1", "kit2"), times = n_samp / 2))

    mat <- matrix(rnorm(n_feat * n_samp), n_feat, n_samp,
                  dimnames = list(paste0("g", 1:n_feat), samples))
    # True condition effect for the first 5 features
    mat[1:5, condition == "B"] <- mat[1:5, condition == "B"] + 1.5
    # Plant influential sample: S07 (group B) shifts g1-g5 strongly
    mat[1:5, "S07"] <- mat[1:5, "S07"] + 8

    col_data <- data.frame(condition = condition, batch = batch,
                           depth = rnorm(n_samp, 1e6, 1e5),
                           row.names = samples)
    list(assay = mat, col_data = col_data)
}

make_sim_res <- function() {
    d <- make_sim_data()
    loo_sensitivity(d$assay, ~ condition, col_data = d$col_data,
                    parallel = FALSE)
}
