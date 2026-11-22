#' @include generics.R class-LOOResult.R class-OmicsInput.R
NULL
# Accessors, coercion, printing ---------------------------------------------

#' Standard methods for LOOResult objects
#'
#' The methods provide compact display, per-term shift summaries, and access
#' to the full-sample coefficient table.
#'
#' @name LOOResult-methods
#' @rdname LOOResult-methods
#' @return \code{show()} prints a summary and returns invisibly; \code{summary()}
#'   returns a data.frame of per-term shift summaries; \code{coef()} returns
#'   the full-sample coefficient data.frame.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(8))
#' x <- matrix(rnorm(10 * 8), nrow = 10,
#'             dimnames = list(paste0("gene", 1:10), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' show(res)
#' head(summary(res))
#' head(coef(res))
NULL

#' @export
setMethod("show", "LOOResult", function(object) {
  cat("LOOResult (engine: ", object@engine, ")\n", sep = "")
  cat("  features:", object@diagnostics$n_features,
      "| samples:", object@diagnostics$n_samples,
      "| terms:", paste(object@diagnostics$terms, collapse = ", "), "\n")
  if (nrow(object@failures) > 0)
    cat("  failed LOO refits:", nrow(object@failures), "\n")
  fl <- object@flags
  if (nrow(fl) > 0) {
    cat("  flagged influential samples:", sum(fl$flagged), "/",
        nrow(fl), " (", fl$method[1], ")\n", sep = "")
  } else cat("  flags: not computed (see flag_influential())\n")
  pr <- object@profiles
  if (length(pr) > 0) cat("  profiles: computed (see profiles())\n")
  else cat("  profiles: not computed (see profile_flags())\n")
  invisible(NULL)
})
#' Display an OmicsInput object
#'
#' Prints a compact summary of the assay dimensions, metadata,
#' engine, transformation, and input provenance.
#'
#' @param object An \code{\linkS4class{OmicsInput}} object.
#' @return Prints a summary and returns invisibly.
#' @export
setMethod("show", "OmicsInput", function(object) {
  cat("OmicsInput\n")
  
  cat("  assay:", nrow(object@assay), "features x",
      ncol(object@assay), "samples\n")
  
  cat("  feature metadata:", ncol(object@row_data), "columns\n")
  
  cat("  sample metadata:", ncol(object@col_data), "columns\n")
  
  cat("  engine:", object@engine, "\n")
  
  transform <- if (nzchar(object@transform)) {
    object@transform
  } else {
    "none"
  }
  cat("  transform:", transform, "\n")
  
  cat("  source class:", object@data_class, "\n")
  
  if (nzchar(object@formula)) {
    cat("  formula:", object@formula, "\n")
  }
  
  invisible(NULL)
})
#' @export
setMethod("summary", "LOOResult", function(object, ...) {
  s <- object@shifts[!object@shifts$error, , drop = FALSE]
  agg <- if (nrow(s) == 0L) {
    data.frame(term = character(), abs_dfbetas.median = numeric(),
               abs_dfbetas.max = numeric(), stringsAsFactors = FALSE)
  } else {
    stats::aggregate(cbind(abs_dfbetas = abs(s$dfbetas)) ~ term,
                     data = s, FUN = function(z)
                       c(median = stats::median(z, na.rm = TRUE),
                         max = max(z, na.rm = TRUE)))
  }
  structure(agg, class = c("summary.LOOResult", "data.frame"))
})

#' @export
print.summary.LOOResult <- function(x, ...) {
  if (nrow(x) == 0L) {
    cat("No successful LOO fits.\n")
  } else {
    cat("Per-term shift summary (|DFBETAS|):\n")
    print.data.frame(x, row.names = FALSE)
  }
  invisible(x)
}

#' Extract the assay matrix from an OmicsInput
#'
#' @param x An \code{\linkS4class{OmicsInput}}.
#' @param i Missing. The OmicsInput contains a single assay, so no index
#'   is required.
#' @param withDimnames Logical. Included for compatibility with the
#'   \code{assay} generic. The stored assay's dimnames are preserved.
#' @param ... Unused.
#' @return The numeric assay matrix (features in rows and samples in columns).
#' @examples
#' mat <- matrix(1:12, 3, 4,
#'               dimnames = list(paste0("f", 1:3), paste0("S", 1:4)))
#' cd <- data.frame(group = factor(c("A", "A", "B", "B")),
#'                  row.names = colnames(mat))
#' input <- as_omics_input(mat, ~ group, col_data = cd)
#' assay(input)
#' @export
setMethod("assay", c(x = "OmicsInput", i = "missing"),
          function(x, i, withDimnames = TRUE, ...) x@assay)

#' Extract coefficients from a LOOResult
#'
#' @param object A \code{\linkS4class{LOOResult}}.
#' @param ... Unused.
#' @return A data.frame of full-sample feature-level coefficients, standard
#'   errors, p-values, and model terms.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(8))
#' x <- matrix(rnorm(10 * 8), nrow = 10,
#'             dimnames = list(paste0("gene", 1:10), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' head(coef(res))
#' @export
setMethod("coef", "LOOResult", function(object, ...) object@coefficients)

#' Tidy a LOOResult into a tibble
#'
#' @param x A \code{\linkS4class{LOOResult}}.
#' @param what One of "shifts" (default), "influence", or "flags".
#' @param ... Unused.
#' @return A tibble containing the requested shifts, influence scores, or flags.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' head(tidy(res, what = "influence"))
#' @export
tidy.LOOResult <- function(x, what = c("shifts", "influence", "flags"), ...) {
  what <- match.arg(what)
  out <- switch(what,
                shifts = x@shifts, influence = x@influence, flags = x@flags)
  tibble::as_tibble(out)
}

#' @export
as.data.frame.LOOResult <- function(x, ...) x@shifts

#' Get the flag profile from a LOOResult
#'
#' @param object A \code{\linkS4class{LOOResult}}.
#' @return The \code{FlagProfile} produced by \code{profile_flags()}, or
#'   \code{NULL} if metadata profiling has not been run.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
#'                  batch = factor(rep(c("batch1", "batch2"), 6)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' res <- flag_influential(res, method = "topk", k = 3)
#' res <- profile_flags(res, multivariable = FALSE)
#' profiles(res)
#' @export
profiles <- function(object) object@profiles$profile

#' Per-term influence scores
#'
#' Returns the per-sample influence aggregates computed separately for each
#' model term. The overall \code{influence} slot pools all terms; this
#' accessor shows which model terms drive a sample's score.
#'
#' @param object A \code{\linkS4class{LOOResult}}.
#' @return A data.frame with columns term, sample, max_abs_shift, n_sig_flips.
#' @examples
#' set.seed(1)
#' samples <- paste0("S", seq_len(12))
#' x <- matrix(rnorm(20 * 12), nrow = 20,
#'             dimnames = list(paste0("gene", 1:20), samples))
#' cd <- data.frame(condition = factor(rep(c("A", "B"), each = 6)),
#'                  row.names = samples)
#' res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
#' head(influence_by_term(res))
#' @export
influence_by_term <- function(object) object@diagnostics$influence_by_term