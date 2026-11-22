\name{LOOResult-methods}
\alias{show,LOOResult-method}
\alias{summary,LOOResult-method}
\docType{methods}
\title{Methods for LOOResult objects}
\value{
  `show` returns invisibly; `summary` returns a per-term data.frame;
  `coef` returns the full-sample coefficient data.frame.
}

\description{
  Standard methods for \code{\linkS4class{LOOResult}} objects.
}
\examples{
set.seed(1)
samples <- paste0("S", seq_len(8))
x <- matrix(rnorm(10 * 8), nrow = 10,
            dimnames = list(paste0("gene", 1:10), samples))
cd <- data.frame(condition = factor(rep(c("A", "B"), each = 4)),
                 row.names = samples)
res <- loo_sensitivity(x, ~ condition, col_data = cd, parallel = FALSE)
show(res)
head(summary(res))
head(coef(res))
}

\section{Methods}{
  \describe{
    \item{\code{signature(object = "LOOResult")}}{
      \describe{
        \item{\code{coef}}{returns the full-sample coefficient table.}
        \item{\code{show}}{prints a compact summary (engine, dimensions,
          flags, profiles).}
        \item{\code{summary}}{prints per-term |DFBETAS| shift summaries.}
      }
    }
  }
}
