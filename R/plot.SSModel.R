#' Diagnostic Plots of State Space Models
#'
#' Diagnostic plots based on standardized residuals for objects of class \code{SSModel}.
#'
#' @export
#' @importFrom graphics plot par hist lines
#' @importFrom stats na.pass density acf
#' @param x Object of class \code{SSModel}.
#' @param nsim The number of independent samples used in importance sampling.
#' Only used for non-Gaussian model. Default is 0, which computes the
#' approximating Gaussian model by \code{\link{approxSSM}} and performs the
#' usual Gaussian filtering/smoothing so that the smoothed state estimates
#' equals to the conditional mode of \eqn{p(\alpha_t|y)}{p(\alpha[t]|y)}.
#' In case of \code{nsim = 0}, the mean estimates and their variances are computed using
#' the Delta method (ignoring the covariance terms).
#' @param zerotol Tolerance parameter for positivity checking in standardization. Default is zero. 
#' The values of D <= zerotol * max(D, 0) are deemed to zero.
#' @param expected Logical value defining the approximation of H_t in case of Gamma 
#' and negative binomial distribution. Default is \code{FALSE} which matches the 
#' algorithm of Durbin & Koopman (1997), whereas \code{TRUE} uses the expected value
#' of observations in the equations, leading to results which match with \code{glm} (where applicable).
#' The latter case was the default behaviour of KFAS before version 1.3.8.
#' Essentially this is the difference between observed and expected information.
#' @param ... Ignored.
#' @examples
#' modelNile <- SSModel(Nile ~ SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
#' modelNile <- fitSSM(inits = c(log(var(Nile)),log(var(Nile))), model = modelNile,
#'  method = "BFGS")$model
#'
#' if (interactive()) {
#'   plot(modelNile)
#' }

plot.SSModel <- function(x, nsim = 0, zerotol = 0, expected = FALSE, ...) {

  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  acf2 <- function(x, acf_args) {
    tsp(x) <- NULL
    tmp <- acf(x, plot = FALSE, na.action = na.pass)
    tmp$n.used <- sum(!is.na(x))
    if (is.null(ncol(x)) || ncol(x) == 1) {
      plot(tmp, main = "ACF of recursive residuals")
    } else {
      plot(tmp)
    }

  }

  par(ask = TRUE)
  if (all(x$distribution == "gaussian")) {
    out <- KFS(x, smoothing = c("mean", "disturbance"), filtering = c("state","mean"))
    res <- rstandard(out, zerotol = zerotol)
    if (attr(x, "p") == 1) {
      plot(cbind(
          observations = x$y,
          recursive = res,
          irregular = rstandard(out, "pearson", zerotol = zerotol)),
          main = "Recursive (one-step-ahead) and irregular (smoothed) residuals")
      acf2(res)
      hist(res, main = "Histogram of recursive residuals", freq = FALSE)
      lines(density(as.numeric(res), na.rm = TRUE))
    } else {

      plot(x$y, main = "Observations")
      plot(res, main = "Recursive (one-step-ahead) residuals")
      plot(rstandard(out, "pearson", zerotol = zerotol), main = "Irregular (smoothed) residuals")
      acf2(res)
      for (i in 1:attr(x, "p")) {
        hist(res[, i], main = paste("Histogram for recursive residuals of ", colnames(res)[i]),
          freq = FALSE)
        lines(density(res[, i], na.rm = TRUE))
      }

    }

  } else {
    out <- KFS(x, smoothing = "mean", filtering = "mean", nsim = nsim,
      expected = expected)
    res <- rstandard(out, zerotol = zerotol)

    if (attr(x, "p") == 1) {
      plot(cbind(
          observations = x$y,
          recursive = res,
          irregular = rstandard(out, "pearson", zerotol = zerotol)),
          main = "recursive (one-step-ahead) and irregular (smoothed) residuals")
     acf2(res)
     hist(res, main = "Histogram of recursive residuals", freq = FALSE)
     lines(density(as.numeric(res), na.rm = TRUE))
    } else {

      plot(x$y, main = "Observations")
      plot(res, main = "Recursive (one-step-ahead) residuals")
      plot(rstandard(out, "pearson", zerotol = zerotol), main = "Irregular (smoothed) residuals")
      acf2(res)
      for (i in 1:attr(x, "p")) {
        hist(res[, i], main = paste("Histogram for recursive residuals of ", colnames(res)[i]),
          freq = FALSE)
        lines(density(res[, i], na.rm = TRUE))
      }
    }
  }


}
