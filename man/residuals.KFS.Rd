\name{residuals.KFS}
\alias{residuals.KFS}
\title{Extract Residuals of KFS output}
\usage{
  \method{residuals}{KFS} (object,
    type = c("recursive", "response", "pearson", "state", "deviance"),
    standardize = TRUE, ...)
}
\arguments{
  \item{object}{KFS object}

  \item{type}{character string defining the type of
  residuals. See details.}

  \item{standardize}{logical. Should the residuals be
  standardized? Not applicable for type \code{'response'}.}

  \item{...}{Ignored.}
}
\description{
  Extract Residuals of KFS output
}
\details{
  For object of class KFS, several types of (possibly
  standardized) residuals can be computed: \itemize{

  \item 'recursive': One-step ahead prediction residuals.
  Standardized version is defined as
  \deqn{v_{t,i})/\sqrt{F_{i,t}},} with residuals being
  undefined in diffuse phase. Only supported for fully
  Gaussian models.

  \item 'response': Data minus fitted values,
  \eqn{y-E(y)}{y-E(y)}.

  \item 'pearson': Without standardization,
  \deqn{(y_{t,i}-\theta_{t,i})/\sqrt{V(\mu)_{t,i}}, \quad
  i=1,\ldots,p,t=1,\ldots,n,}{(y[t,i]-\theta[t,i])V(\mu)[t,i]^(-0.5),
  i=1,\ldots,p, t=1,\ldots,n,} where
  \eqn{V(\mu_{t,i})}{V(\mu[t,i])} is the variance function
  of the model. With standardization, above residuals are
  divided by \eqn{\sqrt{1-h_{t,i}}}{(1-h[t,i])^(0.5)},
  where
  \eqn{h_{t,i}=(\frac{\textrm{d}\hat\mu}{\textrm{d}\theta})^2*V(\hat\mu)*V_\theta}{h[t,i]=(d\hat\mu/d\theta)^2*V(\hat\mu)*V[\theta]}.
  For gaussian models, these coincide with the smoothed
  \eqn{\epsilon} disturbance residuals.

  \item 'state': Residuals based on the smoothed
  disturbance terms \eqn{\eta} are defined as
  \deqn{L^{-1}_t \hat \eta_t, \quad
  t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n,} where
  \eqn{L_t}{L[t]} is the lower triangular matrix from
  Cholesky decomposition of \eqn{V_{\eta,t}}{V[\eta,t]}.

  \item 'deviance': Deviance residuals. }
}

