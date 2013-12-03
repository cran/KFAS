\name{KFS}
\alias{KFS}
\title{Kalman Filter and Smoother with Exact Diffuse Initialization for Exponential Family State Space Models}
\usage{
  KFS(model, filtering, smoothing, simplify = TRUE,
    transform = c("ldl", "augment"), nsim = 0, theta,
    maxiter = 25)
}
\arguments{
  \item{model}{Object of class \code{SSModel}.}

  \item{filtering}{Types of filtering. Possible choices are
  'state', 'signal', 'invlink', and 'none'. For Gaussian
  models the default is 'state', and for non-Gaussian
  models 'none'.  Multiple values are allowed.}

  \item{smoothing}{Types of smoothing. Possible choices are
  'state', 'signal', 'invlink', 'disturbance' and 'none'.
  Default is 'state'. For non-Gaussian models, option
  'disturbance' is not supported, and for Gaussian models
  option 'invlink' is not supported (as it is identical
  with 'signals'). Multiple values are allowed. Note that
  filtering for non-Gaussian models with importance
  sampling can be very slow with large models. Also, in
  approximating invlink filtering only diagonals of P_mu
  are returned.}

  \item{simplify}{If FALSE, KFS returns some generally not
  so interesting variables from filtering and smoothing.
  Default is TRUE.}

  \item{transform}{how to transform the model in case of
  non-diagonal covariance matrix \eqn{H}. Defaults to
  \code{'svd'}. See function \code{\link{transformSSM}} for
  details.}

  \item{nsim}{Number of independent samples.  Only used for
  non-Gaussian model. Default is 0, which computes the
  approximating Gaussian model by \code{\link{approxSSM}}
  and performs the usual Gaussian smoothing so that the
  smoothed state estimates equals to the conditional mode
  of \eqn{p(\alpha_t|y)}{p(\alpha[t]|y)}.}

  \item{theta}{Initial values for conditional mode theta.
  Only used for non-Gaussian model.}

  \item{maxiter}{Maximum number of iterations used in
  linearisation. Default is 25. Only used for non-Gaussian
  model.}
}
\value{
  For Gaussian model, a list with the following components:
  \item{model}{Original state space model.  }
  \item{KFS_transform}{Type of H after possible
  transformation.  } \item{logLik}{Value of the
  log-likelihood function.  } \item{a}{One step predictions
  of states, \eqn{a_t=E(\alpha_t | y_{t-1}, \ldots ,
  y_{1})}{a[t]=E(\alpha[t] | y[t-1], \ldots , y[1])}.  }
  \item{P}{Covariance matrices (or the non-diffuse parts)
  of predicted states, \eqn{P_t=Cov(\alpha_t | y_{t-1},
  \ldots , y_{1})}{P[t]=Cov(\alpha[t] | y[t-1], \ldots ,
  y[1])}.  } \item{Pinf}{Diffuse part of \eqn{P_t}{P[t]}. }
  \item{v}{Prediction errors \eqn{v_{t,i} = y_{t,i} -
  Z_{i,t}a_{t,i}, i=1,\ldots,p}{v[t,i] = y[t,i] -
  Z[i,t]a[t,i], i=1,\ldots,p}, where
  \eqn{a_{t,i}=E(\alpha_t | y_{t,i-1}, \ldots, y_{t,1},
  \ldots , y_{1,1})}{a[t,i]=E(\alpha[t] | y[t,i-1], \ldots,
  y[t,1], \ldots , y[1,1])}.  } \item{F}{Prediction error
  variances \eqn{Var(v_{t,i})}{Var(v[t,i])}.  }
  \item{Finf}{Diffuse part of \eqn{F_t}{F[t]}.  }
  \item{d}{The last index of diffuse phase, i.e. the
  non-diffuse phase began from time \eqn{d+1}.  }
  \item{j}{The index of last \eqn{y_{i,t}} of diffuse
  phase.  } \item{alphahat}{Smoothed estimates of states,
  \eqn{E(\alpha_t | y_1, \ldots , y_n)}{E(\alpha[t] | y[1],
  \ldots , y[n])}. Only computed if
  \code{smoothing='states'}.  } \item{V}{Covariances
  \eqn{Var(\alpha_t | y_1, \ldots , y_n).}{Var(\alpha[t] |
  y[1], \ldots , y[n]).} Only computed if
  \code{smoothing='states'}.  } \item{thetahat}{Smoothed
  estimates of signals, \eqn{E(Z_t\alpha_t | y_1, \ldots ,
  y_n)}{E(Z[t]\alpha[t] | y[1], \ldots , y[n])}. Only
  computed if \code{smoothing='signals'}.  }
  \item{V_theta}{Covariances \eqn{Var(Z[t]\alpha_t | y_1,
  \ldots , y_n).}{Var(Z[t]\alpha[t] | y[1], \ldots ,
  y[n])}. Only computed if \code{smoothing='signals'}.  }
  \item{etahat}{Smoothed disturbance terms \eqn{E(\eta_t |
  y_1, \ldots , y_n)}{E(\eta[t] | y[1], \ldots , y[n])}.
  Only computed if \code{smoothing='disturbances'}.  }
  \item{V_eta}{Covariances \eqn{Var(\eta_t | y_1, \ldots ,
  y_n)}{Var(\eta[t] | y[1], \ldots , y[n])}. Only computed
  if \code{smoothing='disturbances'}.  }
  \item{epshat}{Smoothed disturbance terms
  \eqn{E(\epsilon_{t,i} | y_1, \ldots ,
  y_n)}{E(\epsilon[t,i] | y[1], \ldots , y[n])}. Only
  computed if \code{smoothing='disturbances'}. Note that
  due to possible diagonalization, these are on transformed
  scale. } \item{V_eps}{Diagonal elements of
  \eqn{Var(\epsilon_{t} | y_1, \ldots ,
  y_n)}{Var(\epsilon[t] | y[1], \ldots , y[n])}. Note that
  due to the diagonalization, off-diagonal elements are
  zero. Only computed if \code{smoothing='disturbances'}.
  } In addition, if argument \code{simplify=FALSE}, list
  contains following components: \item{K}{Covariances
  \eqn{Cov(\alpha_{t,i}, y_{t,i} | y_{t,i-1}, \ldots,
  y_{t,1}, y_{t-1}, \ldots , y_{1}), \quad
  i=1,\ldots,p}{Cov(\alpha[t,i], y[t,i] | y[t,i-1], \ldots,
  y[t,1], y[t-1], \ldots , y[1]), i=1,\ldots,p}.  }
  \item{Kinf}{Diffuse part of \eqn{K_t}{K[t]}.  }
  \item{r}{Weighted sums of innovations \eqn{v_{t+1},
  \ldots , v_{n}}{v[t+1], \ldots , v[n]}.  Notice that in
  literature t in \eqn{r_t}{r[t]} goes from \eqn{0, \ldots,
  n}. Here \eqn{t=1, \ldots, n+1}. Same applies to all r
  and N variables.  } \item{r0, r1}{Diffuse phase
  decomposition of \eqn{r_t}{r[t]}.  } \item{N}{Covariances
  \eqn{Var(r_t)}{Var(r[t])} .  } \item{N0, N1, N2}{Diffuse
  phase decomposition of \eqn{N_t}{N[t]}.  }

  For non-Gaussian model, a list with the following
  components: \item{model}{Original state space model.  }
  \item{alphahat}{Smoothed estimates of states,
  \eqn{E(\alpha_t | y_1, \ldots , y_n)}{E(\alpha[t] | y[1],
  \ldots , y[n])}. Only computed if
  \code{smoothing='state'}.  } \item{V}{Covariances
  \eqn{Var(\alpha_t | y_1, \ldots , y_n).}{Var(\alpha[t] |
  y[1], \ldots , y[n]).} Only computed if
  \code{smoothing='state'}.  } \item{thetahat}{Smoothed
  estimates of signals, \eqn{E(Z_t\alpha_t | y_1, \ldots ,
  y_n)}{E(Z[t]\alpha[t] | y[1], \ldots , y[n])}. Only
  computed if \code{smoothing='signal'}.  }
  \item{V_theta}{Covariances \eqn{Var(Z[t]\alpha_t | y_1,
  \ldots , y_n).}{Var(Z[t]\alpha[t] | y[1], \ldots ,
  y[n]).} Only computed if \code{smoothing='signal'}.  }
  \item{muhat}{Smoothed estimates of \eqn{f(\theta_t) |
  y_1, \ldots , y_n)}{f(\theta[t]) | y[1], \ldots , y[n])},
  where \eqn{f} is the inverse link function. Only computed
  if \code{smoothing='invlink'}.  } \item{V_mu}{Covariances
  \eqn{Cov(f(\theta_t)| y_1, \ldots ,
  y_n).}{Cov(f(\theta[t]) | y[1], \ldots , y[n]).} Only
  computed if \code{smoothing='invlink'}. If \code{nsim=0},
  only diagonal elements (variances) are computed, using
  the delta method.  } \item{iterations}{The number of
  iterations used in linearization.}
}
\description{
  Performs Kalman filtering and smoothing with exact
  diffuse initialization using univariate approach for
  exponential family state space models.
}
\details{
  Notice that in case of multivariate observations,
  \code{v}, \code{F}, \code{Finf}, \code{K} and \code{Kinf}
  are usually not the same as those calculated in usual
  multivariate Kalman filter. As filtering is done one
  observation element at the time, the elements of
  prediction error \eqn{v_t}{v[t]} are uncorrelated, and
  \code{F}, \code{Finf}, \code{K} and \code{Kinf} contain
  only the diagonal elemens of the corresponding covariance
  matrices.

  In rare cases of a very long diffuse initialization phase
  with highly correlated states, cumulative rounding errors
  in computing \code{Finf} and \code{Pinf} can sometimes
  cause the diffuse phase end too early. Changing the
  tolerance parameter \code{tol} of the model (see
  \code{\link{SSModel}}) to smaller (or larger) should
  help.

  In case of non-Gaussian models with \code{nsim=0}, the
  smoothed estimates relate the conditional mode of
  \eqn{p(\alpha|y)}, and are comparable with the results
  from generalized linear models. When using importance
  sampling (\code{nsim>0}), results correspond to the
  conditional mean.
}
\references{
  Koopman, S.J. and Durbin J. (2000).  Fast filtering and
  smoothing for non-stationary time series models, Journal
  of American Statistical Assosiation, 92, 1630-38.  \cr

  Koopman, S.J. and Durbin J. (2001).  Time Series Analysis
  by State Space Methods. Oxford: Oxford University Press.
  \cr

  Koopman, S.J. and Durbin J. (2003).  Filtering and
  smoothing of state vector for diffuse state space models,
  Journal of Time Series Analysis, Vol. 24, No. 1.  \cr
}

