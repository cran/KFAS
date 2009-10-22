\name{expflik}
\alias{expflik}
\title{Approximate log-likelihood and gaussian density of exponential family state-space model.}
\description{
  Function expflik computes approximate log-likelihood and approximate gaussian density of univariate exponential family state-space model, based on Durbin and Koopman (1997, 2001).
}
\usage{
  expflik(yt, Zt, Tt, Rt, Qt, a1, P1, P1inf, dist=c("Poisson", "Binomial", "Negative binomial"), offset=1)
}
\arguments{
  \item{yt}{Matrix, array or vector of observations. Note that yt is univariate.}
  \item{Zt}{System matrix or array of observation equation.}
  \item{Tt}{System matrix or array of transition equation.}
  \item{Rt}{System matrix or array of transition equation.}
  \item{Qt}{Variance matrix or array of disturbance terms eta_t.}
  \item{a1}{Initial state vector.}
  \item{P1}{Variance matrix of a1. In diffuse case P1star, the non-diffuse part of P1.}
  \item{P1inf}{Diffuse part of P1.}
  \item{dist}{Distribution of yt.}
  \item{offset}{Vector of length n. See details.}
}

\details{
  
  Function approximates p(alpha|y) with gaussian g(alpha|y) which has same conditional mode (alpha,alpha_n+1|y) as p(alpha|y), and computes approximate log-likelihood\cr


  logL_0(psi) = logL_g(psi) + logE_g(w_T|y),\cr


  where logL_g(psi) is log-likelihood of approximate gaussian distribution and logE_g(w_T|y) is a Taylor-approximation of E_g[p(y|theta)/g(y|theta) | y].\cr


  For details, see J. Durbin and S.J. Koopman (1997).


  The general state space model for exponential family is given by\cr

  p(y_t|theta_t) = exp[theta'_t * y_t - b_t(theta_t) + f_t(y_t)] (observation equation)\cr


  alpha_t+1 = T_t * alpha_t + R_t * eta_t (transition equation)\cr

  where theta_t = Z_t * alpha_t and eta_t ~ N(0,Q_t) \cr
 
  Approximating gaussian model is given by\cr

  y*_t = Z_t * alpha_t + eps_t (observation equation)\cr

  alpha_t+1 = T_t * alpha_t + R_t * eta_t (transition equation)\cr

  where eps_t ~ N(0,H*_t) and eta_t ~ N(0,Q_t) \cr
  
  If yt is Poisson distributed, parameter of Poisson distribution is offset*lambda and theta = log(lambda).\cr
  If yt is from binomial distribution, offset is a vector specifying number of trials at times 1,...,n, and theta = log[pi_t/(1-pi_t)], where pi_t is the probability of success at time t.\cr
  In case of negative binomial distribution, offset is vector of specified number of successes wanted at times 1,...,n, and theta = log(1-pi_t). \cr

  Note that this function works only for univariate observations.
}

\value{
List with output from Kalman smoother and distribution smoother, when model is approximated with gaussian distribution g(y|theta). Note that Ht is H*_t. List also contains following items:
\item{ytilde}{y* of approximating gaussian model.}
\item{theta}{Z_t * alpha_t of approximating gaussian model.}
\item{lik0}{Value of logL_0.}
\item{offset}{}
\item{dist}{Distribution of yt.}
}

\references{
Durbin J. and Koopman, S.J. (1997). Monte Carlo Maximum Likelihood Estimation for Non-Gaussian State Space Models, Biometrica, Vol. 84, No. 3.\cr
Koopman, S.J. and Durbin J. (2001). Time Series Analysis by State Space Methods. Oxford: Oxford University Press\cr
}

