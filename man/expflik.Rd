\name{expflik}
\alias{expflik}
\title{Approximate log-likelihood and gaussian density of exponential family state-space model.}
\description{
  Function expflik computes approximate log-likelihood and approximate gaussian density of exponential family state-space model, based on article of J. Durbin and S.J. Koopman (1997).
}
\usage{
  expflik(yt, Zt, Tt, Rt, Qt, a1, P1, P1inf,dist="Poisson",offset=1)
}
\arguments{
  \item{yt}{Matrix or array of observations.}
  \item{Zt}{System matrix or array of observation equation.}
  \item{Tt}{System matrix or array of transition equation.}
  \item{Rt}{System matrix or array of transition equation.}
  \item{Qt}{Variance matrix or array of disturbance terms eta_t.}
  \item{a1}{Initial state vector.}
  \item{P1}{Variance matrix of a1. In diffuse case P1star, the non-diffuse part of P1.}
  \item{P1inf}{Diffuse part of P1.}
  \item{dist}{Distribution of yt. Currently only "Poisson" is supported.}
  \item{offset}{''Intercept'' part of Poisson distribution. If theta=log(lambda), p(y|theta)= exp[y|theta - offset*exp(theta) + y*log(offset) - log(y!)].}
}

\details{
  
  Function approximates p(alpha|y) with gaussian g(alpha|y) which has same conditional mode (alpha,alpha_n+1|y) as p(alpha|y), and computes approximate log-likelihood\cr
  logL_0(psi) = logL_g(psi) + logE_g(w_T|y),\\


where logL_g(psi) is log-likelihood of approximate gaussian distribution and logE_g(w_T|y) is a Taylor-approximation of E_g[p(y|theta)/g(y|theta) | y].
For details, see J. Durbin and S.J. Koopman (1997).


  The general state space model for exponential family is given by\cr

  p(y_t|theta_t) = exp[theta'_t * y_t - b_t(theta_t) + f_t(y_t)] (observation equation)\cr
  alpha_t+1 = T_t * alpha_t + R_t * eta_t (transition equation)\cr

  where theta_t = Z_t * alpha_t and eta_t ~ N(0,Q_t) \cr
 
  Approximating gaussian model is given by\cr
  y*_t = Z_t * alpha_t + eps_t (observation equation)\cr
  alpha_t+1 = T_t * alpha_t + R_t * eta_t (transition equation)\cr

  where eps_t ~ N(0,H*_t) and eta_t ~ N(0,Q_t) \cr
  
}

\value{
List with output from Kalman smoother and distribution smoother, when model is approximated with gaussian distribution g(y|theta). Note that Ht is H*_t. List also contains following items:
\item{ytilde}{y* of approximating gaussian model.}
\item{theta}{Z_t * alpha_t of approximating gaussian model.}
\item{lik0}{Value of logL_0.}
\item{offset}{''Intercept'' part of Poisson distribution.}
}

\references{
Durbin J. and Koopman, S.J. (1997). Monte Carlo Maximum Likelihood Estimation for Non-Gaussian State Space Models, Biometrica, Vol. 84, No. 3.\cr
}

