\name{eflik0}
\alias{eflik0}
\title{Approximate log-likelihood and gaussian density of exponential family state-space model.}
\description{
  Function eflik0 computes approximate log-likelihood and approximate gaussian density of univariate exponential family state-space model, based on Durbin and Koopman (1997, 2001).
}
\usage{
  eflik0(yt, Zt, Tt, Rt, Qt, a1, P1, P1inf, dist=c("Poisson", "Binomial", "Negative binomial"), offset=1)
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


\examples{

kk <- c(396,  399,  403,  434,  340,  519,  558,  566,  591,  566,  574,  646,  644,  646,  665,
693,  773,  834,  910, 1035, 1002, 1161, 1056, 1097, 1094, 1042, 1194, 1316, 1246, 1503,
1428, 1477, 1490, 1465, 1560, 1860, 2008, 2020, 2167)

vlk <- c(4623785, 4606307, 4612124, 4639657, 4666081, 4690574, 4711440, 4725664, 4738902,
4752528, 4764690, 4779535, 4799964, 4826933, 4855787, 4881803, 4902206, 4918154,
4932123, 4946481, 4964371, 4986431, 5013740, 5041992, 5066447, 5088333, 5107790,
5124573, 5139835, 5153498, 5165474, 5176209, 5188008, 5200598, 5213014, 5228172,
5246096, 5266268, 5288720)


n<-39
Zt<-array(c(1,0),c(1,2))
Tt<-diag(c(1,1))
Tt[1,2]<-1
Rt <- diag(2)
a1 <- array(0,dim=2)
P1 <- diag(0,2)
P1inf <-diag(1,2)

likfn<-function(pars,yt,Zt,Tt,Rt,a1,P1,P1inf,dist,offset=1)
{
Qt<-diag(exp(pars))
lik<-eflik0(yt,Zt,Tt,Rt,Qt,a1,P1,P1inf,dist,offset)$lik0
lik
}

#not run

#opt <- optim(par=c(1,1),yt=array(kk,dim=c(1,39)),Zt=Zt,Tt=Tt,Rt=Rt,a1=a1,P1=P1,P1inf=P1inf,offset=vlk,dist="Poisson",fn=likfn,method="BFGS", control=list(fnscale=-1,trace=1,REPORT=1))

#Qt<-diag(exp(opt$par))
#out<-eflik0(yt=kk,Zt,Tt,Rt,Qt,a1,P1,P1inf,offset=vlk)
#out$Qt*10000
}

