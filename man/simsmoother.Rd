\name{simsmoother}
\alias{simsmoother}
\title{Simulation smoother}
\description{
  Generates random samples of state vector alpha from conditional density p(alpha | y).
}
\usage{
  simsmoother(yt, Zt, Tt, Rt, Ht, Qt, a1, P1, P1inf=0, nsim=1, tol=1e-7)
}
\arguments{
  \item{yt}{Matrix or array of observations.}
  \item{Zt}{System matrix or array of observation equation.}
  \item{Tt}{System matrix or array of transition equation.}
  \item{Rt}{System matrix or array of transition equation.}
  \item{Ht}{Variance matrix or array of disturbance terms eps_t of observation equation.}
  \item{Qt}{Variance matrix or array of disturbance terms eta_t.}
  \item{a1}{Initial state vector.}
  \item{P1}{Variance matrix of a1. In diffuse case P1star, the non-diffuse part of P1.}
  \item{P1inf}{Diffuse part of P1.}
  \item{nsim}{Number of samples. Default is 1.}  
  \item{tol}{Tolerance parameter. Smallest covariance/variance value not counted for zero in diffuse phase. Default is 1e-7.}
}
\details{
  Function simsmoother generates random samples of 
  state vector alpha = (alpha_1, ..., alpha_n) 
  from conditional density p(alpha | y).

  The state space model is given by\cr

  y_t = Z_t * alpha_t + eps_t (observation equation)\cr
  alpha_t+1 = T_t * alpha_t + R_t * eta_t(transition equation)\cr

  where eps_t ~ N(0,H_t) and eta_t ~ N(0,Q_t) \cr

  Simulation smoother algorithm is from article by J. Durbin and S.J. Koopman (2002).
}

\value{
 
\item{}{m*n*nsim array of simulated state vectors alpha.}
}

\references{
Durbin J. and Koopman, S.J. (2002). A simple and efficient simulation smoother for state space time series analysis, Biometrika, Volume 89, Issue 3\cr
}

\examples{

library(KFAS)

#Example of multivariate local level model
#Two temperature datas from David S. Stoffer's webpage

y1<-scan("http://www.stat.pitt.edu/stoffer/tsa2/data/HL.dat")
y2<-scan("http://www.stat.pitt.edu/stoffer/tsa2/data/folland.dat")
yt<-rbind(y1,y2)

#Estimating the variance parameters

likfn<-function(par, yt, a1, P1, P1inf) #Function to optim
{
L<-matrix(c(par[1],par[2],0,par[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(par[4])
lik<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, 
	Ht=H, Qt=q11, a1 = a1, P1=P1, P1inf=P1inf, optcal=c(FALSE,FALSE,FALSE,FALSE))
lik$lik
}

est<-optim(par=c(.1,0,.1,.1), likfn, method="BFGS", control=list(fnscale=-1),
	hessian=TRUE, yt=yt, a1=0, P1=0, P1inf=1) #Diffuse initialisation

pars<-est$par
L<-matrix(c(pars[1],pars[2],0,pars[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(pars[4])

kfd<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, 
	Ht=H, Qt=q11, a1 = 0, P1=0, P1inf=1)

ksd<-ks(kfd)

alpha<-simsmoother(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, 
	Ht=H, Qt=q11, a1 = 0, P1=0, P1inf=1, nsim=1000)

ahat<-NULL
for(i in 1:108) ahat[i]<-mean(alpha[1,i,])

ts.plot(ts(ahat),ts(ksd$ahat[1,]),col=c(1,2))

}
