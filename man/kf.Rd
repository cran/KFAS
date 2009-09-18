\name{kf}
\alias{kf}
\alias{filter}
\alias{Kalman}
\title{Kalman Filter With Exact Diffuse Initialisation}
\description{
  Performs Kalman filtering with exact diffuse initialisation using univariate approach (also known as sequential processing). 
  Written in Fortran, uses subroutines from BLAS and LAPACK.
  See function 'ks' for smoothing.
}
\usage{
  kf(yt, Zt, Tt, Rt, Ht, Qt, a1, P1, P1inf=0, optcal=c(TRUE,TRUE,TRUE,TRUE), tol=1e-7)
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
  \item{P1inf}{Diffuse part of P1. If non-zero, filtering starts with exact diffuse initialisation.}
  \item{optcal}{Vector of length 4. Calculate multivariate vt, Ft, Kt, Lt and their diffuse counterparts Finf, Fstar, Kinf, Kstar, Linf and Lstar. Default is c(TRUE,TRUE,TRUE,TRUE) which calculates all. Note that Kt cannot be calculated without Ft and Lt cannot be calculated without Kt, so even if optcal=c(TRUE,FALSE,TRUE,TRUE), Kt and Lt are not calculated.}
  \item{tol}{Tolerance parameter. Smallest covariance/variance value not counted for zero in diffuse phase. Default is 1e-7.}
}
\details{
  Function kf performs Kalman filtering of gaussian multivariate state space model using the univariate approach from Koopman and Durbin (2000, 2001). Univariate approach is also known as sequential processing, see Anderson and Moore (1979). In case where the distributions of some or all of the elements of initial state vector are not known, kf uses exact diffuse initialisation using univariate approach by Koopman and Durbin (2000, 2003). Note that in univariate approach the prediction error variance matrices Ft, Finf and Fstar does not need to be non-singular, as there is no matrix inversions in univariate approach algorithm. This provides much faster and more general filtering than normal multivariate Kalman filter algorithm.

  Filter can deal partially or totally missing observation vectors. If y_t,i is NA, it is interpreted as missing value, the dimensions of vtuni and Ktuni (or Kstaruni and Kinfuni) are decreased and the corresponding elements of vt, Kt and (and their diffuse counterparts) are marked as NA.

  The state space model is given by\cr

  y_t = Z_t * alpha_t + eps_t (observation equation)\cr
  alpha_t+1 = T_t * alpha_t + R_t * eta_t(transition equation)\cr

  where eps_t ~ N(0,H_t) and eta_t ~ N(0,Q_t) \cr
  Note that error terms eps_t and eta_t are assumed to be uncorrelated.

  When P1inf is non-zero, exact diffuse initialisation is used. Matrices Pt, Kt, Ft and Lt are decomposed to Pinf, Pstar, Kinf, Kstar, Finf, Fstar, Linf and Lstar.
  Diffuse phase is continued until Pinf becomes zero-matrix. See Koopman and Durbin (2000, 2001, 2003) for details for exact diffuse and non-diffuse filtering.

  Notice that vtuni, Ftuni, Fstaruni, Finfuni, Ktuni, Kinfuni and Kstaruni are not the same as calculated in usual multivariate Kalman filter. Also the Lt, Linf and Lstar are not calculated explicitly. If usual vt, Ft, Finf, Fstar, Kt, Kinf, Kstar, Lt, Linf and Lstar are needed, use optcal=c(TRUE, TRUE, TRUE, TRUE).
  When estimating parameters, it is suggested to use optcal=c(FALSE,FALSE,FALSE,FALSE) for maximum speed.

  Dimensions of variables are: \cr
  'yt' p*n\cr
  'Zt' p*m or p*m*n\cr
  'Tt' m*m or m*m*n\cr
  'Rt' m*r or m*r*n\cr
  'Ht' p*p or p*p*n\cr
  'Qt' r*r or r*r*n\cr
  'a1' m*1\cr
  'P1' m*m\cr
  'P1inf' m*m\cr

  where p is dimension of observation vector, m is dimension of state vector and n is number of observations.
}

\value{
  A list with the following elements:

  \item{yt}{p*n matrix or array of observations.}
  \item{ydimt}{array of length n which has dimensions of observation vector at times t=1,...,n.}
  \item{tv}{array of length 4 where tv[i]=0 if i is time-invariant and otherwise tv[i]=1, i is dt, Tt, Rt, Qt.}
  \item{Zt}{system matrix or array of observation equation.}
  \item{Tt}{system matrix or array of transition equation.}
  \item{Rt}{system matrix or array of transition equation.}
  \item{Ht}{variance matrix or array of disturbance terms eps_t of observation equation.}
  \item{Qt}{variance matrix or array of disturbance terms eta_t.}
  \item{a1}{initial state vector.}
  \item{P1}{variance matrix of a1. In diffuse case P1star, the non-diffuse part of P1.}
  \item{at}{m*(n+1) array of E(alpha_t | y_1, y_2, ... , y_t-1)}
  \item{Pt}{m*m*(n+1) array of Var(alpha_t | y_1, y_2, ... , y_t-1)}
  \item{vtuni}{p*1*n array of vt of univariate approach.} 
  \item{Ftuni}{Ft of univariate approach, Var(vtuni).} 
  \item{Ktuni}{m*p*n array of Kalman gain of univariate approach.} 
  \item{Pinf, Pstar}{p*p*d+1 arrays of diffuse phase decomposition of Pt.} 
  \item{Finfuni, Fstaruni}{p*p*d (p*d arrays of diffuse phase decomposition of Ftuni.}
  \item{Kinfuni, Kstaruni}{m*p*d arrays of diffuse phase decomposition of Ktuni.}
  \item{d}{the last index of diffuse phase, ie. the non-diffuse phase began from time d+1}
  \item{j}{the index of last y_t,i of diffuse phase.}
  \item{p}{the dimension of observation vector.}
  \item{m}{the dimension of state vector.}
  \item{r}{the dimension of variance matrix Qt.}
  \item{n}{the number of observations.}
  \item{lik}{Value of the log-likelihood function. If NaN, Ftuni_i,t was zero at some t=d+1,...,n or Finfuni_i,t and Fstaruni_i,t was zero at some t=1,...,d.}
  \item{optcal}{}
  \item{info}{if info != 0, Finf, Fstar or Ft was singular.}
  \item{vt}{p*1*n array of vt = yt - Zt * at.}
  \item{Ft}{p*p*n array of Ft = Var(vt) of Kalman filter.}
  \item{Kt}{m*p*n array of Kalman gain: Kt = Tt * Pt * Zt' * Ft^-1.}
 \item{Lt}{the m*m*n array, Lt = Tt - Kt * Zt.}
  \item{Finf, Fstar}{p*p*d arrays of diffuse phase decomposition of Ft.}
 \item{Kinf, Kstar}{m*p*d arrays of diffuse phase decomposition of Kt.}
  \item{Linf, Lstar}{m*m*d arrays of diffuse phase decomposition of Lt.}
  \item{tol}{Tolerance parameter.}

}

\references{
Koopman, S.J. and Durbin J. (2000). Fast filtering and smoothing for non-stationary time series models, Journal of American Statistical Assosiation, 92, 1630-38\cr

Koopman, S.J. and Durbin J. (2001). Time Series Analysis by State Space Methods. Oxford: Oxford University Press\cr

Koopman, S.J. and Durbin J. (2003). Filtering and smoothing of state vector for diffuse state space models, Journal of Time Series Analysis, Vol. 24, No. 1\cr

Shumway, Robert H. and Stoffer, David S. (2006). Time Series Analysis and Its Applications: With R examples\cr
}

\examples{
library(KFAS)

#Example of local level model
#Using the Nile observations
data(Nile)

yt<-t(data.matrix(Nile))
s2_eps<-15099
s2_eta<-1469.1

f.out<-kf(yt = yt, Zt = 1, Tt=1, Rt=1, Ht= s2_eps, Qt=s2_eta,
	a1 = 0, P1=1e7)

#a1 and P1 are not really estimated, should actually use exact diffuse initialisation:

fd.out<-kf(yt = yt, Zt = 1, Tt=1, Rt=1, Ht= s2_eps, Qt=s2_eta, a1 = 0,
	P1=0, P1inf=1)
#No stationary elements, P1=0, P1inf=1

#Plotting observations, non-diffuse and diffuse at, not plotting the a1=0:

ts.plot(Nile,ts(f.out$at[1,2:length(Nile)],start=1872,end=1970),
	ts(fd.out$at[1,2:length(Nile)],start=1872,end=1970),col=c(1,2,3))

#Looks identical. Actually start of series differs little bit:

f.out$at[1,1:20]
fd.out$at[1,1:20]


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

#Diffuse initialisation
#Notice that diffuse initialisation without univariate approach does not
#work here because Finf is non-singular and non-zero

est<-optim(par=c(.1,0,.1,.1), likfn, method="BFGS", control=list(fnscale=-1),
	hessian=TRUE, yt=yt, a1=0, P1=0, P1inf=1) 

pars<-est$par
L<-matrix(c(pars[1],pars[2],0,pars[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(pars[4])

kfd<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, 
	Ht=H, Qt=q11, a1 = 0, P1=0, P1inf=1)

#Same as non-diffuse, initial values from Shumway and Stoffer (2006):

est<-optim(par=c(.1,0,.1,.1), likfn, method="BFGS", control=list(fnscale=-1),
	hessian=TRUE, yt=yt, a1=-0.35, P1=0.01, P1inf=0)

pars<-est$par
L<-matrix(c(pars[1],pars[2],0,pars[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(pars[4])

kfnd<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, Ht=H, 
	Qt=q11, a1 = -0.35, P1=0.01, P1inf=0)

kfd$Qt
kfnd$Qt
kfd$Ht
kfnd$Ht
#Estimated Qt and Ht differs

ts.plot(ts(y1),ts(y2),ts(kfd$at[1,]),ts(kfnd$at[1,]),col=c(1,2,3,4))
#differs at start


#Example of stationary ARMA(1,1)

n<-1000
ar1<-0.8 
ma1<-0.3
sigma<-0.5
yt<-arima.sim(model=list(ar=ar1,ma=ma1), n=n, sd=sigma)

dt <- matrix(0, nrow = 2)
Zt<-matrix(c(1,0),ncol=2)
Tt<-matrix(c(ar1,0,1,0),ncol=2)
Rt<-matrix(c(1,ma1),ncol=1)
a1<-matrix(c(0,0),ncol=1)
P1<-matrix(0,ncol=2,nrow=2)
P1[1,1]<-1/(1-ar1^2)*(1+ma1^2+2*ar1*ma1)
P1[1,2]<-ma1
P1[2,1]<-ma1
P1[2,2]<-ma1^2
f.out<-kf(yt = array(yt,dim=c(1,n)), Zt = Zt, Tt=Tt, Rt=Rt, Ht= 0,
	Qt=sigma^2, a1 = a1, P1=P1)

}
