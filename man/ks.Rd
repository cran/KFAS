\name{ks}
\alias{ks}
\title{Kalman Smoother With Exact Diffuse Initialisation}
\description{
  Performs Kalman smoothing with exact diffuse phase using univariate
  approach.  Programmed with Fortran, uses subroutines from BLAS and
  LAPACK.  See function 'kf' for smoothing.  
}
\usage{
  ks(out)
}
\arguments{
  \item{out}{Output from function 'kf'}
}
\details{
  Function ks performs Kalman smoothing of gaussian (multivariate) state
  space model using the univariate approach by Koopman and Durbin (2000,
  2001, 2003).  In case where the distributions of some or all elements
  of initial state vector are unknown, ks uses exact diffuse phase using
  univariate approach by Koopman and Durbin (2000).  
}

\value{
  A list with the output elements from Kalman filter and following new elements:

\item{ahat}{m*n array of E(alphat | y_1, y_2, ... , y_n).  }
\item{Vt}{m*m*n array of Var(alphat|y_1, y_2, ... , y_n).  } 
\item{rt}{m*n+1 array of weighted sums of innovations vj, j=t+1,...n.
  Notice that in literature t in rt goes from 0,..n. Here
  t=1,...n+1. Same applies to all r and N variables.  }  
\item{rt0, rt1}{m*d+1 arrays of diffuse phase decomposition of rt.  }  
\item{Nt}{m*m*n+1 array of of Var(rt).  }  
\item{Nt0, Nt1, Nt2}{m*m*d+1 arrays of diffuse phase decomposition of
  Nt.  }
}

\references{
Koopman, S.J. and Durbin J. (2000).  Fast filtering and smoothing for
non-stationary time series models, Journal of American Statistical
Assosiation, 92, 1630-38.  \cr

Koopman, S.J. and Durbin J. (2001). Time Series Analysis by State Space
Methods. Oxford: Oxford University Press.  \cr

Koopman, S.J. and Durbin J. (2003). Filtering and smoothing of state
vector for diffuse state space models, Journal of Time Series Analysis,
Vol. 24, No. 1.  \cr

Shumway, Robert H. and Stoffer, David S. (2006). Time Series Analysis
and Its Applications: With R examples.  \cr
}

\examples{

library(KFAS)

#Example of multivariate local level model
#Two series of average global temperature deviations for years 1880-1987
#See Shumway and Stoffer (2006), p. 327 for details

data(GlobalTemp)
yt<-array(GlobalTemp,c(2,108))

#Estimating the variance parameters

likfn<-function(par, yt, a1, P1, P1inf) #Function to optim
{
L<-matrix(c(par[1],par[2],0,par[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(par[4])
lik<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, Ht=H, Qt=q11, a1 =
a1, P1=P1, P1inf=P1inf, optcal=c(FALSE,FALSE,FALSE,FALSE))
lik$lik
}

est<-optim(par=c(.1,0,.1,.1), likfn, method="BFGS",
control=list(fnscale=-1), hessian=TRUE, yt=yt, a1=0, P1=0, P1inf=1) 
#Diffuse initialisation

pars<-est$par
L<-matrix(c(pars[1],pars[2],0,pars[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(pars[4])

kfd<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, 
	Ht=H, Qt=q11, a1 = 0, P1=0, P1inf=1)

ksd<-ks(kfd)

#Same as non-diffuse, initial values from Shumway and Stoffer (2006):

est<-optim(par=c(.1,0,.1,.1), likfn, method="BFGS",
control=list(fnscale=-1), hessian=TRUE, yt=yt, a1=-0.35, P1=0.01, P1inf=0 )

pars<-est$par
L<-matrix(c(pars[1],pars[2],0,pars[3]),ncol=2)
H<-L\%*\%t(L)
q11<-exp(pars[4])

kfnd<-kf(yt = yt, Zt = matrix(1,nrow=2), Tt=1, Rt=1, Ht=H, Qt=q11, a1 =
-0.35, P1=0.01, P1inf=0)

ksnd<-ks(kfnd)

kfd$Qt
kfnd$Qt
kfd$Ht
kfnd$Ht
#Estimated Qt and Ht differs

#smoothed values:
ts.plot(ts(yt[1,],start=1880),ts(yt[2,],start=1880),ts(ksd$ahat[1,],start=1880),ts(ksnd$ahat[1,],start=1880),col=c(1,2,3,4))

}
