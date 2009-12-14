\name{efsmoother}
\alias{efsmoother}
\title{Kalman Smoother for exponential family models}
\description{
  Performs Kalman smoothing for exponential family models via simulation.
}
\usage{
  efsmoother(out,nsim)
}
\arguments{
  \item{out}{Output from function 'expflik'.}
  \item{nsim}{Number of simulations used for smoothing.}
}
\details{
  Function efsmoother performs Kalman smoothing of univariate exponential family state space model using the simulation algorithms by Koopman and Durbin (2001, 2002).
}

\value{
  A list with the output elements from expflik and following new elements:
\item{nsim}{Number of simulations used for smoothing.}

\item{ahat}{m*n array of E(alphat | y_1, y_2, ... , y_n).}
\item{that}{m*n array of Zt*ahat.}
\item{ahatvar}{m*m*n array of Var(ahat|y_1, y_2, ... , y_n).}
\item{thatvar}{1*1*n array of Zt*Var(ahat|y_1, y_2, ... , y_n)*Zt'.}
}

\references{

Koopman, S.J. and Durbin J. (2001). Time Series Analysis by State Space Methods. Oxford: Oxford University Press\cr

Durbin J. and Koopman, S.J. (2002). A simple and efficient simulation smoother for state space time series analysis, Biometrika, Volume 89, Issue 3\cr
}

\examples{

library(KFAS)



kk <- c(396,  399,  403,  434,  340,  519,  558,  566,  591,  566,  574,  646,  644,  646,  665,
693,  773,  834,  910, 1035, 1002, 1161, 1056, 1097, 1094, 1042, 1194, 1316, 1246, 1503,
1428, 1477, 1490, 1465, 1560, 1860, 2008, 2020, 2167)

vlk <- c(4623785, 4606307, 4612124, 4639657, 4666081, 4690574, 4711440, 4725664, 4738902,
4752528, 4764690, 4779535, 4799964, 4826933, 4855787, 4881803, 4902206, 4918154,
4932123, 4946481, 4964371, 4986431, 5013740, 5041992, 5066447, 5088333, 5107790,
5124573, 5139835, 5153498, 5165474, 5176209, 5188008, 5200598, 5213014, 5228172,
5246096, 5266268, 5288720)

#Model

n <- 39
m <- 2
r <- 2

yt <- array(kk,dim=c(1,n))

offset <- array(vlk,dim=c(1,n))

Zt <- array(c(1,0),dim=c(1,m,n))

Tt <- array(c(1,0,1,1),dim=c(m,m,n))

Rt <- array(c(1,0,0,1),dim=c(m,r,n))

psi_hat <- c(-5.356458, -16.348376)

Qt <- array(c(exp(psi_hat[1]),0,0,exp(psi_hat[2])),dim=c(r,r,n))

a1 <- matrix(0,nrow=m,ncol=1)

P1 <- matrix(0,nrow=m,ncol=m)

P1inf<-diag(1,2)

dist <- "Poisson"

eflikout <- eflik0(yt,Zt,Tt,Rt,Qt,a1,P1,P1inf,dist,offset)

nsim <- 1000

smef <- efsmoother(eflikout,nsim)


}
