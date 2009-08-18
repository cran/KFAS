\name{estmiss}
\alias{estmiss}
\title{Estimate missing observations and mean square error matrix in state space model}
\description{
  Estimates yt and Ft (mean square error matrix) when yt is (partially) missing.
}
\usage{
  estmiss(out)
}
\arguments{
\item{out}{Output from function 'kf' or from 'ks'. See details.}
}

\details{

  Function 'estmiss' estimates mean square error matrices Ft, Fstar and Finf at time t where elements of yt is missing.
  If argument out is output from function 'ks', function also estimates missing observations by formula yt = Zt*ahat.
  
  Does not change values of Ft, Fstar, Finf or yt at times t where there is no elements of observation vector missing.
 
  The state space model is given by\cr

  y_t = Z_t * alpha_t + eps_t (observation equation)\cr
  alpha_t+1 = T_t * alpha_t + R_t * eta_t(transition equation)\cr

  where eps_t ~ N(0,H_t) and eta_t ~ N(0,Q_t) \cr

}
\value{
  A list with the following elements:

  \item{Fstar}{p*p*d array of non-diffuse part of mean square error matrix at diffuse phase.}
  \item{Finf}{p*p*d array of diffuse part of mean square error matrix at diffuse phase.}
  \item{Ft}{p*p*n array of mean square error matrix.}
  \item{yt}{p*n array of observations or estimates of observations.}
}
