\name{distsmoother}
\alias{distsmoother}
\title{Disturbance smoother}
\description{
  Computes smoothed values of disturbance terms eps_t and eta_t, and their variances.
}
\usage{
  distsmoother(out)
}
\arguments{
  \item{out}{Output of function 'ks'. Optcal must be (TRUE,TRUE,TRUE,TRUE) or (TRUE,TRUE,TRUE,FALSE).}
}

\value{
A list with the output elements from Kalman smoother and following new elements: 
\item{epshat}{p*n array of E(eps_t | y_1, ..., y_n).}
\item{epshatvar}{p*p*n array of var(eps_t | y_1, ..., y_n).}
\item{etahat}{r*n array of E(eta_t | y_1, ..., y_n).}
\item{etahatvar}{r*r*n array of var(eta_t | y_1, ..., y_n).}
}

\references{
Koopman, S.J. and Durbin J. (2001). Time Series Analysis by State Space Methods. Oxford: Oxford University Press\cr
}

