\name{forecast}
\alias{forecast}
\title{Forecast state space model}
\description{
  Performs forecasting using output from function 'kf' (Kalman filter).  
}
\usage{
  forecast(out,  fc=1,  Zt.fc=NULL,  Tt.fc=NULL,  Rt.fc=NULL,
Ht.fc=NULL,  Qt.fc=NULL)
}
\arguments{
\item{out}{Output from function 'kf'.  }
\item{fc}{Integer which states how many observations is forecasted.  }
\item{Zt.fc}{In case where matrix Z is not time-invariant, p*m*fc array
  of matrix Zt, t=n+1,...,n+fc.  }
\item{Tt.fc}{In case where matrix T is not time-invariant, m*m*fc array
  of matrix Tt, t=n+1,...,n+fc.  }
\item{Rt.fc}{In case where matrix R is not time-invariant, m*r*fc array
  of matrix Rt, t=n+1,...,n+fc.  }
\item{Ht.fc}{In case where matrix H is not time-invariant, p*p*fc array
  of matrix Ht, t=n+1,...,n+fc.  }
\item{Qt.fc}{In case where matrix Q is not time-invariant, r*r*fc array
  of matrix Qt, t=n+1,...,n+fc.  }
}

\details{
 
  The state space model is given by  

  y_t = Z_t * alpha_t + eps_t (observation equation)  


  alpha_t+1 = T_t * alpha_t + R_t * eta_t(transition equation)  


  where eps_t ~ N(0,H_t) and eta_t ~ N(0,Q_t).  


  Dimensions of variables are: \cr
  'yt' p*n\cr
  'Zt' p*m or p*m*n\cr
  'Tt' m*m or m*m*n\cr
  'Rt' m*r or m*r*n\cr
  'Ht' p*p or p*p*n\cr
  'Qt' r*r or r*r*n\cr
}
\value{
  A list with the following elements:

  \item{yt.fc}{p*fc array of forecasts of observations.  }
  \item{Ft.fc}{p*p*fc array of mean square error matrix.  }
  \item{at.fc}{m*(fc+1) array of E(alpha_t | y_1, y_2, ... , y_n).  }
  \item{Pt.fc}{m*m*(fc+1) array of Var(alpha_t | y_1, y_2, ... , y_n).  }
}
