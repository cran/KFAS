#' Extract Residuals of KFS output
#' @S3method residuals KFS
#' @method residuals KFS
#' @details For object of class KFS, several types of (possibly standardized) residuals can be computed:
#' \itemize{
#' 
#' \item 'recursive': One-step ahead prediction residuals. Standardized version is defined as 
#' \deqn{v_{t,i})/\sqrt{F_{i,t}},}
#' with residuals being undefined in diffuse phase. Only supported for fully Gaussian models.
#' 
#' \item 'response': Data minus fitted values, \eqn{y-E(y)}{y-E(y)}. 
#' 
#' \item 'pearson':  Without standardization, \deqn{(y_{t,i}-\theta_{t,i})/\sqrt{V(\mu)_{t,i}}, \quad i=1,\ldots,p,t=1,\ldots,n,}{(y[t,i]-\theta[t,i])V(\mu)[t,i]^(-0.5), i=1,\ldots,p, t=1,\ldots,n,}
#'                   where \eqn{V(\mu_{t,i})}{V(\mu[t,i])} is the variance function of the model. With standardization, above residuals are divided by 
#'                   \eqn{\sqrt{1-h_{t,i}}}{(1-h[t,i])^(0.5)}, where \eqn{h_{t,i}=(\frac{\textrm{d}\hat\mu}{\textrm{d}\theta})^2*V(\hat\mu)*V_\theta}{h[t,i]=(d\hat\mu/d\theta)^2*V(\hat\mu)*V[\theta]}. 
#'                   For gaussian models, these coincide with the smoothed \eqn{\epsilon} disturbance residuals.
#'
#' \item 'state':  Residuals based on the smoothed disturbance terms \eqn{\eta} are defined as
#' \deqn{L^{-1}_t \hat \eta_t, \quad t=1,\ldots,n,}{L^{-1}[t] \eta[t], t=1,\ldots,n,} where \eqn{L_t}{L[t]} is the lower triangular matrix from Cholesky decomposition of \eqn{V_{\eta,t}}{V[\eta,t]}.
#' 
#' \item 'deviance': Deviance residuals.
#' }
#'
#' @param object KFS object
#' @param type character string defining the type of residuals. See details.
#' @param standardize logical. Should the residuals be standardized? Not applicable for type \code{'response'}.
#' @param ... Ignored.

residuals.KFS <- function(object, type = c("recursive", "response", "pearson", "state", "deviance"), standardize = TRUE, ...) {
    
    type <- match.arg(type)
    
    if ((type == "recursive" || type == "state") && any(object$model$distribution != "gaussian")) 
        stop("Recursive and state residuals are only supported for fully gaussian models.")
    
    
    recursive <- function(object, standardize) {
        if (standardize) {
            series <- object$v/sqrt(t(object$F))
        } else series <- object$v
        series[1:(object$d - 1), ] <- NA
        series[object$d, 1:object$j] <- NA
        series
    }
    
    response <- function(object, ...) {
        if (all(object$model$distribution == "gaussian")) {
            if (is.null(object$thetahat)) 
                stop("KFS object needs to contain smoothed estimates of signal.")
            series <- object$model$y - object$thetahat
            
        } else {
            if (is.null(object$muhat)) 
                stop("KFS object needs to contain smoothed estimates of mu")
            series <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
            for (i in 1:attr(object$model, "p")) series[, i] <- switch(object$model$distribution[i], gaussian = (object$model$y[, 
                i] - object$muhat[, i]), poisson = (object$model$y[, i] - object$model$u[, i] * object$muhat[, i]), binomial = (object$model$y[, 
                i] - object$model$u[, i] * object$muhat[, i]), gamma = (object$model$y[, i] - object$muhat[, i]), `negative binomial` = (object$model$y[, 
                i] - object$muhat[, i]))
            
        }
        
        series
    }
    
    pearson <- function(object, standardize) {
        if (all(object$model$distribution == "gaussian")) {
            if (is.null(object$thetahat)) 
                stop("KFS object needs to contain smoothed estimates of signal.")
            if (standardize) {
                series <- (object$model$y - object$thetahat)/
                  matrix(sqrt(apply(array(object$model$H,c(attr(object$model, "p"),attr(object$model, "p"),attr(object$model, "n"))) 
                                    - object$V_theta, 3, diag)), 
                  attr(object$model, "n"), attr(object$model, "p"), byrow = TRUE)
            } else series <- object$model$y - object$thetahat
            
        } else {
            if (is.null(object$muhat)) 
                stop("KFS object needs to contain smoothed estimates of mu.")
            
            series <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
            for (i in 1:attr(object$model, "p")) series[, i] <- switch(object$model$distribution[i], gaussian = (object$model$y[, 
                i] - object$muhat[, i]), poisson = (object$model$y[, i] - object$model$u[, i] * object$muhat[, i])/sqrt(object$model$u[, 
                i] * object$muhat[, i]), binomial = (object$model$y[, i] - object$model$u[, i] * object$muhat[, i])/sqrt(object$model$u[, 
                i] * object$muhat[, i] * (1 - object$muhat[, i])), gamma = (object$model$y[, i] - object$muhat[, i])/object$muhat[, i], `negative binomial` = (object$model$y[, 
                i] - object$muhat[, i])/sqrt(object$muhat[, i] + object$muhat[, i]^2/object$model$u[, i]))
                      
            if (standardize) {
                w <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
                for (i in 1:attr(object$model, "p")) w[, i] <- switch(object$model$distribution[i], gaussian = 1, poisson = object$model$u[, 
                  i] * object$muhat[, i], binomial = object$model$u[, i] * object$muhat[, i] * (1 - object$muhat[, i]), gamma = 1, `negative binomial` = (object$model$u[, 
                  i] * object$muhat[, i])/(object$model$u[, i] + object$muhat[, i]))
                
                series <- series/sqrt((1 - w * matrix(apply(object$V_theta, 3, diag), attr(object$model, "n"), attr(object$model, 
                  "p"), byrow = TRUE)))
            }
            
        }
        
        series
    }
    
    
    state <- function(object, standardize) {
        if (is.null(object$etahat)) {
            stop("KFS object needs to contain smoothed estimates of state disturbances eta.")
        } else {
            if (standardize) {
                k <- attr(object$model, "k")
                n <- attr(object$model, "n")
                if (dim(object$model$Q)[3] == 1) {
                  z <- which(object$model$Q[, , 1][1 + 0:(k - 1) * (k + 1)] > 0)
                  eta <- array(0, c(n, length(z)))
                  for (i in 1:(n - 1)) {
                    if (!isTRUE(all.equal(object$etahat[i, z], rep(0, length(z))))) {
                      x <- try(chol(solve(object$model$Q[z, z, 1] - object$V_eta[z, z, i])) %*% object$etahat[i, z], TRUE)
                      if (inherits(x, "try-error")) {
                        warning(paste("Could not compute the standardized smoothed state residuals, V_eta[,,", i, "] is not invertible", 
                          sep = ""))
                        break
                      } else eta[i, ] <- x
                    }
                  }
                  
                } else {
                  z <- NULL
                  for (i in 1:k) if (sum(object$model$Q[i, i, ]) > 0) 
                    z <- c(z, i)
                  zlength <- length(z)
                  eta <- array(NA, c(n, zlength))
                  if (zlength > 1) {
                    for (i in 1:(n - 1)) {
                      if (!isTRUE(all.equal(object$etahat[i, z][z2], rep(0, length(z2))))) {
                        z2 <- which(object$V_eta[z, z, i][1 + 0:(zlength - 1) * (zlength + 1)] > 0)
                        x <- try(chol(solve(object$model$Q[z, z, i][z2, z2] - object$V_eta[z, z, i][z2, z2])) %*% object$etahat[i, 
                          z][z2], TRUE)
                        if (inherits(x, "try-error")) {
                          warning(paste("Could not compute the standardized smoothed state residuals, V_eta[,,", i, "] is not invertible", 
                            sep = ""))
                          break
                        } else eta[i, z2] <- x
                      }
                    }
                  } else {
                    for (i in 1:n) {
                      if (!isTRUE(all.equal(object$etahat[i, z], rep(0, length(z))))) 
                        eta[i, 1] <- object$etahat[i, z]/sqrt(object$model$Q[z, z, i] - object$V_eta[z, z, i])
                    }
                  }
                }
                eta[n, ] <- 0
                return(eta)
            } else return(object$etahat)
        }
        
    }
    
    deviance <- function(object, standardize) {
      if (all(object$model$distribution == "gaussian")) {
        if(is.null(object$thetahat))
          stop("KFS object needs to contain smoothed estimates of signal")
        series <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
        for (i in 1:attr(object$model, "p")) series[, i] <- ifelse(object$model$y[,i] > 
                 object$thetahat[, i], 1, -1) * abs(object$model$y[, i] - object$thetahat[, i])
        if (standardize) {
          series <- series/sqrt(matrix(c(apply(object$model$H, 3, diag)) - apply(object$V_theta, 3, diag), attr(object$model, "n"), attr(object$model, "p"), 
                                                byrow = TRUE))
        }                                  
      } else {
      if (is.null(object$muhat)) 
        stop("KFS object needs to contain smoothed estimates of mu.")     
      series <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
        for (i in 1:attr(object$model, "p")) series[, i] <- switch(object$model$distribution[i], gaussian = ifelse(object$model$y[, 
            i] > object$muhat[, i], 1, -1) * sqrt(abs(object$model$y[, i] - object$muhat[, i])), poisson = ifelse(object$model$y[, i] > (object$model$u[, 
            i] * object$muhat[, i]), 1, -1) * sqrt(2 * (object$model$y[, i] * log(ifelse(object$model$y[, i] == 0, 1, object$model$y[, 
            i]/(object$muhat[, i] * object$model$u[, i]))) - (object$model$y[, i] - object$model$u[, i] * object$muhat[, i]))), binomial = ifelse(object$model$y[, 
            i] > (object$model$u[, i] * object$muhat[, i]), 1, -1) * sqrt(2 * (object$model$y[, i] * log(ifelse(object$model$y[, i] == 
            0, 1, object$model$y[, i]/object$muhat[, i])) + (object$model$u[, i] - object$model$y[, i]) * log((object$model$u[, i] - 
            object$model$y[, i])/(object$model$u[, i] - object$muhat[, i])))), gamma = ifelse(object$model$y[, i] > object$muhat[, i], 
            1, -1) * sqrt(-2 * (log(ifelse(object$model$y[, i] == 0, 1, object$model$y[, i]/object$muhat[, i])) - (object$model$y[, 
            i] - object$muhat[, i])/object$muhat[, i])), `negative binomial` = ifelse(object$model$y[, i] > object$muhat[, i], 1, -1) * 
            sqrt(2 * (object$model$y[, i] * log(pmax(1, object$model$y[, i])/object$muhat[, i], 1) - (object$model$y[, i] + object$model$u[, 
                i]) * log((object$model$y[, i] + object$model$u[, i])/(object$muhat[, i] + object$model$u[, i])))))
        if (standardize) {
            w <- matrix(0, attr(object$model, "n"), attr(object$model, "p"))
            for (i in 1:attr(object$model, "p")) w[, i] <- switch(object$model$distribution[i], gaussian = 1, poisson = object$model$u[, 
                i] * object$muhat[, i], binomial = object$model$u[, i] * object$muhat[, i] * (1 - object$muhat[, i]), gamma = 1, `negative binomial` = (object$model$u[, 
                i] * object$muhat[, i])/(object$model$u[, i] + object$muhat[, i]))
            
            series <- series/sqrt((1 - w * matrix(apply(object$V_theta, 3, diag), attr(object$model, "n"), attr(object$model, "p"), 
                byrow = TRUE)))
        }
        
    }
      series
  }
    
    return(ts(drop(do.call(type, list(object, standardize))),start=start(object$model$y),frequency=frequency(object$model$y),names=dimnames(object$model$y)[[2]]))
} 
