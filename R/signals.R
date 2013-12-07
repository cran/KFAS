#' Extracting the Partial Signal Of a State Space Model
#'
#' Function \code{signal} returns the signal of a state space model using only subset of states.
#'
#' If the 
#' @export
#' @param object Object of class \code{KFS}.
#' @param states Which states are combined? Either a numeric vector containing the indices of the corresponding states,
#' or a character vector defining the types of the corresponding states. 
#' Possible choices are \dQuote{all}, \dQuote{arima}, \dQuote{custom}, \dQuote{cycle}, \dQuote{seasonal}, 
#' \dQuote{trend}, or \dQuote{regression}. These can be combined. Default is \dQuote{all}.
#' @param filtered if TRUE, filtered signal is used. Otherwise smoothed signal is used.
#' @return
#'\item{signal}{Time series object of filtered signal \eqn{Z_ta_t}{Z[t]a[t]} or smoothed signal \eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]} using only the defined states.  }
#' \item{variance}{Cov(\eqn{Z_ta_t}{Z[t]a[t]}) or Cov(\eqn{Z_t\hat\alpha_t}{Z[t]\alpha[t]}) using only the defined states. 
#' For the covariance matrices of the filtered signal, only the non-diffuse part P is used.  }
signal <- function(object, states = "all", filtered = FALSE) {
    if(!inherits(object,"KFS"))
      stop("Object must be an output from function KFS.")
    if (is.numeric(states)) {
        states <- as.integer(states)
        if (min(states) < 1 | max(states) > attr(object$model, "m")) 
            stop("Vector states should contain the indices or names of the states which are combined.")
    } else {
        states <- match.arg(arg = states, choices = c("all", "arima", "custom", "cycle", "seasonal", "trend", "regression"), several.ok = TRUE)
        if ("all" %in% states) {
            states <- as.integer(1:attr(object$model, "m"))
        } else states <- which(attr(object$model, "state_types") %in% states)
    }
    
    if (!filtered && !is.null(object$theta) && identical(states, as.integer(1:attr(object$model, "m")))) 
        return(list(signal = object$theta, variance = object$V_theta))
    if (!isTRUE(length(states) > 0)) 
        stop("Selected states not in the model.")
    
    if (filtered) {
        a <- object$a
        P <- object$P       
    } else {
        if (is.null(object$alphahat)) {
            call_KFS <- object$call
            call_KFS$smoothing <- "states"
            call_KFS$model <- object$model
            object <- eval(call_KFS)
        }
        
        a <- object$alphahat
        P <- object$V        
    }
    
    signal <- .Fortran(fsignaltheta, NAOK = TRUE, as.integer(dim(object$model$Z)[3] > 1), 
                       object$model$Z, t(a)[1:attr(object$model, "m"), 1:attr(object$model, "n")], 
                       P[1:attr(object$model, "m"), 1:attr(object$model, "m"), 1:attr(object$model, "n")], 
                       as.integer(attr(object$model, "p")), as.integer(attr(object$model, "n")), 
                       as.integer(attr(object$model, "m")), 
                       theta = array(0, c(attr(object$model, "n"), attr(object$model, "p"))), 
                       V_theta = array(0, c(attr(object$model, "p"), attr(object$model, "p"), attr(object$model, "n"))), 
                       d = as.integer(0), states, as.integer(length(states)))
    # if (d > 0) { signal$theta[1:d, ] <- NA signal$V_theta[, , 1:d] <- NA }
    
    
    attributes(signal$theta) <- attributes(object$model$y)
    list(thetahat = signal$theta, V_theta = signal$V_theta)
} 
