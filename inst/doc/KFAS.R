## ----include=FALSE--------------------------------------------------
library(knitr)
opts_chunk$set(concordance=TRUE,tidy=FALSE)
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library(KFAS)

## ----'glmexample1'--------------------------------------------------
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
d.AD <- data.frame(treatment, outcome, counts)
glmModel1 <- SSModel(counts ~ outcome + treatment, 
                     data = d.AD, distribution = "poisson")

## ----'glmexample2'--------------------------------------------------
glmModel2 <- SSModel(counts ~ 
                       SSMregression(~outcome + treatment, data = d.AD), 
                     distribution = "poisson")

## ----'lmmexample'---------------------------------------------------
suppressWarnings(library("lme4",quietly=TRUE))
# Split data by grouping variable
Y <- split(sleepstudy["Reaction"], sleepstudy["Subject"])
p <- length(Y) # Number of series
Y <- matrix(unlist(Y), ncol = p, 
            dimnames = list(NULL, paste("Subject", names(Y))))
dataf <- split(sleepstudy, sleepstudy["Subject"])

# Assume that we know the covariance structure 
# of random effects and the residual variance

covRandom <- matrix(c(625,36,36,625),2,2)
sigma2 <- 650

# Common fixed effect part and distinct random effect parts for each "series"
# Set P1inf = 0 so diffuse initialization is not used for random effects
lmmModel <- SSModel(Y ~ -1 
                    + SSMregression(rep(list(~ Days), p), type = "common", 
                                    data = dataf, remove.intercept = FALSE)
                    + SSMregression(rep(list(~ Days), p), P1inf = 0,
                                    data = dataf, remove.intercept = FALSE),
                    H = diag(sigma2, p))
# Set covariance structure of the random effects which are states 3 to 38
# One could also use more common way lmmModel$P1[-(1:2),-(1:2)] <- ...
lmmModel["P1", 2+1:(2*p)] <- 
  as.matrix(.bdiag(replicate(p, covRandom, simplify = FALSE)))

## ----'alcoholPlot1', fig.pos = '!ht', fig.cap = 'Alcohol related deaths per 100,000 persons in Finland in 1969--2007.',out.width='\\linewidth'----
data("alcohol")
colnames(alcohol)
ts.plot(window(alcohol[,1:4]/alcohol[,5:8], end = 2007), col = 1:4,
        ylab = "Alcohol related deaths in Finland per 100,000 persons",
        xlab = "Year")
legend("topleft",col = 1:4, lty = 1,
       legend = colnames(alcohol)[1:4])

## ----'alcoholfit1'--------------------------------------------------
# remove the last observations
alcoholPred <- window(alcohol, start = 1969, end = 2007)

model <- SSModel(alcoholPred[,1:4] ~ 
                   SSMtrend(2, Q = list(matrix(NA,4,4), matrix(0,4,4))) +
                   SSMcustom(Z = diag(1,4), T = diag(0,4), 
                             Q = matrix(NA,4,4), P1 = matrix(NA,4,4)),
                 distribution = "poisson", u = alcoholPred[,5:8])

# Model updating function for fitSSM
updatefn <- function(pars, model, ...){
  Q <- diag(exp(pars[1:4]))
  Q[upper.tri(Q)] <- pars[5:10]
  model["Q",etas="level"] <- crossprod(Q)
  Q <- diag(exp(pars[11:14]))
  Q[upper.tri(Q)] <- pars[15:20]
  model["Q",etas=9:12] <- model["P1",states=9:12] <- crossprod(Q)
  model
}
# Initial the covariance structure of the random walks and extra noise
# theta = log(intensity) = log(y/u)
# covariance matrices are parameterized via log-Cholesky in fitSSM
init <- chol(cov(log(alcoholPred[,1:4]/alcoholPred[,5:8]))/10)

fitinit <- fitSSM(model, updatefn = updatefn, 
                  inits = rep(c(log(diag(init)), init[upper.tri(init)]),2),
                  method = "BFGS")
# Now with simulation
fit<-fitSSM(model, updatefn = updatefn, 
            inits = fitinit$optim.out$par, method = "BFGS", nsim = 250)
varcor <- fit$model["Q", etas = "level"]
varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
print(varcor,digits=2) #local level component
varcor <- fit$model["Q", etas = "custom"]
varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
print(varcor,digits=2) #local level component #extra noise component
-fitinit$optim.out$val #log-likelihood without simulation
-fit$optim.out$val     #log-likelihood with simulation  

## ----'KFS'----------------------------------------------------------
out <- KFS(fit$model, nsim = 1000)
out

## ----'states',fig.pos = '!ht', fig.cap = 'Smoothed level and white noise components.',out.width='\\linewidth'----
plot(coef(out,states=c("level","custom")), 
     main = "Smoothed states", yax.flip=TRUE)

## ----'diagnostics1',fig.pos = '!ht', fig.cap = 'Autocorrelations and cross-correlations of recursive residuals.',out.width='\\linewidth'----
res <- rstandard(KFS(fit$model, filtering = "mean", 
                     smoothing = "none", nsim = 1000))
acf(res, na.action = na.pass)

## ----'prediction'---------------------------------------------------
pred<-predict(fit$model, newdata = 
                SSModel(ts(matrix(NA,6,4), start = 2008) ~ -1
                        + SSMcustom(Z = fit$model$Z, T = fit$model$T, 
                                    R = fit$model$R, Q = fit$model$Q), u = 1, 
                        distribution= "poisson"),
              interval = "confidence", nsim = 10000)

## ----'predictplot',fig.pos = '!ht', fig.cap = 'Observed number of alcohol related deaths per 100,000 persons in Finland (black), fitted values (red) and intensity predictions for years 2008--2012 together with 95\\% prediction intervals (green).',out.width='\\linewidth'----
trend <- exp(signal(out, "trend")$signal)
par(mfrow = c(2,2), mar = c(2,2,2,2) + 0.1, oma = c(2,2,0,0))
for(i in 1:4)
  ts.plot(alcohol[,i]/alcohol[,4+i], 
          trend[,i], 
          pred[[i]],
          col = c(1,2,rep(3,3)), xlab = NULL, ylab = NULL, 
          main = colnames(alcohol)[i])
mtext("Number of alcohol related deaths per 100,000 persons in Finland", 
      side = 2, outer = TRUE)
mtext("Year",side=1,outer=TRUE)

