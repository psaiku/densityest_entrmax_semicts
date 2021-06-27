library(EnvStats)

# function code is from https://rstudio-pubs-static.s3.amazonaws.com/205225_5991ce162f504576a84eac9c659aeaaf.html
newton.raphson <- function(f, a, b, tol = 1e-5, n = 1000) {
  require(numDeriv) # Package for computing f'(x)
  
  x0 <- a # Set start value to supplied lower bound
  k <- n # Initialize for iteration results
  
  # Check the upper and lower bounds to see if approximations result in 0
  fa <- f(a)
  if (fa == 0.0) {
    return(a)
  }
  
  fb <- f(b)
  if (fb == 0.0) {
    return(b)
  }
  
  for (i in 1:n) {
    dx <- genD(func = f, x = x0)$D[1] # First-order derivative f'(x0)
    #dx <- derf2(x0)
    x1 <- x0 - (f(x0) / dx) # Calculate next value x1
    if(x1 < a | x1 > b) {
      x1 <- runif(1)
    }
    k[i] <- x1 # Store x1
    # Once the difference between x0 and x1 becomes sufficiently small, output the results.
    if (abs(x1 - x0) < tol) {
      root.approx <- tail(k, n=1)
      res <- list('root approximation' = root.approx, 'iterations' = k)
      return(res)
    }
    # If Newton-Raphson has not yet reached convergence set x1 as x0 and continue
    x0 <- x1
  }
  print('Too many iterations in method')
}

est_mle_tpgamma <- function(y) {
  idx <- which(y>0)
  g <- sum(y==0)/length(y)
  mle_g <- egamma(y[idx])
  k <- mle_g$parameters[1]
  thet <- mle_g$parameters[2]
  entr_g <- k + log(thet) +log(gamma(k)) + (1-k)*digamma(k)
  entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
  return(list("Gamma"=g, "Shape"=k, "Scale"=thet, "MaxEntropy"=entr))
}

est_mle_tpexp <- function(y) {
  idx <- which(y>0)
  g <- sum(y==0)/length(y)
  lambda <- 1/mean(y[idx])
  entr_g <- 1 - log(lambda)
  entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
  return(list("Gamma"=g, "Lambda"=lambda, "MaxEntropy"=entr))
}

est_mle_tplnorm <- function(y) {
  idx <- which(y>0)
  g <- sum(y==0)/length(y)
  mle_g <- elnorm(y[idx], method = "mle")
  meanlog <- mle_g$parameters[1]
  sdlog <- mle_g$parameters[2]
  #entr_g <- sdlog*exp(meanlog + 0.5)*sqrt(2*pi)
  entr_g <- meanlog + 0.5*log(2*pi*exp(1)*sdlog^2)
  entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
  return(list("Gamma"=g, "mean"=meanlog, "sd"=sdlog, "MaxEntropy"=entr))
}

hist.semicts <- function(obj, xlab="X", ylab="Density", main="",
                         cols=c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), legends=c("Intensity", "Proportion of Zeros")) {
  obs.zero.prop <- length(obj[obj == 0]) / length(obj)
  obs.hist <- hist.default(obj, breaks = seq(0, max(obj)+1, by = 0.5), plot = FALSE)
  obs.hist$counts <- obs.hist$counts / sum(obs.hist$counts)
  plot(obs.hist, col = cols[1], xlab = xlab, ylab = ylab, main = main)
  points(x = c(0), y = c(obs.zero.prop), col = cols[2], pch = 19, cex = 1.5)
  #legend("topright", legends, cex = 1.5, lwd = c(10, NA), col = cols, pch = c(NA, 19))
}