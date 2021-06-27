library(semicts)
library(data.table)
library(ggplot2)
library(rootSolve)
library(pracma)

source("semicts_est_utils.R")

#
# Algorithm 1
# Arguments: 
# y - data
# max.iter - maximum number of iterations allowed; tol - stopping criterion tolerance
# gmma
est_scem_tpgamma_old <- function(y, max.iter = 100, tol = 1e-6, g = 0.4) {
  idx <- which(y>0)
  mean_x <- 0
  mean_logx <- 0
  entr_old <- 0
  entrs <- c()
  gammas <- c()
  Hgs <- c()
  
  f1 <- function(x) {
    mean_logx - log(x) - digamma(mean_x / x)
  }
  
  derf1 <- function(x) {
    -(1/x) + trigamma(mean_x / x) * (1/x^2)
  }
  
  for(i in 1:max.iter) {
    print(g)
    gammas[i] <- g
    mean_x <- mean(y)/(1-g)
    mean_logx <- sum(log(y[idx]))/(length(y) * (1-g))
    #thet <- newton.raphson(f1, 0.001, 1000)$`root approximation`
    thet <- uniroot(f1, c(0.00001, 10000))$root
    #thet <- uniroot.all(f1, c(0.00001, 10000))[1]
    k <- mean_x / thet
    entr_g <- k + log(thet) +log(gamma(k)) + (1-k)*digamma(k)
    Hgs[i] <- entr_g
    g <- 1/(1+exp(entr_g))
    entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
    entrs[i] <- entr
    #print(sprintf("Iter: %d, Gamma: %f, Shape: %f, Scale: %f, Entropy: %f", i, g, k, thet, entr))
    if(abs(entr-entr_old) < tol) {
      break
    } else {
      entr_old <- entr
    }
  }
  #print(sprintf("Gamma: %f, Shape: %f, Scale: %f", g, k, thet))
  return(list("Gamma"=g, "Shape"=k, "Scale"=thet, "Gammas" = gammas, "Hgs"=Hgs,"MaxEntropy"=entr, "EntropyValues"=entrs, "IterUsed"=i))
}

est_em_tpgamma <- function(y) {
  idx <- which(y>0)
  mean_x <- 0
  mean_logx <- 0
  entr_old <- 0
  entrs <- c()
  gammas <- c()
  Hgs <- c()
  g <- sum(y==0)/length(y)
  
  f1 <- function(x) {
    mean_logx - log(x) - digamma(mean_x / x)
  }
  
  #mean_x <- mean(y)/(1-g)
  #mean_logx <- sum(log(y[idx]))/(length(y) * (1-g))
  mean_x <- mean(y[idx])
  mean_logx <- mean(log(y[idx]))
  #thet <- uniroot(f1, c(0.00001, 10000))$root
  thet <- uniroot.all(f1, c(0.00001, 10000))[1]
  k <- mean_x / thet
  entr_g <- k + log(thet) +log(gamma(k)) + (1-k)*digamma(k)
  entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
  #print(sprintf("Gamma: %f, Shape: %f, Scale: %f", g, k, thet))
  return(list("Gamma"=g, "Shape"=k, "Scale"=thet, "MaxEntropy"=entr))
}

derHg <- 0
Hg1 <- 0
f2 <- function(x) {
  log((1-x)/x) - Hg1 - (x - 1)*derHg
}
derf2 <- function(x) {
  -1/((1-x)*x) - derHg
}

est_scem_tpgamma <- function(y, max.iter = 100, tol = 1e-6, g0 = 0.5, g1 = 0.4) {
  idx <- which(y>0)
  mean_x <- 0
  mean_logx <- 0
  entr_old <- 0
  entrs <- c()
  gammas <- c()
  Hgs <- c()
  
  f11 <- function(x) {
    t <- mean_logx - log(x) - digamma(mean_x / x)
    if(is.na(t)) {
      return(0)
    } else {
      return(t)
    }
  }
  f1 <- function(x) {
    mean_logx - log(x) - digamma(mean_x / x)
  }
  
  derf1 <- function(x) {
    -(1/x) + trigamma(mean_x / x) * (1/x^2)
  }
  
  for(i in 1:max.iter) {
    #print(i)
    gammas[i] <- g1
    mean_x <- mean(y)/(1-g0)
    mean_logx <- sum(log(y[idx]))/(length(y) * (1-g0))
    #thet <- newton.raphson(f1, 0.01, 100)$`root approximation`
    #thet <- secant(f1, 0.5, 0.51)$root
    #thet <- newtonRaphson(f1, 0.5, dfun = derf1)$root
    thet <- uniroot(f1, c(0.0001, 1000))$root
    #thet <- uniroot.all(f1, c(0.00001, 1))[1]
    k <- mean_x / thet
    Hg0 <- k + log(thet) + log(gamma(k)) + (1-k)*digamma(k)
    mean_x <- mean(y)/(1-g1)
    mean_logx <- sum(log(y[idx]))/(length(y) * (1-g1))
    #thet <- newton.raphson(f1, 0.01, 100)$`root approximation`
    #thet <- secant(f1, 0.5, 0.6)$root
    #thet <- newtonRaphson(f1, 0.5, dfun = derf1)$root
    thet <- uniroot(f1, c(0.0001, 1000))$root
    #thet <- uniroot.all(f1, c(0.00001, 1))[1]
    k <- mean_x / thet
    Hg1 <<- k + log(thet) + log(gamma(k)) + (1-k)*digamma(k)
    #print(sprintf("Iter: %d, Gamma: %f, Shape: %f, Scale: %f, Entropy: %f", i, g1, k, thet, entr))
    derHg <<- (Hg1 - Hg0)/(g1-g0)
    g0 <- g1
    #print(sprintf("Hg1: %s, derHg: %s", Hg1, derHg))
    #g1 <- secant(f2, 0.01, 0.9)$root
    #g1 <- newton.raphson(f2, 0.001, 0.99)$`root approximation`
    #g1 <- newtonRaphson(f2, 0.6, dfun = derf2)$root
    g1 <- uniroot(f2, c(0.001, 0.9999))$root
    #g1 <- uniroot.all(f2, c(0.000001, 0.9999))[1]
    #print(g1)
    #if(is.na(g1)) {
    #  g1 <- runif(1)
    #}
    mean_x <- mean(y)/(1-g1)
    mean_logx <- sum(log(y[idx]))/(length(y) * (1-g1))
    #thet <- newton.raphson(f1, 0.01, 100)$`root approximation`
    #thet <- secant(f1, 0.5, 0.51)$root
    #thet <- newtonRaphson(f1, 0.5, dfun = derf1)$root
    thet <- uniroot(f1, c(0.0001, 1000))$root
    #thet <- uniroot.all(f1, c(0.00001, 10))[1]
    k <- mean_x / thet
    Hg1 <<- k + log(thet) + log(gamma(k)) + (1-k)*digamma(k)
    Hg0 <- Hg1
    Hgs[i] <- Hg1
    entr <- -g1 * log(g1) - (1-g1)*log(1-g1) + (1-g1)*Hg1
    entrs[i] <- entr
    if(abs(entr-entr_old) < tol) {
      break
    } else {
      entr_old <- entr
    }
  }
  mean_x <- mean(y)/(1-g1)
  mean_logx <- sum(log(y[idx]))/(length(y) * (1-g1))
  #thet <- secant(f1, 0.5, 0.51)$root
  #thet <- newtonRaphson(f1, 0.5, dfun = derf1)$root
  thet <- uniroot(f1, c(0.0001, 1000))$root
  #thet <- uniroot.all(f1, c(0.00001, 10000))[1]
  k <- mean_x / thet
  Hg1 <<- k + log(thet) + log(gamma(k)) + (1-k)*digamma(k)
  entr <- -g1 * log(g1) - (1-g1)*log(1-g1) + (1-g1)*Hg1
  #print(sprintf("Gamma: %f, Shape: %f, Scale: %f", g1, k, thet))
  return(list("Gamma"=g1, "Shape"=k, "Scale"=thet, "Gammas" = gammas, "Hgs"=Hgs, 
              "MaxEntropy"=entr, "EntropyValues"=entrs, "IterUsed"=i))
}

n <- 1000
pzero <- 0.4  # Pr of observing a zero
shp <- 5 # k
scl <- 1.5 # theta
exp_val <- (1-pzero)*shp*scl
exp_logx <- (1-pzero)*(digamma(shp) + log(scl))
true_var <- (1-pzero)*(shp*scl^2 + pzero*shp^2*scl^2)
true_Hg <- shp + log(scl) + log(gamma(shp)) + (1-shp)*digamma(shp)
true_Hp <- -pzero*log(pzero)-(1-pzero)*log(1-pzero)+(1-pzero)*true_Hg

M <- 100
max.iter = 100
#gmmas_init <- seq(0.15, 0.9, 0.05)
gmmas_init <- c(0.6)
entrs1 <- matrix(0, ncol=1, nrow=M)
entrVals1 <- matrix(NA, ncol=max.iter, nrow=M)
hj11 <- matrix(0, ncol=1, nrow=M)
hj21 <- matrix(0, ncol=1, nrow=M)
vars1 <- matrix(0, ncol=1, nrow=M)
gammas1 <- matrix(0, ncol=1, nrow=M)
niters <- matrix(0, ncol=1, nrow=M)
entrs2 <- matrix(0, ncol=1, nrow=M)
entrVals2 <- matrix(NA, ncol=max.iter, nrow=M)
hj12 <- matrix(0, ncol=1, nrow=M)
hj22 <- matrix(0, ncol=1, nrow=M)
vars2 <- matrix(0, ncol=1, nrow=M)
gammas2 <- matrix(0, ncol=1, nrow=M)
entrs3 <- matrix(0, ncol=1, nrow=M)
hj13 <- matrix(0, ncol=1, nrow=M)
hj23 <- matrix(0, ncol=1, nrow=M)
vars3 <- matrix(0, ncol=1, nrow=M)
gammas3 <- matrix(0, ncol=1, nrow=M)
entrs4 <- matrix(0, ncol=1, nrow=M)
hj14 <- matrix(0, ncol=1, nrow=M)
hj24 <- matrix(0, ncol=1, nrow=M)
vars4 <- matrix(0, ncol=1, nrow=M)
gammas4 <- matrix(0, ncol=1, nrow=M)
iters_used1 <- 101
iters_used2 <- 101
i <- 1
while(i <= M) {
  print(i)
  y <- rsemicts(n, pzero = pzero, cts.density = "gamma", cts.param = list(shape=shp, scale=scl))
  est <- est_scem_tpgamma(y, max.iter = max.iter, g1 = 0.6, g0 = 0.61)
  entrs1[i] <- est$MaxEntropy
  entrVals1[i, (1:length(est$EntropyValues))] <- est$EntropyValues
  iters_used1 <- min(iters_used1, est$IterUsed)
  g1 <- est$Gamma
  shp1 <- est$Shape
  scl1 <- est$Scale
  gammas1[i] <- g1
  hj11[i] <- (1-g1)*shp1*scl1
  hj21[i] <- (1-g1)*(digamma(shp1) + log(scl1))
  vars1[i] <- (1-g1)*(shp1*scl1^2 + g1*shp1^2*scl1^2)
  est1 <- est_scem_tpgamma_old(y, max.iter = max.iter, g = gmmas_init[1])
  entrs2[i] <- est1$MaxEntropy
  iters_used2 <- min(iters_used2, est1$IterUsed)
  entrVals2[i, (1:length(est1$EntropyValues))] <- est1$EntropyValues
  g1 <- est1$Gamma
  shp1 <- est1$Shape
  scl1 <- est1$Scale
  gammas2[i] <- g1
  hj12[i] <- (1-g1)*shp1*scl1
  hj22[i] <- (1-g1)*(digamma(shp1) + log(scl1))
  vars2[i] <- (1-g1)*(shp1*scl1^2 + g1*shp1^2*scl1^2)
  # tp - em - gamma
  est1 <- est_em_tpgamma(y)
  entrs3[i] <- est1$MaxEntropy
  g1 <- est1$Gamma
  shp1 <- est1$Shape
  scl1 <- est1$Scale
  gammas3[i] <- g1
  hj13[i] <- (1-g1)*shp1*scl1
  hj23[i] <- (1-g1)*(digamma(shp1) + log(scl1))
  vars3[i] <- (1-g1)*(shp1*scl1^2 + g1*shp1^2*scl1^2)
  # tp - mle
  est1 <- est_mle_tpgamma(y)
  entrs4[i] <- est1$MaxEntropy
  g1 <- est1$Gamma
  shp1 <- est1$Shape
  scl1 <- est1$Scale
  gammas4[i] <- g1
  hj14[i] <- (1-g1)*shp1*scl1
  hj24[i] <- (1-g1)*(digamma(shp1) + log(scl1))
  vars4[i] <- (1-g1)*(shp1*scl1^2 + g1*shp1^2*scl1^2)
  i <- i + 1
}

print(sprintf("Entropy-new: %f; Entropy-old: %f; Entropy-tpem: %f; Entropy-mle: %f", 
              mean(entrs1, na.rm = TRUE), mean(entrs2), mean(entrs3), mean(entrs4)))
print(sprintf("MSE(hj1)-new: %f; MSE(hj1)-old: %f; MSE(hj1)-tpem: %f; MSE(hj1)-mle: %f", 
              mean(hj11-exp_val)^2, mean(hj12-exp_val)^2, mean(hj13-exp_val)^2, mean(hj14-exp_val)^2))
print(sprintf("MSE(hj2)-new: %f; MSE(hj2)-old: %f; MSE(hj2)-tpem: %f; MSE(hj2)-mle: %f", 
              mean(hj21-exp_logx)^2, mean(hj22-exp_logx)^2, mean(hj23-exp_logx)^2, mean(hj24-exp_logx)^2))
print(sprintf("Var-new: %f; Var-old: %f; Var-tpem: %f; Var-mle: %f", mean(vars1), mean(vars2), mean(vars3), mean(vars4)))
print(sprintf("Gamma-new: %f; Gamma-old: %f; Gamma-tpem: %f; Gamma-mle: %f", 
              mean(gammas1, na.rm = TRUE), mean(gammas2), mean(gammas3), mean(gammas4)))

entrVals1 <- entrVals1[,1:iters_used1]
entrVals2 <- entrVals2[,1:iters_used2]
entrVals1_arr <- rowMeans(entrVals1, na.rm = TRUE)
entrVals2_arr <- rowMeans(entrVals2, na.rm = TRUE)
conv_dta <- data.table("val" = entrVals1_arr[which(!is.na(entrVals1_arr))], "iter" = seq(1:length(entrVals1_arr[which(!is.na(entrVals1_arr))])), type = "new")
conv_dta <- rbindlist(list(conv_dta, data.table("val" = entrVals2_arr[which(!is.na(entrVals2_arr))], "iter" = seq(1:length(entrVals2_arr[which(!is.na(entrVals2_arr))])), type = "old")))
p <- ggplot(conv_dta, aes(iter, val)) + geom_line(aes(colour = type))
p <- p + ggtitle("Convergence of the new and the old algorithms")
p
if(0) {
dt1 <- data.table(sd = sqrt(vars1), method="new")
sd_dt <- data.table(sd = sqrt(vars2), method="old")
sd_dt <- rbindlist(list(dt1, sd_dt))
dt1 <- data.table(sdev = sqrt(vars1)-sqrt(vars2))
plt <- ggplot(data = dt1, aes(dt1$sdev.V1)) + 
         geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
         geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=1) +
         labs(x="value", y="Density") + ggtitle("main")

plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    #geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))
}
plot_multi_histogram(sd_dt, 'sd.V1', 'method')
}
