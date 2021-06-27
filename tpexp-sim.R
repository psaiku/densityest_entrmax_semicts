library(semicts)
library(data.table)
library(ggplot2)
library(reshape2)

source("semicts_est_utils.R")

est_scem_tpexp_old <- function(y, max.iter = 100, tol = 1e-6, g = 0.8) {
  entr_old <- 0
  entrs <- c()
  for(i in 1:max.iter) {
    mn <- mean(y)/(1-g)
    entr_g <- 1 - log(1/mn)
    g <- 1/(1+exp(entr_g))
    entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
    entrs[i] <- entr
    #print(sprintf("Iter: %d, Gamma: %f, Lambda: %f, Entropy: %f", i, g, (1/mn), entr))
    
    if(abs(entr-entr_old) < tol) {
      break
    } else {
      entr_old <- entr
    }
  }
  return(list("Gamma"=g, "Lambda"=(1/mn), "MaxEntropy"=entr, "EntropyValues"=entrs, "IterUsed"=i))
}

est_em_tpexp <- function(y) {
  idx <- which(y>0)
  g <- sum(y==0)/length(y)
  mn <- mean(y[idx])
  entr_g <- 1 - log(1/mn)
  entr <- -g * log(g) - (1-g)*log(1-g) + (1-g)*entr_g
  return(list("Gamma"=g, "Lambda"=(1/mn), "MaxEntropy"=entr))
}

derHg <- 0
Hg1 <- 0
f1 <- function(x) {
  log((1-x)/x) - Hg1 - (x - 1)*derHg
}

est_scem_tpexp <- function(y, max.iter = 100, tol = 1e-6, g0 = 0.8, g1 = 0.4, Hg0 = 0) {
  entr_old <- 0
  entrs <- c()
  for(i in 1:max.iter) {
    #mn <- mean(y)/(1-g1)
    mn <- exp_val/(1-g1)
    Hg1 <<- 1 - log(1/mn)
    Hg0 <- 1 - log((1-g0)/exp_val)
    #print(sprintf("Iter: %d, Gamma: %f, Lambda: %f, Entropy: %f", i, g1, (1/mn), entr))
    derHg <<- (Hg1 - Hg0)/(g1-g0)
    g0 <- g1
    g1 <- uniroot(f1, c(0.0001, 0.9999))$root
    entr <- -g1 * log(g1) - (1-g1)*log(1-g1) + (1-g1)*Hg1
    entrs[i] <- entr
    Hg0 <- Hg1
    if(abs(entr-entr_old) < tol) {
      break
    } else {
      entr_old <- entr
    }
  }
  lambda <- mean(y)/(1-g1)
  Hg1 <- 1 - log(1/lambda)
  entr <- -g1 * log(g1) - (1-g1)*log(1-g1) + (1-g1)*Hg1
  return(list("Gamma"=g1, "Lambda"=1/lambda, "MaxEntropy"=entr, "EntropyValues"=entrs, "IterUsed"=i))
}

n <- 1000
pzero <- 0.1  # Pr of observing a zero
lmbda <- 0.5
exp_val <- (1-pzero)*(1/lmbda)
#true_Hp <- -pzero*log(pzero)-(1-pzero)*log(1-pzero)+(1-pzero)*(1-log(lmbda))
#pz <- ((1-pzero) - sqrt((1-pzero)^2 - 4*lmbda^2))/(2*lmbda)
pz <- ((1-pzero+2*lmbda) - sqrt((1-pzero+2*lmbda)^2 - 4*lmbda^2))/(2*lmbda)
lmbda1 <- (1-pz)*lmbda/(1-pzero)
true_Hp <- -pz*log(pz)-(1-pz)*log(1-pz)+(1-pz)*(1-log(lmbda1))
true_var <- (1-pz^2)*(1/lmbda1^2)
M <- 100
max.iter = 80
#gmmas_init <- seq(0.15, 0.9, 0.05)
gmmas_init <- c(0.7)
entrs1 <- matrix(0, ncol=1, nrow=M)
entrVals1 <- matrix(NA, ncol=max.iter, nrow=M)
gammas1 <- matrix(0, ncol=1, nrow=M)
hj11 <- matrix(0, ncol=1, nrow=M)
vars1 <- matrix(0, ncol=1, nrow=M)
niters <- matrix(0, ncol=1, nrow=M)
entrs2 <- matrix(0, ncol=1, nrow=M)
entrVals2 <- matrix(NA, ncol=max.iter, nrow=M)
gammas2 <- matrix(0, ncol=1, nrow=M)
hj12 <- matrix(0, ncol=1, nrow=M)
vars2 <- matrix(0, ncol=1, nrow=M)
entrs3 <- matrix(0, ncol=1, nrow=M)
gammas3 <- matrix(0, ncol=1, nrow=M)
hj13 <- matrix(0, ncol=1, nrow=M)
vars3 <- matrix(0, ncol=1, nrow=M)
entrs4 <- matrix(0, ncol=1, nrow=M)
gammas4 <- matrix(0, ncol=1, nrow=M)
hj14 <- matrix(0, ncol=1, nrow=M)
vars4 <- matrix(0, ncol=1, nrow=M)
entrs5 <- matrix(0, ncol=1, nrow=M)
gammas5 <- matrix(0, ncol=1, nrow=M)
hj15 <- matrix(0, ncol=1, nrow=M)
vars5 <- matrix(0, ncol=1, nrow=M)
entrs6 <- matrix(0, ncol=1, nrow=M)
gammas6 <- matrix(0, ncol=1, nrow=M)
hj16 <- matrix(0, ncol=1, nrow=M)
vars6 <- matrix(0, ncol=1, nrow=M)
iters_used1 <- 101
iters_used2 <- 101
i <- 1
while(i <= M) {
  print(i)
  y <- rsemicts(n, pzero = pzero, cts.density = "exp", cts.param = list(rate=lmbda))
  est <- est_scem_tpexp(y, max.iter = max.iter, g0 = 0.3, g1 = gmmas_init[1], Hg0 = 1)
  entrs1[i] <- est$MaxEntropy
  entrVals1[i, (1:length(est$EntropyValues))] <- est$EntropyValues
  iters_used1 <- min(iters_used1, est$IterUsed)
  g1 <- est$Gamma
  gammas1[i] <- g1
  lambda1 <- est$Lambda
  hj11[i] <- (1-g1)*(1/lambda1)
  vars1[i] <- (1-g1^2)*(1/lambda1^2)
  est <- est_scem_tpexp_old(y, max.iter = max.iter, g = gmmas_init[1])
  entrs2[i] <- est$MaxEntropy
  iters_used2 <- min(iters_used2, est$IterUsed)
  entrVals2[i, (1:length(est$EntropyValues))] <- est$EntropyValues
  g1 <- est$Gamma
  gammas2[i] <- g1
  lambda1 <- est$Lambda
  hj12[i] <- (1-g1)*(1/lambda1)
  vars2[i] <- (1-g1^2)*(1/lambda1^2)
  est <- est_em_tpexp(y)
  entrs3[i] <- est$MaxEntropy
  g1 <- est$Gamma
  gammas3[i] <- g1
  lambda1 <- est$Lambda
  hj13[i] <- (1-g1)*(1/lambda1)
  vars3[i] <- (1-g1^2)*(1/lambda1^2)
  est <- est_mle_tpexp(y)
  entrs4[i] <- est$MaxEntropy
  g1 <- est$Gamma
  gammas4[i] <- g1
  lambda1 <- est$Lambda
  hj14[i] <- (1-g1)*(1/lambda1)
  vars4[i] <- (1-g1^2)*(1/lambda1^2)
  est <- est_mle_tpgamma(y)
  entrs5[i] <- est$MaxEntropy
  g1 <- est$Gamma
  gammas5[i] <- g1
  shp <- est$Shape
  scl <- est$Scale
  hj15[i] <- (1-g1)*shp*scl
  vars5[i] <- (1-g1)*(shp*scl^2 + g1*shp^2*scl^2)
  est <- est_mle_tplnorm(y)
  entrs6[i] <- est$MaxEntropy
  g1 <- est$Gamma
  gammas6[i] <- g1
  meanlog <- est$mean
  sdlog <- est$sd
  hj16[i] <- (1-g1)*exp(meanlog + sdlog^2/2)
  vars6[i] <- (1-g1)*(exp(2*meanlog + 2*sdlog^2) - (1-g1)*exp(2*meanlog + sdlog^2))
  i <- i + 1
}

print(sprintf("MSE-Entropy-new: %f; Entropy-old: %f; Entropy-tpem: %f; Entropy-mle(exp): %f; Entropy-mle(gamma): %f; Entropy-mle(lnorm): %f", 
              mean((entrs1-true_Hp)/true_Hp), mean((entrs2-true_Hp)/true_Hp), mean((entrs3-true_Hp)/true_Hp), mean(entrs4-true_Hp)^2, mean(entrs5-true_Hp)^2, mean(entrs6-true_Hp)^2))
print(sprintf("MSE(hj1)-new: %f; MSE(hj1)-old: %f; MSE(hj1)-tpem: %f; MSE(hj1)-mle(exp): %f; MSE(hj1)-mle(gamma): %f; MSE(hj1)-mle(lnorm): %f", 
              mean(hj11-exp_val)^2, mean(hj12-exp_val)^2, mean(hj13-exp_val)^2, mean(hj14-exp_val)^2, mean(hj15-exp_val)^2, mean(hj16-exp_val)^2))
print(sprintf("MSE-Var-new: %f; Var-old: %f; Var-tpem: %f; Var-mle(exp): %f; Var-mle(gamma): %f; Var-mle(lnorm): %f", 
              mean((vars1-true_var)/true_var), mean((vars2-true_var)/true_var), mean((vars3-true_var)/true_var), mean((vars4-true_var)^2), mean((vars5-true_var)^2), mean((vars6-true_var)^2)))
print(sprintf("Gamma-new: %f; Gamma-old: %f; Gamma-tpem: %f; Gamma-mle(exp): %f; Gamma-mle(gamma): %f; Gamma-mle(lnorm): %f", 
              mean(gammas1), mean(gammas2), mean(gammas3), mean(gammas4), mean(gammas5), mean(gammas6)))

entr_vals1 <- colMeans(entrVals1)
entr_vals2 <- colMeans(entrVals2)
tp_em_entr_val <- mean(entrs3)
n_iter <- max(sum(!is.na(entr_vals1)), sum(!is.na(entr_vals2)))
entr_vals1 <- entr_vals1[1:n_iter]
entr_vals2 <- entr_vals2[1:n_iter]
plot_dt <- data.table(iters = seq(1:n_iter), AEM = entr_vals1, AEM_p = entr_vals2, 
                      true_Hp = true_Hp, tp_EM = tp_em_entr_val)
plot_dt <- melt(plot_dt, id.vars = 'iters', variable.names = 'model')
p <- ggplot(plot_dt, aes(iters, value)) +  geom_line(aes(linetype = variable, colour = variable), size=1.25)
p <- p + scale_linetype_manual(values=c("solid", "twodash", "longdash", "dotdash"))
p <- p + geom_point(aes(shape = variable), size = 1.5)
p <- p + ggtitle(bquote("Convergence for " ~ gamma == ~.(pzero) ~ " and " ~ lambda == ~. (lmbda)))
p <- p + xlab("Number of Iterations") + ylab("Entropy")
if(0) {
entrVals1 <- entrVals1[1:iters_used1,]
entrVals2 <- entrVals2[1:iters_used2,]
entrVals1_arr <- rowMeans(entrVals1, na.rm = TRUE)
entrVals2_arr <- rowMeans(entrVals2, na.rm = TRUE)
conv_dta <- data.table("val" = entrVals1_arr[which(!is.na(entrVals1_arr))], "iter" = seq(1:length(entrVals1_arr[which(!is.na(entrVals1_arr))])), type = "new")
conv_dta <- rbindlist(list(conv_dta, data.table("val" = entrVals2_arr[which(!is.na(entrVals2_arr))], "iter" = seq(1:length(entrVals2_arr[which(!is.na(entrVals2_arr))])), type = "old")))
p <- ggplot(conv_dta, aes(iter, val)) + geom_line(aes(colour = type))
p <- p + ggtitle("Convergence of the new and the old algorithms")
p
##########################3
dt1 <- data.table(sdev = sqrt(vars1)-sqrt(vars2))
plt <- ggplot(data = dt1, aes(dt1$sdev.V1)) + 
  geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
  geom_vline(aes(xintercept=0), color="black", linetype="dashed", size=1) +
  labs(x="value", y="Density") + ggtitle("main")
plt
dt1 <- data.table(sd = sqrt(vars1), method="new")
sd_dt <- data.table(sd = sqrt(vars2), method="old")
sd_dt <- rbindlist(list(dt1, sd_dt))

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
