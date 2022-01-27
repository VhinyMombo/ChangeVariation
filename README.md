# ChangeVariation

---
title: "R Notebook"
output: html_notebook
---

---
title: "Binary Segmentation"
author: "Ennatiqi"
date: "1/25/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
data_generator <- function(n, positions, variance_values){
  library(gfpop)
  data = dataGenerator(n, positions, variance_values, type = 'variance')
  return(data)
}

library(sys)

positions = c(0.1,0.4,0.6, 0.8, 1)
variances = c(4, 1, 8, 5, 2)
n = 1000
y = data_generator(n, positions, variances)


css_stat <- function(y){
  stat = c()
  n = length(y)
  for(t in 1:n){
    a = sum(y[1:t]**2)/sum(y**2)
    b = t/n
    stat <- c(stat, sqrt(n/2)*abs(a-b))
  }
  return(stat)
}

binSeg <- function(s,e,y,pen=1.358){
  if (e-s==1) {
    return(c())
  } else {
    css_stats <- css_stat(y[s:e])
    cp <- which.max(css_stats) + s
    message("cp: ",cp)
    message("pen: ",pen)
    message("max(css_stats): ",max(css_stats))
    if (max(css_stats)> pen) {
      return(c(binSeg(s,cp,y,pen),cp,binSeg(cp+1,e,y,pen)))
    } else {
      return(c())
    }
  }
  
}

```

SIMULATION IN TIME
```{r}
N=10000
time_rcpp = c()
library(variance)
for(i in seq(from=2, to=N, by=80)){
  n = 500+i
  y = data_generator(n, positions, variances)
  start_time <- Sys.time()
    BS(0, length(y), y, 1.358)
  time_rcpp = c(time_rcpp, Sys.time()-start_time)
}
plot(time_rcpp, type = "s")
```

COMPARAISON AVEC LE PACKAGE CPT
```{r}
library(changepoint)
cpt = cpt.var(y,method="BinSeg",Q=5,class=FALSE) 
my_cpt = unlist(BS(0,length(y), y, 1.358))

```



```{r}
cost_function <- function(y){
  n = length(y)
  #log(sum((y - mean(y))**2)/n)
  #var(y)
  m = 0
  if (n==0) {
    C = 0
  }
  else{
    C=n * (log(2*pi) + log(mean(((y-m)**2)))+ 1)}
  return(C)
}

penalty_function <- function(n, params){
  beta = params  *  log(n)
  return(beta)
}

penalty_function2 <- function(data){
  n = length(data)
  return (2*var(data)*log(n))
}

minimisation <- function(data, F_comp,beta){
  t_star = length(data)
  #print(t_star)
  F_temp = cost_function(data)
  #print(F_temp)
  i_temp = 0
  for (i in (2:(t_star))) {
    #print(i)
    F_tau = F_comp[i-1] + cost_function(data[i:t_star])  + beta

    if (F_tau < F_temp) {
      F_temp = F_tau
      i_temp = i-1}
    #print(F_temp)
  }
  return(list(F_tau_star = F_temp, tau_prime = i_temp))
}

OP <- function(data,params=1,K=0){
  #"if cost function to be minus the log-likelihood then the constant K = 0
  #if we take -penalised log-likelihood then K =  penalisation factor."
  n = length(data)
  beta = penalty_function(n,params)
  #beta = penalty_function2(data)
  #print(beta)
  cp <- rep(0,n)
  F_comp <- rep(0,n)
  F_comp[1] = -beta
  #cp = c(Inf)
  for (tau in c(2:n)) {
    #print(tau)
    res_minimisation <- minimisation(data = data[1:tau],
                                     beta = beta,
                                     F_comp = F_comp)
    #print(res_minimisation)
    F_comp[tau] <- res_minimisation$F_tau_star
    cp[tau] <- res_minimisation$tau_prime
    #print(F_comp)
    #print(cp)
  }

  v <- cp[n]
  P <- cp[n]

  while (v > 1)
  {
    P <- c(P, cp[v])
    v <- cp[v]
  }
  P <- rev(P)[-1]

  return(unique(P))
}


minimisation_PELT <- function(data, t_star ,F_comp,beta, R ,K = 0){
  t_star = length(data)
  #print(t_star)
  F_temp = cost_function(data)
  #print(F_temp)
  F_tau_K = c()
  i_temp = 0
  for (i in c(R)) {
    #print(i)
    F_tau = F_comp[i-1] + cost_function(data[i:t_star])  + beta

    F_tau_K = c(F_tau_K,F_comp[i-1] + cost_function(data[i:t_star])  + K)

    if (F_tau < F_temp) {
      F_temp = F_tau
      i_temp = i-1}
    #print(F_temp)
  }
  R_new <- R[which(F_tau_K < F_temp)]

  return(list(F_tau_star = F_temp, tau_prime = i_temp, R_new = R_new))
}


PELT <- function(data,params=1,K=0){
  #"if cost function to be minus the log-likelihood then the constant K = 0
  #if we take -penalised log-likelihood then K =  penalisation factor."
  n = length(data)
  beta = penalty_function(n,params)
  #print(beta)
  cp <- rep(0,n)
  F_comp <- rep(0,n)
  F_comp[1] = -beta
  R = c(2)
  #cp = c(Inf)
  for (tau in c(2:n)) {
    #print(tau)
    res_minimisation <- minimisation_PELT(data = data[1:tau],
                                          beta = beta,
                                          F_comp = F_comp,
                                          R = R,
                                          K = 0)
    #print(res_minimisation)
    F_comp[tau] <- res_minimisation$F_tau_star
    cp[tau] <- res_minimisation$tau_prime
    R <- c(res_minimisation$R_new, tau)
    #print(F_comp)
    #print(cp)
  }
  v <- cp[n]
  P <- cp[n]

  while (v > 1)
  {
    P <- c(P, cp[v])
    v <- cp[v]
  }
  P <- rev(P)[-1]

  return(unique(P))
}
```
COMPARAISON ENTRE ALGO
```{r}
library(variance)
positions = c(0.1,0.4,0.6, 0.8, 1)
variances = c(4, 1, 8, 5, 2)

first_true = c()
first_BS = c()
first_OP = c()
first_PELT = c()
N = 1000
for(i in seq(from = 100, to = N, by = 100)){
  y = data_generator(i, positions, variances)
  first_BS = c(first_BS, unlist(BS(0, i, y, 1.358))[2])
  first_OP = c(first_OP, OP(y, params = 2)[2])
  first_PELT = c(first_PELT, OP(y, params = 2)[2])
  first_true = c(first_true, i*positions[2])
}
```

```{r}
plot(first_true,col="green")
par(new=TRUE)

plot(first_BS,col="orange")

```

