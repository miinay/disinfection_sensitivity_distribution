############################################################
###Procedure
#0. Package 
#1. Data set
#2. Stan models
#3. Summary of model fits 
#4. Preparation for visualization
#5. Visualization 
############################################################

# R.version 
# _                                
# platform       x86_64-w64-mingw32               
# arch           x86_64                           
# os             mingw32                          
# crt            ucrt                             
# system         x86_64, mingw32                  
# status                                          
# major          4                                
# minor          4.0                              
# year           2024                             
# month          04                               
# day            24                               
# svn rev        86474                            
# language       R                                
# version.string R version 4.4.0 (2024-04-24 ucrt)
# nickname       Puppy Cup

#### 0. Package ----
### packages
library(openxlsx)
library(tidyverse)
library(ssdtools)
library(ggplot2)
library(dplyr)
library(ggbeeswarm)
library(patchwork)
library(rstan)
library(loo)
library(cowplot)

### rstan option
rstan_options(auto_write=TRUE)

#### 1. Data set ----
### Import rate constant data
## UV
raw_data_uv <- read.xlsx("Summary_k_UV.xlsx")
raw_data_uv[is.na(raw_data_uv)] <- 0 
raw_data_uv <- raw_data_uv %>%
  mutate("method" = "UV")

## Ozone
raw_data_ozone <- read.xlsx("Summary_k_Ozone.xlsx")
raw_data_ozone[is.na(raw_data_ozone)] <- 0 
raw_data_ozone <- raw_data_ozone %>%
  mutate("method" = "Ozone")
raw_data_ozone <- raw_data_ozone[,-c(5:6)] # remove min and max

## Free chlorine
raw_data_fc <- read.xlsx("Summary_k_FreeChlorine.xlsx")
raw_data_fc[is.na(raw_data_fc)] <- 0 
raw_data_fc <- raw_data_fc %>%
  mutate("method" = "FC")

### Prepare guideline values
## UV
gl_UV <- 46.5 # Required UV dose for n log inactivation = 46.5n (mJ cm^-2)
## Ozone
gl_OZ <- 0.113 # Required ozone CT for n log inactivation = 0.113*n
## Free chlorine
gl_FC <- 0.75 # Required free chlorine CT for n log inactivation = 0.75n

#### 2. Stan models -----
### 2-1. UV ----
## UV data
input_data_UV <- list(N=length(raw_data_uv$`Species`), 
                   tMean= raw_data_uv$`mean_k`,
                   tSD= 2*raw_data_uv$`sd_k`
)

## Model fitting for 2-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 2 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_UV, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_UV_2 <- models
fit_UV_2 <- fit
pos_UV_2 <- pos

#save(models_UV_2, fit_UV_2, pos_UV_2, file = "Stan_UV_n2.RData")

## Model fitting for 4-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 4 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k]
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_UV, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_UV_4 <- models
fit_UV_4 <- fit
pos_UV_4 <- pos

#save(models_UV_4, fit_UV_4, pos_UV_4, file = "Stan_UV_n4.RData")

## Model fitting for 6-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 6 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_UV, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_UV_6 <- models
fit_UV_6 <- fit
pos_UV_6 <- pos

#save(models_UV_6, fit_UV_6, pos_UV_6, file = "Stan_UV_n6.RData")

### 2-2. Ozone----
## Ozone data
input_data_OZ <- list(N=length(raw_data_ozone$`Species`), 
                      tMean= raw_data_ozone$`mean_k`,
                      tSD= 2*raw_data_ozone$`sd_k`
)

## Model fitting for 2-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 2 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_OZ, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_OZ_2 <- models
fit_OZ_2 <- fit
pos_OZ_2 <- pos

#save(models_OZ_2, fit_OZ_2, pos_OZ_2, file = "Stan_OZ_n2.RData")

## Model fitting for 4-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 4 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_OZ, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_OZ_4 <- models
fit_OZ_4 <- fit
pos_OZ_4 <- pos

#save(models_OZ_4, fit_OZ_4, pos_OZ_4, file = "Stan_OZ_n4.RData")

## Model fitting for 6-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 6 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_OZ, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_OZ_6 <- models
fit_OZ_6 <- fit
pos_OZ_6 <- pos

#save(models_OZ_6, fit_OZ_6, pos_OZ_6, file = "Stan_OZ_n6.RData")

### 2-3. Free Chlorine (FC) ----
## FC data
input_data_FC <- list(N=length(raw_data_fc$`Species`), 
                      tMean= raw_data_fc$`mean_k`,
                      tSD= 2*raw_data_fc$`sd_k`
)

## Model fitting for 2-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 2 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_FC, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_FC_2 <- models
fit_FC_2 <- fit
pos_FC_2 <- pos

#save(models_FC_2, fit_FC_2, pos_FC_2, file = "Stan_FC_n2.RData")

## Model fitting for 4-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 4 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k] 
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_FC, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_FC_4 <- models
fit_FC_4 <- fit
pos_FC_4 <- pos

#save(models_FC_4, fit_FC_4, pos_FC_4, file = "Stan_FC_n4.RData")

## Model fitting for 6-log reduction
distributions <- c("lognormal", "gamma", "weibull") 
code <- sprintf("
  data{
    int<lower=1> N; // number of data
    vector[N] tMean;
    vector[N] tSD;
  }
  parameters{
    real par[2]; // 2-para model
    vector<lower=-1, upper=1>[N] uE;	// Uniform value for sampling between lower and upper 95%% CI
  }
  transformed parameters{
    vector[N] tk; 	//  transformed parameter (the corresponding k value for each n-log reduction)
    tk = 6 * log(10) ./ (tMean + uE .* tSD); // [n-log reduction] * log(10) / [k]
  }
  model{
    // Contribution to likelihood of SSD
    target += %s_lpdf(tk  | par[1], par[2]);
  }
  generated quantities {
    // likelihood for calculation of looIC
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = %s_lpdf(tk[i]  | par[1], par[2]);
    }
  }
", distributions, distributions)
names(code) <- distributions

models <- mapply(stan_model, model_code=code)
seed_value <- 12345
fit <- mapply(function(model) {
  sampling(model, data = input_data_FC, iter = 10000, warmup = 3000, chains = 4, seed = seed_value)
}, models)
pos <- mapply(function(z) rstan::extract(z)$par, fit, SIMPLIFY=FALSE)

models_FC_6 <- models
fit_FC_6 <- fit
pos_FC_6 <- pos

#save(models_FC_6, fit_FC_6, pos_FC_6, file = "Stan_FC_n6.RData")

#### 3. Summary of model fits  -----
### 3-1. UV ----
# Load RData if necessary
# load('Stan_UV_n2.RData')
# load('Stan_UV_n4.RData')
# load('Stan_UV_n6.RData')

# Make lists to store the results
means <- list()
sds <- list()
res1 <- list()
res2 <- list()
ll <- list()
waic <- list()
cens_w <- list()
cens_g <- list()
cens_ln <- list()
percentile <- list()
para_list <- list()

for (i in c(2,4,6)){
  # get posterior samples
  pos_name <- paste0("pos_UV_",i)
  fit_name <- paste0("fit_UV_",i)
  pos <- get(pos_name)
  fit <- get(fit_name)
  
  # calculate mean and sd
  means[[i]] <- data.frame(weibull = pos$weibull[,2]*gamma(1+1/pos$weibull[,1]),
                              gamma = pos$gamma[,1] / pos$gamma[,2],
                              lognormal = exp(pos$lognormal[,1]+pos$lognormal[,2]^2/2)) 
  sds[[i]] <- data.frame(weibull = pos$weibull[,2]*sqrt((gamma(1 + 2/pos$weibull[,1]) - (gamma(1 + 1/pos$weibull[,1]))^2)),
                            gamma = sqrt(pos$gamma[,1]) / pos$gamma[,2],
                            lognormal = sqrt((exp(pos$lognormal[,2]^2) - 1) * exp(2 * pos$lognormal[,1] + pos$lognormal[,2]^2))) 
  
  # calculate credible intervals for mean and SD, log-likelihood (ll), WAIC
  a_percentile <- c(0.025, 0.5, 0.975)
  res1[[i]] <- data.frame(mean = apply(means[[i]], 2, quantile, a_percentile)) %>%
    mutate(nlog = rep(i, n()))
  res2[[i]] <- data.frame(sd = apply(sds[[i]], 2, quantile, a_percentile)) %>%
    mutate(nlog = rep(i, n()))
  ll[[i]] <- data.frame(ll = mapply(function(z) loo(extract_log_lik(z))$looic, fit)) %>%
    mutate(nlog = rep(i, n()))
  waic[[i]] <- data.frame(waic=mapply(function(z) waic(extract_log_lik(z))$waic, fit) ) %>%
    mutate(nlog = rep(i, n()))
  
  # calculate percentile values for estimated distributions
  cens_w[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                           function(p) quantile(qweibull(p = p, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_g[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                           function(p) quantile(qgamma(p = p, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_ln[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                            function(p) quantile(qlnorm(p = p, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_w[[i]] <- cbind(i,"w", round(cens_w[[i]], 1))
  cens_g[[i]] <- cbind(i, "g", round(cens_g[[i]], 1))
  cens_ln[[i]] <- cbind(i, "ln", round(cens_ln[[i]], 1))
  percentile[[i]] <- rbind(cens_w[[i]], cens_g[[i]], cens_ln[[i]])

  # obtain parameters of each estimated distribution
  calculate_stats <- function(parameters) {
    median_val <- apply(parameters, 2, median)
    sd_val <- apply(parameters, 2, function(x) sd(x))
    return(data.frame(median = median_val, sd = sd_val))
  }
  weibull_stats <- calculate_stats(pos$weibull)
  gamma_stats <- calculate_stats(pos$gamma)
  lognormal_stats <- calculate_stats(pos$lognormal)
  para_list[[i]] <- data.frame(
    weibull = cbind(i, "weibull", weibull_stats),
    gamma = cbind(i, "gamma", gamma_stats),
    lognormal = cbind(i, "lognormal", lognormal_stats)
  )
 }

# store results
percentile_UV <- do.call(rbind, percentile)
colnames(percentile_UV) <- c("i", "dist", 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
res1_UV <- do.call(rbind, res1)
res2_UV <- do.call(rbind, res2)
ll_UV <- do.call(rbind, ll)
waic_UV <- do.call(rbind, waic)
MeanSD_UV <- round(cbind(res1_UV, res2_UV),1)
aicc_UV <- round(cbind(ll_UV,waic_UV),1)
para_list_UV <- do.call(rbind, para_list)

# export
write.csv(percentile_UV, "percentile_UV.csv") 
write.csv(MeanSD_UV, "MeanSD_UV.csv")  # for Table S1
write.csv(aicc_UV, "aicc_UV.csv")   # for Table S1
write.csv(para_list_UV, "parameters_UV.csv")

### 3-2. Ozone ----
# Load RData if necessary
# load('Stan_OZ_n2.RData')
# load('Stan_OZ_n4.RData')
# load('Stan_OZ_n6.RData')

# Make lists to store the results
means <- list()
sds <- list()
res1 <- list()
res2 <- list()
ll <- list()
waic <- list()
cens_w <- list()
cens_g <- list()
cens_ln <- list()
percentile <- list()
para_list <- list()

for (i in c(2, 4, 6)){
  # get posterior samples
  pos_name <- paste0("pos_OZ_",i)
  fit_name <- paste0("fit_OZ_",i)
  pos <- get(pos_name)
  fit <- get(fit_name)
  
  # calculate mean and sd
  means[[i]] <- data.frame(weibull = pos$weibull[,2]*gamma(1+1/pos$weibull[,1]),
                           gamma = pos$gamma[,1] / pos$gamma[,2],
                           lognormal = exp(pos$lognormal[,1]+pos$lognormal[,2]^2/2)) 
  sds[[i]] <- data.frame(weibull = pos$weibull[,2]*sqrt((gamma(1 + 2/pos$weibull[,1]) - (gamma(1 + 1/pos$weibull[,1]))^2)),
                         gamma = sqrt(pos$gamma[,1]) / pos$gamma[,2],
                         lognormal = sqrt((exp(pos$lognormal[,2]^2) - 1) * exp(2 * pos$lognormal[,1] + pos$lognormal[,2]^2))) 
  
  # calculate credible intervals for mean and SD, log-likelihood (ll), WAIC
  a_percentile <- c(0.025, 0.5, 0.975)
  res1[[i]] <- data.frame(mean = apply(means[[i]], 2, quantile, a_percentile)) %>%
    mutate(nlog = rep(i, n()))
  res2[[i]] <- data.frame(sd = apply(sds[[i]], 2, quantile, a_percentile)) %>%
    mutate(nlog = rep(i, n()))
  ll[[i]] <- data.frame(ll = mapply(function(z) loo(extract_log_lik(z))$looic, fit)) %>%
    mutate(nlog = rep(i, n()))
  waic[[i]] <- data.frame(waic=mapply(function(z) waic(extract_log_lik(z))$waic, fit) ) %>%
    mutate(nlog = rep(i, n()))
  
  # calculate percentile values for estimated distributions
  cens_w[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                        function(p) quantile(qweibull(p = p, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_g[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                        function(p) quantile(qgamma(p = p, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_ln[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                         function(p) quantile(qlnorm(p = p, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_w[[i]] <- cbind(i,"w", round(cens_w[[i]], 1))
  cens_g[[i]] <- cbind(i, "g", round(cens_g[[i]], 1))
  cens_ln[[i]] <- cbind(i, "ln", round(cens_ln[[i]], 1))
  percentile[[i]] <- rbind(cens_w[[i]], cens_g[[i]], cens_ln[[i]])
  
  # obtain parameters of each estimated distribution
  calculate_stats <- function(parameters) {
    median_val <- apply(parameters, 2, median)
    sd_val <- apply(parameters, 2, function(x) sd(x))
    return(data.frame(median = median_val, sd = sd_val))
  }
  weibull_stats <- calculate_stats(pos$weibull)
  gamma_stats <- calculate_stats(pos$gamma)
  lognormal_stats <- calculate_stats(pos$lognormal)
  para_list[[i]] <- data.frame(
    weibull = cbind(i, "weibull", weibull_stats),
    gamma = cbind(i, "gamma", gamma_stats),
    lognormal = cbind(i, "lognormal", lognormal_stats)
  )
}

# store results
percentile_OZ <- do.call(rbind, percentile)
colnames(percentile_OZ) <- c("i", "dist", 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
res1_OZ <- do.call(rbind, res1)
res2_OZ <- do.call(rbind, res2)
ll_OZ <- do.call(rbind, ll)
waic_OZ <- do.call(rbind, waic)
MeanSD_OZ <- round(cbind(res1_OZ, res2_OZ),4)
aicc_OZ <- round(cbind(ll_OZ,waic_OZ),1)
para_list_OZ <- do.call(rbind, para_list)

# export
write.csv(percentile_OZ, "percentile_OZ.csv")
write.csv(MeanSD_OZ, "MeanSD_OZ.csv")  # for Table S1
write.csv(aicc_OZ, "aicc_OZ.csv")  # for Table S1
write.csv(para_list_OZ, "parameters_OZ.csv")

### 3-3. Free chlorine ----
# Load RData if necessary
# load('Stan_FC_n2.RData')
# load('Stan_FC_n4.RData')
# load('Stan_FC_n6.RData')

# Make lists to store the results
means <- list()
sds <- list()
res1 <- list()
res2 <- list()
ll <- list()
waic <- list()
cens_w <- list()
cens_g <- list()
cens_ln <- list()
percentile <- list()
para_list <- list()

for (i in c(2, 4, 6)){
  # get posterior samples
  pos_name <- paste0("pos_FC_",i)
  fit_name <- paste0("fit_FC_",i)
  pos <- get(pos_name)
  fit <- get(fit_name)
  
  # calculate mean and sd
  means[[i]] <- data.frame(weibull = pos$weibull[,2]*gamma(1+1/pos$weibull[,1]),
                           gamma = pos$gamma[,1] / pos$gamma[,2],
                           lognormal = exp(pos$lognormal[,1]+pos$lognormal[,2]^2/2)) 
  sds[[i]] <- data.frame(weibull = pos$weibull[,2]*sqrt((gamma(1 + 2/pos$weibull[,1]) - (gamma(1 + 1/pos$weibull[,1]))^2)),
                         gamma = sqrt(pos$gamma[,1]) / pos$gamma[,2],
                         lognormal = sqrt((exp(pos$lognormal[,2]^2) - 1) * exp(2 * pos$lognormal[,1] + pos$lognormal[,2]^2))) 
  a_percentile <- c(0.025, 0.5, 0.975)
  res1[[i]] <- data.frame(mean = apply(means[[i]], 2, quantile, a_percentile)) %>%
    mutate(nlog = rep(i, n()))
  res2[[i]] <- data.frame(sd = apply(sds[[i]], 2, quantile, a_percentile)) %>%
    mutate(nlog = rep(i, n()))
  ll[[i]] <- data.frame(ll = mapply(function(z) loo(extract_log_lik(z))$looic, fit)) %>%
    mutate(nlog = rep(i, n()))
  waic[[i]] <- data.frame(waic=mapply(function(z) waic(extract_log_lik(z))$waic, fit) ) %>%
    mutate(nlog = rep(i, n()))
  
  # calculate credible intervals for mean and SD, log-likelihood (ll), WAIC
  cens_w[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                        function(p) quantile(qweibull(p = p, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_g[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                        function(p) quantile(qgamma(p = p, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_ln[[i]] <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                         function(p) quantile(qlnorm(p = p, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025, 0.5, 0.975)))
  cens_w[[i]] <- cbind(i,"w", round(cens_w[[i]], 1))
  cens_g[[i]] <- cbind(i, "g", round(cens_g[[i]], 1))
  cens_ln[[i]] <- cbind(i, "ln", round(cens_ln[[i]], 1))
  percentile[[i]] <- rbind(cens_w[[i]], cens_g[[i]], cens_ln[[i]])
  
  # obtain parameters of each estimated distribution
  calculate_stats <- function(parameters) {
    median_val <- apply(parameters, 2, median)
    sd_val <- apply(parameters, 2, function(x) sd(x))
    return(data.frame(median = median_val, sd = sd_val))
  }
  weibull_stats <- calculate_stats(pos$weibull)
  gamma_stats <- calculate_stats(pos$gamma)
  lognormal_stats <- calculate_stats(pos$lognormal)
  para_list[[i]] <- data.frame(
    weibull = cbind(i, "weibull", weibull_stats),
    gamma = cbind(i, "gamma", gamma_stats),
    lognormal = cbind(i, "lognormal", lognormal_stats)
  )
}

# store results
percentile_FC <- do.call(rbind, percentile)
colnames(percentile_FC) <- c("i", "dist", 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
res1_FC <- do.call(rbind, res1)
res2_FC <- do.call(rbind, res2)
ll_FC <- do.call(rbind, ll)
waic_FC <- do.call(rbind, waic)
MeanSD_FC <- round(cbind(res1_FC, res2_FC),2)
aicc_FC <- round(cbind(ll_FC,waic_FC),1)
para_list_FC <- do.call(rbind, para_list)

# export
write.csv(percentile_FC, "percentile_FC.csv")
write.csv(MeanSD_FC, "MeanSD_FC.csv")  # for Table S1
write.csv(aicc_FC, "aicc_FC.csv")  # for Table S1
write.csv(para_list_FC, "parameters_FC.csv")

#### 4.  Preparation for visualization ----
### 4-1. UV----
# load data
# load('Stan_UV_n2.RData')
# load('Stan_UV_n4.RData')
# load('Stan_UV_n6.RData')

# make lists to store the results
pt95 <- list()
pt99 <- list()
pt999 <- list()
pt95_long <- list()
pt99_long <- list()
pt999_long <- list()
gl_pred <- list()
Gam_plot <- list()
Wei_plot <- list()
ln_plot <- list()
ecdf_plot <- list()

# guideline value
gl_UV 
#46.5

for (i in c(2,4,6)){
  pos_name <- paste0("pos_UV_",i)
  pos <- get(pos_name)

  # for violin plot
  pt95[[i]] <- data.frame(
  Gamma = qgamma(p = 0.95, shape = pos$gamma[,1], rate = pos$gamma[,2]),
  Weibull = qweibull(p = 0.95, shape = pos$weibull[,1], scale = pos$weibull[,2]),
  Lognormal = qlnorm(p = 0.95, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]))%>%
    mutate(nlog = rep(i, n()))
  
  pt99[[i]] <- data.frame(
    Gamma = qgamma(p = 0.99, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.99, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.99, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]))%>%
    mutate(nlog = rep(i, n()))
  
  pt999[[i]] <- data.frame(
    Gamma = qgamma(p = 0.999, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.999, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.999, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2])) %>%
    mutate(nlog = rep(i, n()))

  gl_pred[[i]] <- data.frame(
    Gamma = pgamma(q = 4* as.numeric(gl_UV) , shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = pweibull(q = 4* as.numeric(gl_UV) , shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = plnorm(q = 4* as.numeric(gl_UV) , meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2])) %>%
    mutate(nlog = rep(i, n()))
  
  # for ssd plot
  logx <- seq(-4,5, by=0.1)
  x_plot <- 10^(logx)
  
  Gam_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                               pred= sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.5))),
                               low = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025))),
                               upp = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  Wei_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                 pred= sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.5))),
                                 low = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025))),
                                 upp = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  ln_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                pred= sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.5))),
                                low = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025))),
                                upp = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  
  ecdf_plot[[i]] <- data.frame(
    #Take mean-k values to draw empirical CDF
    Mean = i*log(10)/(raw_data_uv$mean_k),
    lower = i*log(10)/(raw_data_uv$mean_k + 2*raw_data_uv$sd_k),
    upper = ifelse(raw_data_uv$mean_k < 2*raw_data_uv$sd_k, 
                   NA, i*log(10)/(raw_data_uv$mean_k - 2*raw_data_uv$sd_k)))%>%
    mutate(nlog = rep(i, n()))
  
  l <- ecdf(ecdf_plot[[i]]$lower)
  u <- ecdf(ecdf_plot[[i]]$upper)
  m <- ecdf(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$lower1 <- l(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$upper1 <- u(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$mean1 <- m(ecdf_plot[[i]]$Mean)
}

# prepare for fig 2b
gl_pred_UV <- do.call(rbind, gl_pred)
gl_pred_UV_long <- pivot_longer(gl_pred_UV, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )

# prepare for fig3
pt95_UV <- do.call(rbind, pt95)
pt99_UV <- do.call(rbind, pt99)
pt999_UV <- do.call(rbind, pt999)

pt95_UV_long <- pivot_longer(pt95_UV, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" ) %>%
  mutate(pt = rep(95, n()))
pt99_UV_long <- pivot_longer(pt99_UV, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )%>%
  mutate(pt = rep(99, n()))
pt999_UV_long <- pivot_longer(pt999_UV, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )%>%
  mutate(pt = rep(99.9, n()))

violin_UV <- bind_rows(pt95_UV_long,pt99_UV_long,pt999_UV_long)

# prepare for fig 1 and 3
Gam_plot_UV <- do.call(rbind, Gam_plot)
Wei_plot_UV <- do.call(rbind, Wei_plot)
ln_plot_UV <- do.call(rbind, ln_plot)
ecdf_plot_UV <- do.call(rbind, ecdf_plot)
ecdf_plot_UV_2 <- ecdf_plot_UV %>%
  mutate(Species = rep(raw_data_uv$Species,3),
         EV = rep(raw_data_uv$EntericVirus,3),
         upper_adj = ifelse(is.na(upper), Mean, upper),
         Group = ifelse(EV == TRUE, "Waterborne viruses", "Other viruses"))
ecdf_plot_UV_2$Group <- factor(ecdf_plot_UV_2$Group, levels = c("Waterborne viruses","Other viruses"))

# save
write.csv(gl_pred_UV_long, "gl_pred_UV_long.csv",row.names = FALSE)
write.csv(violin_UV, "violin_UV.csv",row.names = FALSE)
write.csv(Gam_plot_UV, "Gam_plot_UV.csv",row.names = FALSE)
write.csv(Wei_plot_UV, "Wei_plot_UV.csv",row.names = FALSE)
write.csv(ln_plot_UV, "ln_plot_UV.csv",row.names = FALSE)
write.csv(ecdf_plot_UV_2, "ecdf_plot_UV_2.csv",row.names = FALSE)

### 4-2. Ozone----
# load data
# load('Stan_OZ_n2.RData')
# load('Stan_OZ_n4.RData')
# load('Stan_OZ_n6.RData')

# make lists to store the results
pt95 <- list()
pt99 <- list()
pt999 <- list()
pt95_long <- list()
pt99_long <- list()
pt999_long <- list()
gl_pred <- list()
Gam_plot <- list()
Wei_plot <- list()
ln_plot <- list()
ecdf_plot <- list()

# guideline value
gl_OZ 
#0.113

for (i in c(2,4,6)){
  pos_name <- paste0("pos_OZ_",i)
  pos <- get(pos_name)
  
  # for violin plot
  pt95[[i]] <- data.frame(
    Gamma = qgamma(p = 0.95, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.95, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.95, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]))%>%
    mutate(nlog = rep(i, n()))
  
  pt99[[i]] <- data.frame(
    Gamma = qgamma(p = 0.99, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.99, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.99, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]))%>%
    mutate(nlog = rep(i, n()))
  
  pt999[[i]] <- data.frame(
    Gamma = qgamma(p = 0.999, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.999, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.999, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2])) %>%
    mutate(nlog = rep(i, n()))
  
  gl_pred[[i]] <- data.frame(
    Gamma = pgamma(q = 4* as.numeric(gl_OZ) , shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = pweibull(q = 4* as.numeric(gl_OZ) , shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = plnorm(q = 4* as.numeric(gl_OZ) , meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2])) %>%
    mutate(nlog = rep(i, n()))
  
  # for ssd plot
  logx <- seq(-12,4, by=0.1)
  x_plot <- 10^(logx)
  
  Gam_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                      pred= sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.5))),
                                      low = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025))),
                                      upp = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  Wei_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                      pred= sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.5))),
                                      low = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025))),
                                      upp = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  ln_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                     pred= sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.5))),
                                     low = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025))),
                                     upp = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  
  ecdf_plot[[i]] <- data.frame(
    #Take mean-k values to draw empirical CDF
    Mean = i*log(10)/(raw_data_ozone$mean_k),
    lower = i*log(10)/(raw_data_ozone$mean_k + 2*raw_data_ozone$sd_k),
    upper = ifelse(raw_data_ozone$mean_k < 2*raw_data_ozone$sd_k, 
                   NA, i*log(10)/(raw_data_ozone$mean_k - 2*raw_data_ozone$sd_k)))%>%
    mutate(nlog = rep(i, n()))
  
  l <- ecdf(ecdf_plot[[i]]$lower)
  u <- ecdf(ecdf_plot[[i]]$upper)
  m <- ecdf(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$lower1 <- l(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$upper1 <- u(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$mean1 <- m(ecdf_plot[[i]]$Mean)
}

# prepare for fig 2b
gl_pred_OZ <- do.call(rbind, gl_pred)
gl_pred_OZ_long <- pivot_longer(gl_pred_OZ, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )

# prepare for fig3
pt95_OZ <- do.call(rbind, pt95)
pt99_OZ <- do.call(rbind, pt99)
pt999_OZ <- do.call(rbind, pt999)

pt95_OZ_long <- pivot_longer(pt95_OZ, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" ) %>%
  mutate(pt = rep(95, n()))
pt99_OZ_long <- pivot_longer(pt99_OZ, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )%>%
  mutate(pt = rep(99, n()))
pt999_OZ_long <- pivot_longer(pt999_OZ, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )%>%
  mutate(pt = rep(99.9, n()))

violin_OZ <- bind_rows(pt95_OZ_long,pt99_OZ_long,pt999_OZ_long)

# prepare for fig 1 and 3
Gam_plot_OZ <- do.call(rbind, Gam_plot)
Wei_plot_OZ <- do.call(rbind, Wei_plot)
ln_plot_OZ <- do.call(rbind, ln_plot)
ecdf_plot_OZ <- do.call(rbind, ecdf_plot)

ecdf_plot_OZ_2 <- ecdf_plot_OZ %>%
  mutate(Species = rep(raw_data_ozone$Species,3),
         EV = rep(raw_data_ozone$EntericVirus,3),
         upper_adj = ifelse(is.na(upper), Mean, upper),
         Group = ifelse(EV == TRUE, "Waterborne viruses", "Other viruses"))
ecdf_plot_OZ_2$Group <- factor(ecdf_plot_OZ_2$Group, levels = c("Waterborne viruses","Other viruses"))

# save
write.csv(gl_pred_OZ_long, "gl_pred_OZ_long.csv",row.names = FALSE)
write.csv(violin_OZ, "violin_OZ.csv",row.names = FALSE)
write.csv(Gam_plot_OZ, "Gam_plot_OZ.csv",row.names = FALSE)
write.csv(Wei_plot_OZ, "Wei_plot_OZ.csv",row.names = FALSE)
write.csv(ln_plot_OZ, "ln_plot_OZ.csv",row.names = FALSE)
write.csv(ecdf_plot_OZ_2, "ecdf_plot_OZ_2.csv",row.names = FALSE)

### 4-3. Free chlorine ----
# load data
# load('Stan_FC_n2.RData')
# load('Stan_FC_n4.RData')
# load('Stan_FC_n6.RData')

# make lists to store the results
pt95 <- list()
pt99 <- list()
pt999 <- list()
pt95_long <- list()
pt99_long <- list()
pt999_long <- list()
gl_pred <- list()
Gam_plot <- list()
Wei_plot <- list()
ln_plot <- list()
ecdf_plot <- list()

# guideline value
gl_FC 
#0.75

for (i in c(2,4,6)){
  pos_name <- paste0("pos_FC_",i)
  pos <- get(pos_name)
  
  # for violin plot
  pt95[[i]] <- data.frame(
    Gamma = qgamma(p = 0.95, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.95, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.95, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]))%>%
    mutate(nlog = rep(i, n()))
  
  pt99[[i]] <- data.frame(
    Gamma = qgamma(p = 0.99, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.99, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.99, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]))%>%
    mutate(nlog = rep(i, n()))
  
  pt999[[i]] <- data.frame(
    Gamma = qgamma(p = 0.999, shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = qweibull(p = 0.999, shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = qlnorm(p = 0.999, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2])) %>%
    mutate(nlog = rep(i, n()))
  
  gl_pred[[i]] <- data.frame(
    Gamma = pgamma(q = 4* as.numeric(gl_FC) , shape = pos$gamma[,1], rate = pos$gamma[,2]),
    Weibull = pweibull(q = 4* as.numeric(gl_FC) , shape = pos$weibull[,1], scale = pos$weibull[,2]),
    Lognormal = plnorm(q = 4* as.numeric(gl_FC) , meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2])) %>%
    mutate(nlog = rep(i, n()))
  
  # for ssd plot
  logx <- seq(-6,4, by=0.1)
  x_plot <- 10^(logx)
  
  Gam_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                      pred= sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.5))),
                                      low = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.025))),
                                      upp = sapply(x_plot, function(q) quantile(pgamma(q = q, shape = pos$gamma[,1], rate = pos$gamma[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  Wei_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                      pred= sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.5))),
                                      low = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.025))),
                                      upp = sapply(x_plot, function(q) quantile(pweibull(q = q, shape = pos$weibull[,1], scale = pos$weibull[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  ln_plot[[i]] <- as.data.frame(list(dose= x_plot, 
                                     pred= sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.5))),
                                     low = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.025))),
                                     upp = sapply(x_plot, function(q) quantile(plnorm(q = q, meanlog = pos$lognormal[,1], sdlog= pos$lognormal[,2]), probs = c(0.975)))))%>%
    mutate(nlog = rep(i, n()))
  
  ecdf_plot[[i]] <- data.frame(
    #Take mean-k values to draw empirical CDF
    Mean = i*log(10)/(raw_data_fc$mean_k),
    lower = i*log(10)/(raw_data_fc$mean_k + 2*raw_data_fc$sd_k),
    upper = ifelse(raw_data_fc$mean_k < 2*raw_data_fc$sd_k, 
                   NA, i*log(10)/(raw_data_fc$mean_k - 2*raw_data_fc$sd_k)))%>%
    mutate(nlog = rep(i, n()))
  
  l <- ecdf(ecdf_plot[[i]]$lower)
  u <- ecdf(ecdf_plot[[i]]$upper)
  m <- ecdf(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$lower1 <- l(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$upper1 <- u(ecdf_plot[[i]]$Mean)
  ecdf_plot[[i]]$mean1 <- m(ecdf_plot[[i]]$Mean)
}

# prepare for fig 2b
gl_pred_FC <- do.call(rbind, gl_pred)
gl_pred_FC_long <- pivot_longer(gl_pred_FC, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )

# prepare for fig3
pt95_FC <- do.call(rbind, pt95)
pt99_FC <- do.call(rbind, pt99)
pt999_FC <- do.call(rbind, pt999)

pt95_FC_long <- pivot_longer(pt95_FC, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" ) %>%
  mutate(pt = rep(95, n()))
pt99_FC_long <- pivot_longer(pt99_FC, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )%>%
  mutate(pt = rep(99, n()))
pt999_FC_long <- pivot_longer(pt999_FC, cols = c("Gamma", "Weibull", "Lognormal"), names_to = "Distribution", values_to = "value" )%>%
  mutate(pt = rep(99.9, n()))

violin_FC <- bind_rows(pt95_FC_long,pt99_FC_long,pt999_FC_long)

# prepare for fig 1 and 3
Gam_plot_FC <- do.call(rbind, Gam_plot)
Wei_plot_FC <- do.call(rbind, Wei_plot)
ln_plot_FC <- do.call(rbind, ln_plot)
ecdf_plot_FC <- do.call(rbind, ecdf_plot)

ecdf_plot_FC_2 <- ecdf_plot_FC %>%
  mutate(Species = rep(raw_data_fc$Species,3),
         EV = rep(raw_data_fc$EntericVirus,3),
         upper_adj = ifelse(is.na(upper), Mean, upper),
         Group = ifelse(EV == TRUE, "Waterborne viruses", "Other viruses"))

ecdf_plot_FC_2$Group <- factor(ecdf_plot_FC_2$Group, levels = c("Waterborne viruses","Other viruses"))

# save
write.csv(gl_pred_FC_long, "gl_pred_FC_long.csv",row.names = FALSE)
write.csv(violin_FC, "violin_FC.csv",row.names = FALSE)
write.csv(Gam_plot_FC, "Gam_plot_FC.csv",row.names = FALSE)
write.csv(Wei_plot_FC, "Wei_plot_FC.csv",row.names = FALSE)
write.csv(ln_plot_FC, "ln_plot_FC.csv",row.names = FALSE)
write.csv(ecdf_plot_FC_2, "ecdf_plot_FC2.csv",row.names = FALSE)

#### 5. Visualization  ----
### Figure S1 ----

hist_UV <- raw_data_uv %>%
  mutate(inv_k = 4*log(10)/mean_k) %>%
  ggplot(aes(x = inv_k)) +
  geom_histogram(binwidth = 0.1, fill = "lightpink3", color = "black", alpha=0.2)+
  theme_bw(base_size = 15)+
  scale_x_log10()+
  ylim(0,6) +
  labs(title = "UV", x = expression("Dose [mJ" ~ cm^-2 * "]"), y = "Frequency")

hist_OZ <- raw_data_ozone %>%
  mutate(inv_k = 4*log(10)/mean_k) %>%
  ggplot(aes(x = inv_k)) +
  geom_histogram(binwidth = 0.1, fill = "blue3", color = "black", alpha=0.2)+
  theme_bw(base_size = 15)+
  scale_x_log10()+
  ylim(0,6) +
  labs(title = "Ozone", x = expression("Dose [mg min" ~ L^-1 * "]"), y = "Frequency")

hist_FC <- raw_data_fc %>%
  mutate(inv_k = 4*log(10)/mean_k) %>%
  ggplot(aes(x = inv_k)) +
  geom_histogram(binwidth = 0.1, fill = "olivedrab3", color = "black", alpha=0.2)+
  theme_bw(base_size = 15)+
  scale_x_log10()+
  ylim(0,6) +
  labs(title = "Free chlorine", x = expression("Dose [mg min" ~ L^-1 * "]"), y = "Frequency")

figS1 <- cowplot::plot_grid(hist_UV, hist_OZ, hist_FC, nrow=1, labels="AUTO")

ggsave("FigS1_v1.png", plot = figS1,  width=9, height=3, dpi = 1200)

### Figure 1 (only for lognormal)----
# nlog = 4

fig1_ln_UV_n4 <-  ggplot() +
  geom_point(data =ecdf_plot_UV_2 %>% filter(nlog == 4),
             aes(x = Mean, y = mean1, color = Group, shape = Group), size = 4) +
  scale_color_manual(values = c("Waterborne viruses" = "springgreen2", "Other viruses" = "yellow4")) + 
  scale_shape_manual(values = c("Waterborne viruses" = 17, "Other viruses" = 16)) + 
  geom_text(data =ecdf_plot_UV_2 %>% filter(nlog == 4),
            aes(x = Mean, y = mean1, label = Species,fontface = "italic", vjust = -0.2, hjust = 1.2), size = 3)+
  geom_errorbarh(data =ecdf_plot_UV_2 %>% filter(nlog == 4) %>% filter(Mean != lower & !is.na(upper)),
                  aes(xmin = lower, xmax = upper, y =mean1))+
  geom_line(data =ln_plot_UV %>% filter(nlog == 4),aes(x=dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], linewidth=1) +
  geom_ribbon(data =ln_plot_UV %>% filter(nlog == 4),aes(x=dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  geom_vline(xintercept = gl_UV*4, linetype = "dashed", color = "grey", size=1) +
  theme_bw(base_size = 15)+
  theme(legend.position = 'bottom')+
  scale_x_continuous(trans = 'log10',
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x=expression("Dose [mJ" ~ cm^-2 * "]"), y = "Potentially inactivated fraction")+
  ggtitle("UV (4-log reduction)")

fig1_ln_OZ_n4 <-  ggplot() +
  geom_point(data =ecdf_plot_OZ_2 %>% filter(nlog == 4), 
             aes(x = Mean, y = mean1, color = Group, shape = Group), size = 4) +
  scale_color_manual(values = c("Waterborne viruses" = "springgreen2", "Other viruses" = "yellow4")) + 
  scale_shape_manual(values = c("Waterborne viruses" = 17, "Other viruses" = 16)) + 
  geom_text(data =ecdf_plot_OZ_2 %>% filter(nlog == 4),
            aes(x = Mean, y = mean1, label = Species,fontface = "italic", vjust = -0.2, hjust = 1.2), size = 3)+
  geom_errorbarh(data =ecdf_plot_OZ_2 %>% filter(nlog == 4) %>% filter(Mean != lower & !is.na(upper)),
                 aes(xmin = lower, xmax = upper, y =mean1))+
  geom_line(data =ln_plot_OZ %>% filter(nlog == 4),aes(x=dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], linewidth=1) +
  geom_ribbon(data =ln_plot_OZ %>% filter(nlog == 4),aes(x=dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  geom_vline(xintercept = gl_OZ*4, linetype = "dashed", color = "grey", size=1) +
  theme_bw(base_size = 15)+
  scale_x_continuous(trans = 'log10',
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction")+
  ggtitle("Ozone (4-log reduction)")

fig1_ln_FC_n4 <-  ggplot() +
  geom_point(data =ecdf_plot_FC_2 %>% filter(nlog == 4), 
             aes(x = Mean, y = mean1, color = Group, shape = Group), size = 4) +
  scale_color_manual(values = c("Waterborne viruses" = "springgreen2", "Other viruses" = "yellow4")) + 
  scale_shape_manual(values = c("Waterborne viruses" = 17, "Other viruses" = 16)) + 
  geom_text(data =ecdf_plot_FC_2 %>% filter(nlog == 4),
            aes(x = Mean, y = mean1, label = Species,fontface = "italic", vjust = -0.2, hjust = 1.2), size = 3)+
  geom_errorbarh(data =ecdf_plot_FC_2 %>% filter(nlog == 4) %>% filter(Mean != lower & !is.na(upper)),
                 aes(xmin = lower, xmax = upper, y =mean1))+
  geom_line(data =ln_plot_FC %>% filter(nlog == 4),aes(x=dose, y=pred), color=RColorBrewer::brewer.pal(11, "RdBu")[11], linewidth=1) +
  geom_ribbon(data =ln_plot_FC %>% filter(nlog == 4),aes(x=dose,ymin=low,ymax=upp), fill = RColorBrewer::brewer.pal(11, "RdBu")[11], alpha=0.1) +
  geom_vline(xintercept = gl_FC*4, linetype = "dashed", color = "grey", size=1) +
  theme_bw(base_size = 15)+
  scale_x_continuous(trans = 'log10',
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction")+
  ggtitle("Free chlorine (4-log reduction)")

legend_fig1 <- cowplot::get_plot_component(fig1_ln_UV_n4, 'guide-box-bottom', return_all = TRUE)
fig1a <- fig1_ln_UV_n4 + theme(legend.position = "none")
fig1b <- fig1_ln_OZ_n4 + theme(legend.position = "none")
fig1c <- fig1_ln_FC_n4 + theme(legend.position = "none")
fig1abc <- cowplot::plot_grid(fig1a, fig1b, fig1c, nrow=1, labels="AUTO")

# Figure 1
Figure1 <- cowplot::plot_grid(fig1abc,legend_fig1,nrow = 2,  rel_heights = c(5, 0.4))
ggsave("Fig1_v1.png", plot = Figure1,  width=16, height=6, dpi = 1200)

### Figure S2 (UV)  ----

fig_ln_UV_n246 <- ln_plot_UV %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("snow3", "lightpink3", "orchid4")) + 
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  theme_bw(base_size = 15)+
  theme(legend.position = 'bottom')+
  scale_x_log10()+
  labs(x=expression("Dose [mJ" ~ cm^-2 * "]"), y = "Potentially inactivated fraction",
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  #guides(linetype = FALSE)+
  coord_cartesian(xlim = c(0.1, 10000)) +
  ggtitle("Lognormal")

fig_wei_UV_n246 <- Wei_plot_UV %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("snow3", "lightpink3", "orchid4")) + 
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  theme_bw(base_size = 15)+
  scale_x_log10()+
  labs(x=expression("Dose [mJ" ~ cm^-2 * "]"), y = "Potentially inactivated fraction",
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  #guides(linetype = FALSE)+
  coord_cartesian(xlim = c(0.1, 10000)) +
  ggtitle("Weibull")

fig_gam_UV_n246 <- Gam_plot_UV %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("snow3", "lightpink3", "orchid4")) + 
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  theme_bw(base_size = 15)+
  scale_x_log10()+
  labs(x=expression("Dose [mJ" ~ cm^-2 * "]"), y = "Potentially inactivated fraction",
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  #guides(linetype = FALSE)+
  coord_cartesian(xlim = c(0.1, 10000)) +
  ggtitle("Gamma")

legend_figS2 <- cowplot::get_plot_component(fig_ln_UV_n246, 'guide-box-bottom', return_all = TRUE)

figS2a <- fig_ln_UV_n246 + theme(legend.position = "none")
figS2b <- fig_wei_UV_n246 + theme(legend.position = "none")
figS2c <- fig_gam_UV_n246 + theme(legend.position = "none")
figS2abc <- cowplot::plot_grid(figS2a, figS2b, figS2c, nrow=1, labels="AUTO")

# Figure S2
FigureS2 <- cowplot::plot_grid(figS2abc,legend_figS2, nrow=2,  rel_heights = c(5, 0.4))
ggsave("FigS2_v1.png", plot = FigureS2,  width=16, height=6, dpi = 1200)

### Figure S3 (ozone)----
fig_ln_OZ_n246 <- ln_plot_OZ %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("darkturquoise", "blue3", "ivory4")) + 
  scale_fill_manual(values = c("darkturquoise", "blue3", "ivory4")) +
  theme_bw(base_size = 15)+
  theme(legend.position = 'bottom')+
  scale_x_log10()+
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  coord_cartesian(xlim = c(10^(-12), 10000)) +
  ggtitle("Lognormal")

fig_wei_OZ_n246 <- Wei_plot_OZ %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("darkturquoise", "blue3", "ivory4")) + 
  scale_fill_manual(values = c("darkturquoise", "blue3", "ivory4")) +
  theme_bw(base_size = 15)+
  scale_x_log10()+
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  coord_cartesian(xlim = c(10^(-12), 10000)) +
  ggtitle("Weibull")

fig_gam_OZ_n246 <- Gam_plot_OZ %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("darkturquoise", "blue3", "ivory4")) + 
  scale_fill_manual(values = c("darkturquoise", "blue3", "ivory4")) +
  theme_bw(base_size = 15)+
  scale_x_log10()+
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  coord_cartesian(xlim = c(10^(-12), 10000)) +
  ggtitle("Gamma")

legend_figS3 <- cowplot::get_plot_component(fig_ln_OZ_n246, 'guide-box-bottom', return_all = TRUE)

figS3a <- fig_ln_OZ_n246 + theme(legend.position = "none")
figS3b <- fig_wei_OZ_n246 + theme(legend.position = "none")
figS3c <- fig_gam_OZ_n246 + theme(legend.position = "none")
figS3abc <- cowplot::plot_grid(figS3a, figS3b, figS3c, nrow=1, labels="AUTO")

# Figure S3
FigureS3 <- cowplot::plot_grid(figS3abc,legend_figS3,nrow = 2,  rel_heights = c(5, 0.4))
ggsave("FigS3_v1.png", plot = FigureS3,  width=16, height=6, dpi = 1200)

### Figure S4 (FC) ----

fig_ln_FC_n246 <- ln_plot_FC %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("khaki", "olivedrab3", "seagreen")) + 
  scale_fill_manual(values = c("khaki", "olivedrab3", "seagreen")) +
  theme_bw(base_size = 15)+
  theme(legend.position = 'bottom')+
  scale_x_log10()+
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  coord_cartesian(xlim = c(10^(-6), 10000)) +
  ggtitle("Lognormal")

fig_wei_FC_n246 <- Wei_plot_FC %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("khaki", "olivedrab3", "seagreen")) + 
  scale_fill_manual(values = c("khaki", "olivedrab3", "seagreen")) +
  theme_bw(base_size = 15)+
  scale_x_log10()+
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  coord_cartesian(xlim = c(10^(-6), 10000)) +
  ggtitle("Weibull")

fig_gam_FC_n246 <- Gam_plot_FC %>% 
  ggplot() +
  geom_line(aes(x=dose, y=pred, color=as.factor(nlog), linetype = as.factor(nlog)), size=1.5) +
  geom_ribbon(aes(x=dose,ymin=low,ymax=upp, fill=as.factor(nlog)), alpha=0.2) +
  scale_color_manual(values = c("khaki", "olivedrab3", "seagreen")) + 
  scale_fill_manual(values = c("khaki", "olivedrab3", "seagreen")) +
  theme_bw(base_size = 15)+
  scale_x_log10()+
  labs(x=expression("Dose [mg min" ~ L^-1 * "]"), y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction",linetype = "n-log reduction")+
  coord_cartesian(xlim = c(10^(-6), 10000)) +
  ggtitle("Gamma")

legend_figS4 <- cowplot::get_plot_component(fig_ln_FC_n246, 'guide-box-bottom', return_all = TRUE)

figS4a <- fig_ln_FC_n246 + theme(legend.position = "none")
figS4b <- fig_wei_FC_n246 + theme(legend.position = "none")
figS4c <- fig_gam_FC_n246 + theme(legend.position = "none")
figS4abc <- cowplot::plot_grid(figS4a, figS4b, figS4c, nrow=1, labels="AUTO")

# Figure S4
FigureS4 <- cowplot::plot_grid(figS4abc,legend_figS4,nrow = 2,  rel_heights = c(5, 0.4))
ggsave("FigS4_v1.png", plot = FigureS4,  width=16, height=6, dpi = 1200)

### Figure 2 ----

gl_UV_plot_ln <- gl_pred_UV_long %>%
  filter(Distribution == "Lognormal") %>%
  ggplot(aes(x = as.factor(nlog), y = value,fill=as.factor(nlog)))+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  theme_bw(base_size = 15)+
  labs(x="n-log reduction", y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction")+
  coord_cartesian(ylim = c(0.4, 1)) +
  ggtitle("UV (n-log reduction)")

gl_OZ_plot_ln <- gl_pred_OZ_long %>%
  filter(Distribution == "Lognormal") %>%
  ggplot(aes(x = as.factor(nlog), y = value,fill=as.factor(nlog)))+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  theme_bw(base_size = 15)+
  labs(x="n-log reduction", y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction")+
  coord_cartesian(ylim = c(0.4, 1)) +
  ggtitle("Ozone (n-log reduction)")

gl_FC_plot_ln <- gl_pred_FC_long %>%
  filter(Distribution == "Lognormal") %>%
  ggplot(aes(x = as.factor(nlog), y = value,fill=as.factor(nlog)))+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  theme_bw(base_size = 15)+
  labs(x="n-log reduction", y = "Potentially inactivated fraction", 
       color = "n-log reduction", fill = "n-log reduction")+
  coord_cartesian(ylim = c(0.4, 1)) +
  ggtitle("Free chlorine (n-log reduction)")

legend_fig2 <- cowplot::get_plot_component(fig_ln_UV_n246, 'guide-box-bottom', return_all = TRUE)

fig2a <- fig_ln_UV_n246 + theme(legend.position = "none")

fig2b <- gl_UV_plot_ln + theme(legend.position = "none",
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank()) +
  ggtitle("B) UV (n-log reduction)")
fig2c <- gl_OZ_plot_ln + theme(legend.position = "none",
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank()) +
  ggtitle("C) Ozone (n-log reduction)")
fig2d <- gl_FC_plot_ln + theme(legend.position = "none",
                               axis.title.x = element_blank(),
                               axis.title.y = element_blank()) +
  ggtitle("D) Free chlorine (n-log reduction)")

labeled_fig2a <- fig2a + labs(title = "A) UV sensitivity distributions for n-log reduction")
combined_plot <- plot_grid(fig2b, fig2c, fig2d, 
                           ncol = 1, align = "v", rel_heights = c(1, 1, 1))
x_label <- ggdraw() + draw_label("n-log reduction", hjust = 0.5, vjust = 1, size = 14)
y_label <- ggdraw() + draw_label("Potentially inactivated fraction", angle = 90, hjust = 0.5, vjust = 1, size = 14)
fig2cbd <- plot_grid(
  y_label,
  plot_grid(
    combined_plot,
    x_label,
    ncol = 1,
    rel_heights = c(1, 0.06)
  ),
  ncol = 2,
  rel_widths = c(0.05, 1))

fig2abcd_all <- (
  (labeled_fig2a | fig2cbd) /
    legend_fig2)

# Figure 2
Figure2 <- fig2abcd_all + plot_layout(widths = c(1, 1), heights = c(5,0.5)) 
ggsave("Fig2_v1.png", plot = Figure2,  width=16, height=8, dpi = 1200)

### Figure 3 ----
violin_UV_plot_ln <- violin_UV  %>%
  filter(Distribution == "Lognormal") %>%
  mutate(across(nlog, ~factor(., levels = c("2", "4", "6")))) %>%
  mutate(across(pt, ~factor(., levels = c("95", "99", "99.9")))) %>%
  ggplot(aes(x = nlog, y = value, fill=as.factor(nlog)))+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  facet_wrap(~ pt, labeller = labeller(pt = c("95" = "95th", "99" = "99th", "99.9" = "99.9th"))) +
  scale_y_log10()+
  theme_bw(base_size = 18)+
  theme(legend.position = 'bottom')+
  labs(x="n-log reduction", y = expression("Dose [mJ" ~ cm^-2 * "]"), 
       color = "n-log reduction", fill = "n-log reduction")+
  coord_cartesian(ylim = c(0.01, 10^12)) +
  ggtitle("UV (n-log reduction)")+
  theme(strip.background = element_rect(
      color="black", fill="white"))

violin_OZ_plot_ln <- violin_OZ  %>%
  filter(Distribution == "Lognormal") %>%
  mutate(across(nlog, ~factor(., levels = c("2", "4", "6")))) %>%
  mutate(across(pt, ~factor(., levels = c("95", "99", "99.9")))) %>%
  ggplot(aes(x = nlog, y = value, fill=as.factor(nlog)))+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  facet_wrap(~ pt, labeller = labeller(pt = c("95" = "95th", "99" = "99th", "99.9" = "99.9th"))) +
  scale_y_log10()+
  theme_bw(base_size = 18)+
  labs(x="n-log reduction", y = expression("Dose [mg min" ~ L^-1 * "]"), 
       color = "n-log reduction", fill = "n-log reduction")+
  coord_cartesian(ylim = c(0.01, 10^12)) +
  ggtitle("Ozone (n-log reduction)")+
  theme(strip.background = element_rect(
    color="black", fill="white"))

violin_FC_plot_ln <- violin_FC  %>%
  filter(Distribution == "Lognormal") %>%
  mutate(across(nlog, ~factor(., levels = c("2", "4", "6")))) %>%
  mutate(across(pt, ~factor(., levels = c("95", "99", "99.9")))) %>%
  ggplot(aes(x = nlog, y = value, fill=as.factor(nlog)))+
  geom_violin(alpha=0.5)+
  scale_fill_manual(values = c("snow3", "lightpink3", "orchid4")) +
  facet_wrap(~ pt, labeller = labeller(pt = c("95" = "95th", "99" = "99th", "99.9" = "99.9th"))) +
  scale_y_log10()+
  theme_bw(base_size = 18)+
  labs(x="n-log reduction", y = expression("Dose [mg min" ~ L^-1 * "]"), 
       color = "n-log reduction", fill = "n-log reduction")+
  coord_cartesian(ylim = c(0.01, 10^12)) +
  ggtitle("Free chlorine (n-log reduction)")+
  theme(strip.background = element_rect(
    color="black", fill="white"))

legend_fig3 <- cowplot::get_plot_component(violin_UV_plot_ln, 'guide-box-bottom', return_all = TRUE)

fig3a <- violin_UV_plot_ln + theme(legend.position = "none")
fig3b <- violin_OZ_plot_ln + theme(legend.position = "none")
fig3c <- violin_FC_plot_ln + theme(legend.position = "none")
fig3abc <- cowplot::plot_grid(fig3a, fig3b, fig3c, nrow=1, labels="AUTO")

# Figure 3
Figure3 <- cowplot::plot_grid(fig3abc,legend_fig3,nrow = 2,  rel_heights = c(5, 0.4))
ggsave("Fig3_v1.png", plot = Figure3,  width=16, height=6, dpi = 1200)


