#import libraries
library(tidyverse)
library(mlmRev)
library(GGally)
library(rstanarm)
library(bayesplot) 
library(loo)
library(caret)
library(rstan)

#import and prepare data
seed=3
data = Mmmec
data
deaths=data$deaths
groupingL3=data$nation
groupingL2=data$region
groupingL1=data$county
#preparing data transforming factors in integers
data2=data
data2$region_id = as.integer(factor(data$region, levels = unique(data$region)))
data2$nation_id = as.integer(factor(data$nation, levels = unique(data$nation)))
nationLookupVec <- as.integer(unique(data2[c("region","nation")])[,"nation"])

stan_data =
  with(data2,
       list(
         N = length(deaths),
         L = length(unique(groupingL2)),
         J = length(unique(groupingL3)),
         deaths = deaths,
         uvb = uvb,
         population = log(expected),
         region_id = region_id,
         nation_id = nation_id,
         nationWithinRegion = nationLookupVec 
       )
  )


#import  and compile model

model_a_stan = rstan::stan_model("models/modelA.stan")
model_b_stan = rstan::stan_model("models/modelB.stan")
#fit models
fitted_a_stan <- rstan::sampling(model_a_stan, data = stan_data,
                                 chains = 4, cores = 4, iter = 4000, control = list(adapt_delta = 0.80, max_treedepth=10),verbose=TRUE,seed=seed)

fitted_b_stan <- rstan::sampling(model_b_stan, data = stan_data,
                                 chains = 4, cores = 4, iter = 4000, control = list(adapt_delta = 0.80, max_treedepth=10),verbose=TRUE,seed=5)
saveRDS(fitted_a_stan,"models/modelA.stanModel")
saveRDS(fitted_b_stan,"models/modelB.stanModel")

#loo analysis
library(loo)

#function for loo from stan model
custom_loo = function(fitted_stan_model){
  log_lik_slopes = extract_log_lik(fitted_stan_model)
  loo_slopes <- loo(log_lik_slopes)
  return(loo_slopes)
}

loo_A  = custom_loo(fitted_a_stan)
loo_B  = custom_loo(fitted_b_stan)
print(loo_compare(loo_A,loo_B)) #model B wins

#std residuals modelA
fit = fitted_a_stan
y_rep = as.matrix(fit, pars = "y_rep")
y_rep_values = y_rep[1:length(data$deaths),]
mean_inv_phi <- mean(as.matrix(fit, pars = "inv_phi"))
mean_y_rep <- colMeans(y_rep)
std_resid <- (data$deaths - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
plotData=data.frame(std_resid,mean_y_rep)
ggplot(plotData,aes(x=mean_y_rep,y=std_resid))+
  geom_point()+
  geom_hline(yintercept = 2,color="red")+
  geom_hline(yintercept = -2,color = "red")+
  geom_hline(yintercept = 3,color="green")+
  geom_hline(yintercept = -3,color = "green")+
  ggtitle("Std residuals for model A ")

fit = fitted_b_stan
y_rep = as.matrix(fit, pars = "y_rep")
y_rep_values = y_rep[1:length(data$deaths),]
mean_inv_phi <- mean(as.matrix(fit, pars = "inv_phi"))
mean_y_rep <- colMeans(y_rep)
std_resid <- (data$deaths - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
plotData=data.frame(std_resid,mean_y_rep)
ggplot(plotData,aes(x=mean_y_rep,y=std_resid))+
  geom_point()+
  geom_hline(yintercept = 2,color="red")+
  geom_hline(yintercept = -2,color = "red")+
  geom_hline(yintercept = 3,color="green")+
  geom_hline(yintercept = -3,color = "green")+
  ggtitle("Std residuals for model B ")

#posterior checks model A
fit = fitted_a_stan
paper_results=data.frame(var_nat=0.105,var_reg=0.0147)
mcmc_areas(fit,regex_pars = "var",prob = 0.95)+
  geom_vline(data = paper_results, xintercept = paper_results$var_nat ,color="green")+
  geom_vline(data = paper_results, xintercept = paper_results$var_reg,color="red")+
  annotate("text", x = paper_results$var_nat, y = 1.5, label = "var_nat paper")+
  annotate("text", x = paper_results$var_reg, y = 1.8, label = "var_reg paper")

mcmc_intervals(fit,pars = "alpha")
mcmc_intervals(fit,regex_pars = "sigma")
mcmc_intervals(fit,regex_pars = "^alpha_nat")
mcmc_intervals(fit,regex_pars = "^alpha_reg")

list_of_draws=extract(fit)
dev=data.frame(
  alpha=mean(list_of_draws$alpha),
  alpha_sd=sd(list_of_draws$alpha),
  beta=mean(list_of_draws$beta),
  beta_sd=sd(list_of_draws$beta),
  sigma_nat=mean(list_of_draws$sigma_nat),
  sigma_nat_sd=sd(list_of_draws$sigma_nat),
  sigma_reg=mean(list_of_draws$sigma_reg),
  sigma_reg_sd=sd(list_of_draws$sigma_reg)
)

dev 


#posterior checks model B 
fit = fitted_b_stan
paper_results=data.frame(var_nat=0.105,var_reg=0.0147)
mcmc_areas(fit,regex_pars = "var",prob = 0.95)+
  geom_vline(data = paper_results, xintercept = paper_results$var_nat ,color="green")+
  geom_vline(data = paper_results, xintercept = paper_results$var_reg,color="red")+
  annotate("text", x = paper_results$var_nat, y = 1.5, label = "var_nat paper")+
  annotate("text", x = paper_results$var_reg, y = 1.8, label = "var_reg paper")

mcmc_intervals(fit,pars = "alpha")
mcmc_intervals(fit,regex_pars = "sigma")
mcmc_intervals(fit,regex_pars = "^beta_nat")
mcmc_intervals(fit,regex_pars = "^alpha_nat")
mcmc_intervals(fit,regex_pars = "^alpha_reg")
mcmc_intervals(fit,regex_pars = "^beta_reg")

list_of_draws=extract(fit)
dev=data.frame(
  alpha=mean(list_of_draws$alpha),
  alpha_sd=sd(list_of_draws$alpha),
  beta=mean(list_of_draws$beta),
  beta_sd=sd(list_of_draws$beta),
  sigma_nat=mean(list_of_draws$sigma_nat),
  sigma_nat_sd=sd(list_of_draws$sigma_nat),
  sigma_reg=mean(list_of_draws$sigma_reg),
  sigma_reg_sd=sd(list_of_draws$sigma_reg),
  sigma_beta_reg=mean(list_of_draws$sigma_beta_reg),
  sigma_beta_reg_sd=sd(list_of_draws$sigma_beta_reg),
  sigma_beta_nat=mean(list_of_draws$sigma_beta_nat),
  sigma_beta_nat_sd=sd(list_of_draws$sigma_beta_nat),
  cor_sigma_nat_beta = cor(list_of_draws$sigma_nat,list_of_draws$sigma_beta_nat),
  cor_sigma_reg_beta = cor(list_of_draws$sigma_reg,list_of_draws$sigma_beta_reg)
)

dev 
#correlations
cor(list_of_draws$sigma_nat,list_of_draws$sigma_beta_nat)
cor(list_of_draws$sigma_reg,list_of_draws$sigma_beta_reg)


ppc_dens_overlay(data$deaths,y_rep_values)+
  ggtitle("Replicated vs real")
y=data$deaths
ppc_stat(y,y_rep_values,stat="mean")
ppc_stat(y,y_rep_values,stat="sd")
ppc_stat_grouped(y,y_rep_values,group=data$nation,stat="mean")+
  ggtitle("Using exposures means comparison")
ppc_stat_grouped(y,y_rep_values,group=data$nation,stat="sd")+
  ggtitle("Using exposures sd comparison")
ppc_stat_grouped(y,y_rep_values,group=data$nation,stat="min")+
  ggtitle("Using exposures min comparison")
ppc_stat_grouped(y,y_rep_values,group=data$nation,stat="max")+
  ggtitle("Using exposures max comparison")
ppc_ecdf_overlay(y,y_rep)
ppc_intervals(y,y_rep)
