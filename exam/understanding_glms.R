#UNDERSTANDING GLM'S
library(tidyverse)
library(mlmRev)
library(GGally)
library(rstanarm)
library(bayesplot) 
library(loo)
library(caret)
library(rstan)
#import and arrange data
#prepare data
seed=pi
data = Mmmec
data
deaths=data$deaths
groupingL3=data$nation
groupingL2=data$region
groupingL1=data$county
#super simple glm

a=glm(data=data, formula = deaths ~ 1 + uvb + offset(log(expected)),family=poisson(link="log"))



a=glm(data=data, formula = deaths ~ 1 + uvb + offset(log(expected)),family=MASS::negative.binomial(theta=10,link = "log"))
a=glmer(data=data, formula = deaths ~ (1|groupingL3)+(0 + uvb|groupingL3) + offset(log(expected)),family=MASS::negative.binomial(theta=10,link = "log"))

#model A of the paper
#https://rpsychologist.com/r-guide-longitudinal-lme-lmer#three-level-models three level nesting
#results are compatible with paper, but we have one more dof.
#model A
mod_a=glmer(data=data, formula = deaths ~ 1+(1|groupingL3:groupingL2) + (1|groupingL3)+uvb + offset(log(expected)),family=poisson(link = "log"))
#model B
mod_b=glmer(data=data, formula = deaths ~ 1+(1 + uvb|groupingL3:groupingL2) + (1 + uvb|groupingL3)+uvb + offset(log(expected)),family=poisson(link = "log"))
#model A using negative Binomial
mod_a_neg_bin=glmer(data=data, formula = deaths ~ 1+(1|groupingL3:groupingL2) + (1|groupingL3)+uvb + offset(log(expected)),family=MASS::negative.binomial(theta=10,link = "log"))
#model B using negative Binomial
mod_b_neg_bin=glmer(data=data, formula = deaths ~ 1+(1 + uvb|groupingL3:groupingL2) + (1 + uvb|groupingL3)+uvb + offset(log(expected)),family=MASS::negative.binomial(theta=10,link = "log"))
#same models using bayesian style
#model A using negative Binomial
mod_a_neg_bin_bayes=stan_glmer(data=data, formula = deaths ~ 1+(1|groupingL3:groupingL2) + (1|groupingL3)+uvb + offset(log(expected)),family=neg_binomial_2(link = "log"),cores=4)
#model B using negative Binomial
mod_b_neg_bin_bayes=stan_glmer(data=data, formula = deaths ~ 1+(1 + uvb|groupingL3:groupingL2) + (1 + uvb|groupingL3)+uvb + offset(log(expected)),family=neg_binomial_2(link = "log"),cores=4)
#model A using Stan indpired from https://biologyforfun.wordpress.com/2016/12/08/crossed-and-nested-hierarchical-models-with-stan-and-r/

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
readRDS("models/modelB.stanModel")
model_b_stan = rstan::stan_model("models/modelB.stan")

model_a_poisson_stan = rstan::stan_model("models/modelA_poisson.stan")
#run model
#fitted_a_stan <- rstan::sampling(model_a_stan, data = stan_data,
#                          chains = 4, cores = 4, iter = 4000, control = list(adapt_delta = 0.99, max_treedepth=50),verbose=TRUE,seed=seed)
fitted_a_stan = readRDS("models/modelA.stanModel")

fitted_b_stan <- rstan::sampling(model_b_stan, data = stan_data,
                          chains = 4, cores = 4, iter = 4000, control = list(adapt_delta = 0.90, max_treedepth=40),verbose=TRUE,seed=4)
#fitted_b_stan = readRDS("models/modelB.stanModel")

#fitted_a_poisson_stan <- rstan::sampling(model_a_poisson_stan, data = stan_data,
#                          chains = 4, cores = 4, iter = 1000, control = list(adapt_delta = 0.99, max_treedepth=10),verbose=TRUE,seed=seed)


launch_shinystan(fitted_b_stan)
#saveRDS(fitted_b_stan,"models/modelB.stanModel")
a=as.data.frame(fit,pars=c())


fit=fitted_b_stan
#retrive parameters of interest (sd)
list_of_draws=extract(fit)
sigma_nat=mean(list_of_draws$sigma_nat)
sigma_nat_sd=sd(list_of_draws$sigma_nat)
sigma_reg=mean(list_of_draws$sigma_reg)
sigma_reg_sd=sd(list_of_draws$sigma_reg)
sigma_beta_reg=mean(list_of_draws$sigma_beta_reg)
sigma_beta_reg_sd=sd(list_of_draws$sigma_beta_reg)
sigma_beta_nat=mean(list_of_draws$sigma_beta_nat)
sigma_beta_nat_sd=sd(list_of_draws$sigma_beta_nat)
#fitted parameters
alpha=mean(list_of_draws$alpha)
beta=mean(list_of_draws$beta)
#variations on nation level
beta_nat=colMeans(list_of_draws$beta_nat)
beta_nat_sd = apply(list_of_draws$beta_nat, 2, sd)
alpha_nat=colMeans(list_of_draws$alpha_nat)
alpha_nat_sd = apply(list_of_draws$alpha_nat, 2, sd)
#rate
eta = colMeans(list_of_draws$eta_rep)
#putting all together
coeff = data.frame(nation=unique(data$nation),beta=beta_nat,beta_nat_sd,alpha = alpha_nat,alpha_nat_sd)
#plot figure 2
#plot slopes
ggplot(data=data,aes(x=uvb,y=log(eta),color=nation))+
  geom_point()+
  geom_abline(data=coeff,aes(slope= beta,intercept=alpha,colour=nation))+
  ylim(-1,1)+
  xlim(-9,15)+
  ggtitle("Random slopes and intercepts like in figure 1")

#plots
mcmc_intervals(fit,regex_pars = "sigma")
mcmc_intervals(fit,regex_pars = "^beta_nat")
mcmc_intervals(fit,regex_pars = "^alpha_nat")

cor(as.matrix(list_of_draws$beta_nat),as.matrix(list_of_draws$alpha_nat))

aa=as.matrix(fit,pars="sigma")
summary(fit)
posterior_dev=rstan::get_posterior_mean(fit)
#basic plots
y_rep = as.matrix(fit, pars = "y_rep")
rate = as.matrix(fit, pars = "eta_n")
y_rep_values = y_rep[1:length(data$deaths),]
ppc_dens_overlay(data$deaths,y_rep_values)+
  ggtitle("Rock and Roll")
y=data$deaths
ppc_stat(y,y_rep_values,stat="mean")
ppc_stat(y,y_rep_values,stat="sd")
ppc_stat_grouped(y,y_rep_values,group=data$nation,stat="mean")+
  ggtitle("Using exposures means comparison")
ppc_stat_grouped(y,y_rep_values,group=data$nation,stat="sd")+
  ggtitle("Using exposures sd comparison")
#variance:
std_final=mean(sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi))

#std residuals
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
  ggtitle("Std residuals for fit ")





summary(a)
resid(fit,type="deviance")
ranef(a)
fixef(a)
se(a)
coef(a)
launch_shinystan(a)
AIC(logLik(a))
#getCoeff
coef=coefficients(a)
coef=rownames_to_column(coef$groupingL3,var="nation")
colnames(coef)[colnames(coef)=="uvb"] <- "beta"
colnames(coef)[3] <- "alpha"

data=inner_join(coef,data,by="nation")


ggplot(data=data,aes(x=nation ,y=beta))+
  geom_point()



data$predicted=predict(a,type="response")
data$predictedLambda=predict(a)
data$residuals = data$predicted-data$deaths
data$stdrResiduals = (data$predicted-data$deaths)/data$deaths

#plot slopes
ggplot(data=data,aes(x=uvb,y=log(predictedLambda),color=nation))+
  geom_point()+
  geom_abline(data=coef,aes(slope= beta,intercept=alpha,colour=nation))+
  ylim(-1,1)+
  xlim(-9,15)+
  ggtitle("Random slopes and intercepts like in figure 1")

ggplot(data,aes(x=deaths,y=predicted))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)

ggplot(data,aes(x=deaths,y=residuals))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)

ggplot(data,aes(x=deaths,y=stdrResiduals))+
  geom_point()+
  geom_abline(slope=0,intercept = 0)


#same procedure with best model
a=glmer(data=data, formula = deaths ~ (1|groupingL3)+(0 + uvb|groupingL3) + offset(log(expected)),family=MASS::negative.binomial(theta=10,link = "log"))

a=glmer(data=data, formula = deaths ~ (1|groupingL3:groupingL2) +(0+uvb|groupingL3:groupingL2) + offset(log(expected)),family=MASS::negative.binomial(theta=3*18))

summary(a)
resid(a,type="deviance")
AIC(logLik(a))
#getCoeff
coef=coefficients(a)
coef=rownames_to_column(coef$groupingL3,var="nation")
colnames(coef)[colnames(coef)=="uvb"] <- "beta"
colnames(coef)[3] <- "alpha"

data=inner_join(coef,data,by="nation")


ggplot(data=data,aes(x=nation ,y=beta))+
  geom_point()



data$predicted=predict(a,type="response")
data$predictedLambda=predict(a)
data$residuals = data$predicted-data$deaths
data$stdrResiduals = (data$predicted-data$deaths)/data$deaths

#plot slopes
ggplot(data=data,aes(x=uvb,y=log(predictedLambda),color=nation))+
  geom_point()+
  geom_abline(data=coef,aes(slope= beta,intercept=alpha,colour=nation))+
  ylim(-1,1)+
  xlim(-9,15)+
  ggtitle("Random slopes and intercepts like in figure 2")

ggplot(data,aes(x=deaths,y=predicted))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)

ggplot(data,aes(x=deaths,y=residuals))+
  geom_point()+
  geom_abline(slope=1,intercept = 0)

ggplot(data,aes(x=deaths,y=stdrResiduals))+
  geom_point()+
  geom_abline(slope=0,intercept = 0)


#bayesian variance component analysis
fit6 <- stan_glmer(formula = deaths ~ (1|groupingL3:groupingL2) +uvb + offset(log(expected)), data=data, family=neg_binomial_2(),adapt_delta=0.80,cores=4)
