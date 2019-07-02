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
#import model
model_a_stan = rstan::stan_model("models/modelA.stan")
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

#run model
fitted_a_stan <- rstan::sampling(model_a_stan, data = stan_data,
                          chains = 4, cores = 4, iter = 4000, verbose=TRUE,seed=seed)
launch_shinystan(fitted_a_stan)
saveRDS(fitted,"models/fitted_a_model.stanModel")

mod=rstan::get_stanmodel(mod_a_neg_bin_bayes)

a=mod_a_neg_bin_bayes
summary(a)
resid(a,type="deviance")
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
