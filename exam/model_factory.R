#import libraries
library(tidyverse)
library(mlmRev)
library(GGally)
library(rstanarm)
library(bayesplot) 
library(loo)


#prepare data
data = Mmmec
populationbycountry19802010millions <- read_csv("populationbycountry19802010millions.csv", 
                                                col_types = cols(`1980` = col_double()))
#rename Western Germany to make in compatible with our data
filter(populationbycountry19802010millions,X1 =="Germany, West")
populationbycountry19802010millions = populationbycountry19802010millions%>%mutate(X1 =replace(X1,X1=="Germany, West","W.Germany"))


nation = data%>%group_by(nation)%>%summarise(tot_expected=sum(expected))
nation_complete = inner_join(pop,nation,by=c("nation"="nation"))
male_nation_complete = nation_complete%>%mutate(male_population=population/2)%>%select(-population)
data_complete2 = inner_join(male_nation_complete,data,by=c("nation"="nation"))
data_complete2 = data_complete2%>%mutate(local_male_population=(male_population*expected)/tot_expected)
data_complete2 = data_complete2%>%mutate(local_rate = (deaths)/(local_male_population))


#fuck the models

fit1 <- stan_glm(deaths ~ uvb , data=data_complete2, family=neg_binomial_2(),cores=2)
fit2 <- stan_glm(deaths ~ uvb + nation , data=data_complete2, family=neg_binomial_2(),cores=4)
fit3 <- stan_glm(deaths ~ uvb + nation + offset(local_male_population) , data=data_complete2, family=neg_binomial_2(),cores=4)
fit4 <- stan_glm(deaths ~ uvb + nation + region, data=data_complete2, family=neg_binomial_2(),cores=4)
fit5 <- stan_glm(deaths ~ uvb + nation + region  + offset(local_male_population), data=data_complete2, family=neg_binomial_2(),cores=4)
fit6 <- stan_glm(deaths ~ uvb + region , data=data_complete2, family=neg_binomial_2(),cores=4)
fit7 <- stan_glm(deaths ~ uvb + region + offset(local_male_population), data=data_complete2, family=neg_binomial_2(),cores=4)
saveRDS(fit1,"models/fit1.stanModel")
saveRDS(fit2,"models/fit2.stanModel")
saveRDS(fit3,"models/fit3.stanModel")
saveRDS(fit4,"models/fit4.stanModel")
saveRDS(fit5,"models/fit5.stanModel")
saveRDS(fit6,"models/fit6.stanModel")
saveRDS(fit7,"models/fit7.stanModel")
compare_fits = function (fit,title,data){
  posterior <- as.array(fit)
  yrep_poisson <- posterior_predict(fit, draws = 500)
  y=data$deaths
  mean_y_rep <- colMeans(yrep_poisson)
  std_resid <- (y - mean_y_rep) / sqrt(mean_y_rep)
  qplot(mean_y_rep, std_resid,main=title) + hline_at(2) + hline_at(-2)+
    
    labs(x="Mean of y_rep", y= "Stand. residuals")+
    xaxis_text(on =TRUE, size=22)+
    yaxis_text(on =TRUE, size=22)+
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16))
  }
compare_fits(fit1,"1 uvb",data=data_complete2)
compare_fits(fit2,"2 uvb + nation",data=data_complete2)
compare_fits(fit3,"3 uvb + nation + offset(local_male_population)",data=data_complete2)
compare_fits(fit4,"4 deaths ~ uvb + nation + region",data=data_complete2)
compare_fits(fit5,"5 deaths ~ uvb + nation + region + offset(local_male_population)",data=data_complete2)
compare_fits(fit6,"6 deaths ~ uvb + region",data=data_complete2)
compare_fits(fit7,"7 deaths ~ uvb + region + offset(local_male_population)",data=data_complete2)

#MODEL COMPARISON
#loo
loo_1  = loo(fit1)
loo_2  = loo(fit2)
loo_3  = loo(fit3)
loo_4  = loo(fit4)
loo_5  = loo(fit5)
loo_6  = loo(fit6)
loo_7  = loo(fit7)
loo_diff = compare(loo_1,loo_2,loo_3,loo_4,loo_5,loo_6,loo_7)
loo_diff
#waic

waic_1  = waic(fit1)
waic_2  = waic(fit2)
waic_3  = waic(fit3)
waic_4  = waic(fit4)
waic_5  = waic(fit5)
waic_6  = waic(fit6)
waic_7  = waic(fit7)
waic_diff = compare(waic_1,waic_2,waic_3,waic_4,waic_5,waic_6,waic_7)
waic_diff

#hierarchical models! =) let's check if we can do better, hierarchy is nation>region
comp_model_NB_hier <- stan_model('hierarchical_NB_regression.stan')


N_nations=8 #DEVI TRASFORMARE IL FACTOR NAZIONE IN INTEGER
#nation_data=
nation_idx=seq(1,N_nations)
data3=data_complete2%>%mutate(
  nation_fac = factor(nation, levels = unique(nation)),
  nation_idx = as.integer(nation_fac)
)
  #prepare data
stan_dat_hier =
  with(data3,
       list(deaths = deaths,
            uvb = uvb,
            N = length(deaths),
            J = N_nations,
            K = 2,
            nation_idx = nation_idx
       )
  )

#prepare data
stan_dat_hier_9 =
  with(data3,
       list(deaths = deaths,
            uvb = uvb,
            N = length(deaths),
            J = length(unique(data3$nation)),
            K = 2,
            nation_idx = nation_idx,
            population = log(local_male_population)
       )
  )
#RUN MODEL
model9 = stan_model("fit9.stan")S

fitted_9 <- sampling(model9, data = stan_dat_hier_9,
                     chains = 4, cores = 4, iter = 4000)
samps_hier_9 <- rstan::extract(fitted_9)
saveRDS(fitted_9,"models/fit9.stanModel")

print(fitted_9, pars = c('sigma_mu','beta','alpha','phi','mu'))
plot(fitted_9)



y_rep <- as.matrix(fitted_9, pars = "y_rep")
y=data_complete2$deaths
ppc_stat(y,y_rep,stat="sd")

ppc_stat_grouped(y,y_rep,group=data_complete$nation,stat="mean")+
  ggtitle("Using exposures means comparison")



#y_rep <- posterior_predict(fit3, draws = 500)
ppc_stat_grouped(y,y_rep,group=data_complete$nation,stat="sd")+
  ggtitle("Using exposures sd comparison")

mean_y_rep <- colMeans(y_rep)

mean_inv_phi <- mean(as.matrix(fitted_9, pars = "inv_phi"))
std_resid <- (stan_dat_hier$deaths - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)


#RUN MODEL 10
#prepare data
stan_dat_hier_10 =
  with(data3,
       list(deaths = deaths,
            uvb = uvb,
            N = length(deaths),
            J = length(unique(region)),
            K = 2,
            #region_idx = as.integer(unique(region)),
            region_fac = factor(region, levels = unique(region)),
            region_idx = as.integer(factor(region, levels = unique(region))),
            population = log(local_male_population)
       )
  )
model10 = stan_model("fit10.stan")

fitted_10 <- sampling(model10, data = stan_dat_hier_10,
                     chains = 1, cores = 4, iter = 4000, verbose=TRUE)
samps_hier_10 <- rstan::extract(fitted_10)
saveRDS(fitted_10,"models/fit10.stanModel")

print(fitted_10, pars = c('sigma_mu','beta','alpha','phi','mu'))
plot(fitted_10)



y_rep <- as.matrix(fitted_10, pars = "y_rep")
y=data_complete2$deaths
ppc_stat(y,y_rep,stat="sd")

ppc_stat_grouped(y,y_rep,group=data_complete$nation,stat="mean")+
  ggtitle("Using exposures means comparison")



#y_rep <- posterior_predict(fit3, draws = 500)
ppc_stat_grouped(y,y_rep,group=data_complete$nation,stat="sd")+
  ggtitle("Using exposures sd comparison")

mean_y_rep <- colMeans(y_rep)

mean_inv_phi <- mean(as.matrix(fitted_10, pars = "inv_phi"))
std_resid <- (stan_dat_hier$deaths - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)



