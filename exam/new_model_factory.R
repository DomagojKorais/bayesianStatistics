#well made model factory
library(tidyverse)
library(mlmRev)
library(GGally)
library(rstanarm)
library(bayesplot) 
library(loo)
#import and arrange data
#prepare data
seed=pi
data = Mmmec
populationbycountry19802010millions <- read_csv("populationbycountry19802010millions.csv", 
                                                col_types = cols(`1980` = col_double()))
#rename Western Germany to make in compatible with our data
filter(populationbycountry19802010millions,X1 =="United Kingdom")
populationbycountry19802010millions = populationbycountry19802010millions%>%mutate(X1 =replace(X1,X1=="Germany, West","W.Germany"))
populationbycountry19802010millions = populationbycountry19802010millions%>%mutate(X1 =replace(X1,X1=="United Kingdom","UK"))
filter(populationbycountry19802010millions,X1 =="UK")
filter(data,nation =="UK")
pop = populationbycountry19802010millions%>%select(nation=X1,population = '1980')

nation = data%>%group_by(nation)%>%summarise(tot_expected=sum(expected))
nation_complete = inner_join(pop,nation,by=c("nation"="nation"))
male_nation_complete = nation_complete%>%mutate(male_population=population/2)%>%select(-population)
data_complete2 = inner_join(male_nation_complete,data,by=c("nation"="nation"))
data_complete2 = data_complete2%>%mutate(local_male_population=(male_population*expected)/tot_expected)
data_complete2 = data_complete2%>%mutate(local_rate = (deaths)/(local_male_population))

N_nations = length(unique(data_complete2$region))
#prepare data for stan
nation_idx=seq(1,N_nations)
data3=data_complete2%>%mutate(
  nation_fac = factor(nation, levels = unique(nation)),
  nation_idx = as.integer(nation_fac)
)
#prepare group
group_data = data.frame(data3)%>%
  group_by(nation_idx)%>%
  summarise(uvb=mean(uvb))%>%
  select(uvb)
#prepare data
stan_dat_hier =
  with(data3,
       list(deaths = deaths,
            uvb = uvb,
            N = length(deaths),
            J = N_nations,
            K = 1,
            group_idx = nation_idx,
            population = log(local_male_population)
       )
  )
#import model
model = rstan::stan_model("models/fit16.stan")
#run model
fitted <- rstan::sampling(model, data = stan_dat_hier,
                             chains = 2, cores = 4, iter = 4000, verbose=TRUE,seed=seed)
saveRDS(fitted,"models/fit17.stanModel")

color_scheme_set("darkgray")
np_cp <- nuts_params(fitted)
head(np_cp)
posterior <- as.array(fitted)
mcmc_parcoord(posterior, np = np_cp,regex_pars = c("mu"))

