#import libraries
library(tidyverse)
library(mlmRev)
library(GGally)
library(rstanarm)
library(bayesplot) 


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
compare_fits(fit1,"uvb",data=data_complete2)
compare_fits(fit2,"uvb + nation",data=data_complete2)
compare_fits(fit3,"uvb + nation + offset(local_male_population)",data=data_complete2)
compare_fits(fit4,"deaths ~ uvb + nation + region",data=data_complete2)
compare_fits(fit5,"deaths ~ uvb + nation + region + offset(local_male_population)",data=data_complete2)
compare_fits(fit6,"deaths ~ uvb + region",data=data_complete2)
compare_fits(fit7,"deaths ~ uvb + region + offset(local_male_population)",data=data_complete2)

