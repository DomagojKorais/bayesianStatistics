## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.align = 'center', warning=FALSE, message=FALSE, fig.asp=0.625, dev='png', global.par = TRUE, dev.args=list(pointsize=10), fig.path = 'figs/',  eval=TRUE)
library(MASS)

## ----setup, include=FALSE------------------------------------------------
library(knitr)
local({
    hook_plot = knit_hooks$get('plot')
    knit_hooks$set(plot = function(x, options) {
        paste0('\n\n----\n\n', hook_plot(x, options))
    })
})


## ----eval=TRUE,echo=TRUE, warning=FALSE, results='hide',message=FALSE----
library(rstan)
library(dplyr)
library(lubridate)
library(ggplot2)
library(bayesplot)

theme_set(bayesplot::theme_default())

# seed for R's pseudo-RNGs, not Stan's
set.seed(1123) 


## ----eval=TRUE,echo=TRUE-------------------------------------------------
pest_data <- readRDS('pest_data.RDS')
str(pest_data)
summary(pest_data)

##explore the data
ggplot(data = pest_data,aes(x=complaints , y=traps))+
  geom_point()
library(GGally)
ggpairs(pest_data[c(14,5)], aes( alpha = 0.4))
##number of buildings
N_buildings <- length(unique(pest_data$building_id))
N_buildings


## ----eval=TRUE,echo=FALSE------------------------------------------------
ggplot(pest_data, aes(x = complaints)) + 
  geom_bar()

ggplot(pest_data, aes(x = traps, y = complaints, color = live_in_super == TRUE)) + 
    geom_jitter()

## ---- data-plots-ts, fig.width = 6, fig.height = 8-----------------------
ggplot(pest_data, aes(x = date, y = complaints, color = live_in_super == TRUE)) + 
  geom_line(aes(linetype = "Number of complaints")) + 
  geom_point(color = "black") + 
  geom_line(aes(y = traps, linetype = "Number of traps"), color = "black", size = 0.25) + 
  facet_wrap(~building_id, scales = "free", ncol = 2, labeller = label_both) + 
  scale_x_date(name = "Month", date_labels = "%b") + 
  scale_y_continuous(name = "", limits = range(pest_data$complaints)) + 
  scale_linetype_discrete(name = "") + 
  scale_color_discrete(name = "Live-in super")


## ------------------------------------------------------------------------
## compile the model
comp_model_P <- stan_model('simple_poisson_regression2.stan')

## ----stan-data-----------------------------------------------------------
## arrange data into a list (data has to be in a list to pass them to Stan)
stan_dat_simple <- list(
  N = nrow(pest_data), 
  complaints = pest_data$complaints,
  traps = pest_data$traps
)
str(stan_dat_simple)

## ----fit_P_real_data, cache=TRUE-----------------------------------------
## fit the model
## Nel warmup stan esplora lo stpazio dei parametri
##Nel metropolis o nel Gibbs non c'è questa differenziazione
fit_P_real_data <- sampling(comp_model_P, data = stan_dat_simple)


## ----results_simple_P----------------------------------------------------
## print the parameters 
print(fit_P_real_data, pars = c('alpha','beta'))
##n_eff è un valore che indica quanto è buono in termini di indipendenza la simulazione

## ----hist_simple_P-------------------------------------------------------
mcmc_hist(as.matrix(fit_P_real_data, pars = c('alpha','beta')))
mcmc_scatter(as.matrix(fit_P_real_data, pars = c('alpha','beta')), alpha = 0.2)
# il fatto che alpha e beta siano correlate è la prova del fatto che complaints e traps sono correlate, perchè?

## ------------------------------------------------------------------------
## posterior predictive checking
y_rep <- as.matrix(fit_P_real_data, pars = "y_rep")
ppc_dens_overlay(y = stan_dat_simple$complaints, y_rep[1:200,])
#da qui vediamo che il modello simulato presenta una varianza inferiore alla reale

## ------------------------------------------------------------------------
## standardised residuals of the observed vs predicted number of complaints
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$complaints - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)
#residui per lo più positivi=> modello sottostima

## ------------------------------------------------------------------------
ggplot(pest_data, aes(x = log(total_sq_foot), y = log1p(complaints))) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)


## ------------------------------------------------------------------------
## add the two variables to the list of the data
stan_dat_simple$area <- log(pest_data$total_sq_foot/1e4)
stan_dat_simple$super <- pest_data$live_in_super
str(stan_dat_simple)


## ----compmultPDGP--------------------------------------------------------
## compile the model
comp_model_P_mult <- stan_model('multiple_poisson_regression.stan')


## ----fit_mult_P_real_dat-------------------------------------------------
fit_model_P_mult_real <- sampling(comp_model_P_mult, data = stan_dat_simple)
y_rep <- as.matrix(fit_model_P_mult_real, pars = "y_rep")
ppc_dens_overlay(stan_dat_simple$complaints, y_rep[1:200,])
#più o meno abbiamo gli stessi problemi di prima(dobbiamo rilassare l'ipotesi sulla varianza = media)
## ------------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_simple$complaints, 
  yrep = y_rep,
  x = stan_dat_simple$traps
) + 
  labs(x = "Number of traps", y = "Number of complaints")


## ---- cache=TRUE, results="hide", message=FALSE--------------------------
comp_model_NB <- stan_model('multiple_NB_regression.stan')


## ----runNB---------------------------------------------------------------
fitted_model_NB <- sampling(comp_model_NB, data = stan_dat_simple)
samps_NB <- rstan::extract(fitted_model_NB)


## ----ppc-full------------------------------------------------------------
## predictions vs. the data
y_rep <- samps_NB$y_rep
ppc_dens_overlay(stan_dat_simple$complaints, y_rep[1:200,])


## ------------------------------------------------------------------------
## standardisd residuals
mean_inv_phi <- mean(samps_NB$inv_phi)
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$complaints - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)


## ------------------------------------------------------------------------
## prediction by number of traps
ppc_intervals(
  y = stan_dat_simple$complaints, 
  yrep = y_rep,
  x = stan_dat_simple$traps
) + 
  labs(x = "Number of traps", y = "Number of complaints")


## ----ppc-group_means-----------------------------------------------------
ppc_stat_grouped(
  y = stan_dat_simple$complaints, 
  yrep = y_rep, 
  group = pest_data$building_id, 
  stat = 'mean',
  binwidth = 0.2
)

