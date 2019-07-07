data {

int<lower=1> N;//number of records (level 1 counties)
int<lower=1> L;//numero di regioni (level 2 groups)
int<lower=1> J;//numero di nazioni (level 3 groups)

int<lower=0> deaths[N]; //number of deaths
 vector<lower=-10>[N] uvb;//uvb value
 vector[N] population;  //exposure

int<lower=0> region_id[N]; //region
int<lower=0> nation_id[N]; //nation
int<lower=1,upper=J> nationWithinRegion[L]; 

}

parameters {
real beta;
real<lower=0>  alpha; //the intercept


//the deviation from the intercept at the different levels
real dev_reg[L]; //deviation between the regions within a nation
real dev_nat[J]; //deviation between the nations

//the standard deviation for the deviations
real<lower=0> sigma_reg;
real<lower=0> sigma_nat;

}

transformed parameters {
//varying intercepts
  vector[J] alpha_nat; //nations level
  vector[L] alpha_reg;//regions level
  real var_reg = sigma_reg^2;
  real var_nat = sigma_nat^2;


//real<lower=0> lambda[N];

//compute the varying intercept at the nation level
for(j in 1:J){
alpha_nat[j] = alpha + dev_nat[j];
}

//compute varying intercept at the population within region level

for(l in 1:L){
alpha_reg[l] = alpha_nat[nationWithinRegion[l]] + dev_reg[l];
}



}


model {

alpha ~ gamma(3,1);
beta ~ normal(0,1);
//weakly informative prior on the standard deviation
sigma_reg ~ gamma(9,20);
sigma_nat ~ gamma(9,20);
//distribution of the varying intercept
dev_reg ~ normal(0,sigma_reg);
dev_nat ~ normal(0,sigma_nat);

//response
deaths ~ poisson( alpha_reg[region_id] + beta * uvb + population);
} 

generated quantities {
int y_rep[N];
real log_lik[N];

for (n in 1:N) {
real eta_n = alpha_reg[region_id[n]] + beta * uvb[n] + population[n];
y_rep[n] = poisson_rng(eta_n);
log_lik[n] = poisson_log_lpmf(deaths[n]| eta_n);

}
}


