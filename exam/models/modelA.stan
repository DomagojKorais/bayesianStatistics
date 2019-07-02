//add exposure
//fit in cui vario l'intercetta sulle nazioni e il coeff + non centered parametrization since we had convergence problem
functions {
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real gamma_rate = gamma_rng(phi, phi / exp(eta));
    if (gamma_rate >= exp(20.79))
      return -9;
  return poisson_rng(gamma_rate);
}

}
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

parameters {//default is weakly informative
  real<lower=0> inv_phi;   
  real beta;
  real alpha; //the intercept
  vector[J] alpha_nat; //nations level
  
  //the deviation from the intercept at the different levels
  real dev_reg[L]; //deviation between the regions within a nation
  real dev_nat[J]; //deviation between the nations

  //the standard deviation for the deviations
  real<lower=0> sigma_reg;
  real<lower=0> sigma_nat;
}

transformed parameters {
    //varying intercepts

  real alpha_reg[L];//regions level
  
  real phi = inv(inv_phi);
  //real<lower=0> lambda[N];
  
  //compute the varying intercept at the nation level
    /*
  for(j in 1:J){
    alpha_nat[j] = alpha + dev_nat[j];
  }
  //compute varying intercept at the population within region level

  for(l in 1:L){
     alpha_reg[l] = alpha_nat[nationWithinRegion[l]] + dev_reg[l];
  }
*/

    
}
 

model {
  inv_phi ~ normal(0,1);
  alpha_nat ~ normal(4.5,1);
  beta ~ normal(0,1);
  //weakly informative prior on the standard deviation
  sigma_reg ~ cauchy(0,2.5);
  sigma_nat ~ cauchy(0,2.5);
  //distribution of the varying intercept
  dev_reg ~ normal(0,2);
  dev_nat ~ normal(0,2);

  //response
  deaths ~ neg_binomial_2_log( alpha_nat[nation_id] + beta * uvb + population, phi);
} 

generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  
  for (n in 1:N) {
    real eta_n = alpha_nat[nation_id[n]] + beta * uvb[n] + population[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(deaths[n]| eta_n, phi);

  }
}


