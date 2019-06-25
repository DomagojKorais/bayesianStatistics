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

  int<lower=1> N;//number of records
  int<lower=0> deaths[N]; //number of deaths
  vector<lower=-10>[N] uvb;//uvb value
  vector[N] population;  //exposure
  // region-level data
  int<lower=1> J;//numero di gruppi (es n nazioni)
  int<lower=0> group_idx[N];

}
parameters {//default is weakly informative
  real<lower=0> inv_phi;   
  //intercette di mu e kappa
  real beta;   
  real alpha;
  //intercette e coeff raggruppati
  //riparametrizzazione
  vector[J] kappa_raw; 
    real<lower=0> sigma_kappa; 
  vector[J] mu_raw;  
    real<lower=0> sigma_mu; 
}
transformed parameters {
  real phi = inv(inv_phi);
  vector[J] mu = alpha  + sigma_mu * mu_raw;
  vector[J] kappa = beta + sigma_kappa * kappa_raw;
}
model {
deaths ~ neg_binomial_2_log(mu[group_idx] + kappa[group_idx] .* uvb + population , phi);
} 
generated quantities {
int y_rep[N];
vector[N] log_lik;
for (n in 1:N) {
real eta_n = mu[group_idx[n]] + kappa[group_idx[n]] * uvb[n]+ population[n];
y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
log_lik[n] = neg_binomial_2_log_lpmf(deaths[n]| eta_n, phi);

}
}


