//add exposure
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
  
  vector[N] population;  
  
  
  // building-level data
  int<lower=1> K;//number of columns model matrix (1+numero variabili)
  int<lower=1> J;//numero di gruppi (es n nazioni)
  matrix[J,K] group_data;
  int<lower=0> region_idx[N];
  //int<lower=1, upper=J> nation_idx[N];//indice della nazione
  
}
parameters {
  real<lower=0> inv_phi;   
  real beta;     
  real<lower=0> sigma_mu; 
  real alpha;
  vector[J] mu_raw;        
  vector[K] zeta;          
  vector[J] kappa_raw;       
  real<lower=0> sigma_kappa; 
  vector[K] gamma;   
  
}
transformed parameters {
  real phi = inv(inv_phi);
  vector[J] mu = alpha + group_data * zeta + sigma_mu * mu_raw;
  vector[J] kappa = beta + group_data * gamma + sigma_kappa * kappa_raw;
}
model {
  //deaths ~ neg_binomial_2_log(mu[nation_idx] + beta * uvb + population , phi);
  //mu ~ normal(alpha, sigma_mu);
  //sigma_mu ~ normal(0, 1);
  
  deaths ~ neg_binomial_2_log(
    mu[region_idx] + kappa[region_idx] .* uvb + population ,
    phi
  );
  
} 
generated quantities {
  int y_rep[N];
  vector[N] log_lik;
  for (n in 1:N) {
    real eta_n = mu[region_idx[n]] + kappa[region_idx[n]] * uvb[n]+ population[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
    log_lik[n] = neg_binomial_2_log_lpmf(deaths[n]| eta_n, phi);
  }
}


