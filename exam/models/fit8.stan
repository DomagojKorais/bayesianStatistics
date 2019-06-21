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
    int<lower=0> nation_idx[J];


  // building-level data
    int<lower=1> K;//number of columns model matrix (1+numero variabili)
    int<lower=1> J;//numero di gruppi (es n nazioni)
    //int<lower=1, upper=J> nation_idx[N];//indice della nazione
   
}
parameters {
  real<lower=0> inv_phi;   
  real beta;     
  vector[J] mu;  
  real<lower=0> sigma_mu; 
  real alpha;
  //vector[K] zeta;  

}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
    deaths ~ neg_binomial_2_log(mu[nation_idx] + beta * uvb , phi);
    mu ~ normal(alpha, sigma_mu);
    sigma_mu ~ normal(0, 1);
} 
generated quantities {
  int y_rep[N];
  for (n in 1:N) {
    real eta_n = mu[nation_idx[n]] + beta * uvb[n];
    y_rep[n] = neg_binomial_2_log_safe_rng(eta_n, phi);
  }
}
 

 