functions {
  /*
    * Alternative to poisson_log_rng() that 
  * avoids potential numerical problems during warmup
  */
    int poisson_log_safe_rng(real eta) {
      real pois_rate = exp(eta);
      if (pois_rate >= exp(20.79))
        return -9;
      return poisson_rng(pois_rate);
    }
}
data {
  //input data
  int<lower=1> N;
  int<lower=0> complaints[N]; //another syntax is with vector
  vector<lower=0>[N] traps;
  vector<lower=0>[N] super;
  vector<lower=0>[N] area;
  
}
parameters {
  //parameters of the model
  //in this case alpha and beta
  real alpha;
  real beta;
  real betaSuper;
}
model {
  //model
  //priors and likelihood
  //se non mettiamo niente come prior prenderà una uniforme
  beta ~ normal(-0.25,1);
  
  alpha ~ normal(log(4),1);
  complaints ~ poisson_log(alpha + beta * traps + betaSuper * super + area);
} 
generated quantities {
  // sample predicted values from the model for posterior predictive checks
  int y_rep[N];
  
  for (n in 1:N) {
    real eta_n = alpha + beta * traps[n] + betaSuper * super[n] + area [n];
    y_rep[n] = poisson_log_safe_rng(eta_n); //se il parametro è troppo grande lo riportiamo a un generico valore perchè c'è un bug in poisson_log_rng 
  }
}
//riga da lasciare vuota
