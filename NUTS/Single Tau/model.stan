data {
  int<lower=0> t; // total time of timeseries
  int N; // number of intervals
  vector[t] x; // tension points
  vector[N] alpha; // alpha value for dirichlet prior
}

parameters {
  real<lower=0> tau; // characteristic time
  real<lower=0> amplitude; // maximum amplitude between start and finish tension
  simplex[N] x_eq; // max tension-value of each t_switch note: simplex data type
  simplex[N] t_switch; // note: simplex data type
  real<lower=0> sigma_x; // fluctuation
}

transformed parameters {
  vector<lower=0>[N] t_sw; // actual t_switch output
  vector<lower=0>[N] x_hat;

  t_sw = cumulative_sum((t-1)*t_switch);
  x_hat = cumulative_sum(amplitude*x_eq);
}

model {
  int m;
  real t1;
  real t2;

  // Priors
  tau ~ normal(50, 250);
  sigma_x ~ normal(2, 10);
  amplitude ~ normal(9, 5);

  t_switch ~ dirichlet(alpha);
  x_eq ~ dirichlet(alpha);

  x[1] ~ normal(x_hat[1], sigma_x);
  m = 1;
  for (i in 2:t) {
    if(i - 1 <= t_sw[m]) {
      x[i] ~ normal((x[i - 1] - x_hat[m]) * exp(-1 / tau) + x_hat[m], sigma_x * sqrt(1 - exp(-2 / tau)));
    } else {
      t1 = t_sw[m] - (i - 2);
      t2 = (i - 1) - t_sw[m];
      x[i] ~ normal(((((x[i - 1]- x_hat[m]) * exp(-t1 / tau) + x_hat[m]) - x_hat[m + 1]) * exp(-t2 / tau)) + x_hat[m + 1], sigma_x * sqrt(1 - exp(-2 / tau)));
      m = m + 1;
    }
  }
  
}

// generated quantities {
//   vector[t] log_lik;
  
//   int m;
//   real t1;
//   real t2;

//   log_lik[1] = normal_lpdf(x[1] | x_hat[1], sigma_x);
//   m = 1;
//   for (i in 2:t) {
//     if(i - 1 <= t_sw[m]) {
//       log_lik[i] = normal_lpdf(x[i] | (x[i - 1] - x_hat[m]) * exp(-1 / tau) + x_hat[m], sigma_x * sqrt(1 - exp(-2 / tau)));
//     } else {
//       t1 = t_sw[m] - (i - 2);
//       t2 = (i - 1) - t_sw[m];
//       log_lik[i] = normal_lpdf(x[i] | ((((x[i - 1]- x_hat[m]) * exp(-t1 / tau) + x_hat[m]) - x_hat[m + 1]) * exp(-t2 / tau)) + x_hat[m + 1], sigma_x * sqrt(1 - exp(-2 / tau)));
//       m = m + 1;
//     }
//   }
// }
