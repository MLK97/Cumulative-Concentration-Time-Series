data {
  int<lower=0> t;
  vector[t] x;
}

parameters {
  real<lower=0> tau; // characteristic time
  real<lower=0> sigma_x2; // fluctuation; x2 means sigma_x^2
  
  // x_2 - x_1 is the amplitude
  real<lower=0> x1;
  real<lower=x1> x2;
  
  real<lower=0, upper=t> t_switch; //time where tension graph starts increasing
}


model {
  real tsw;
  real dummy;
  // Priors
  tau ~ normal(700, 50);
  sigma_x2 ~ normal(10, 3);
  t_switch ~ uniform(0, t);
  x1 ~ normal(0.1, 1);
  x2 ~ normal(3, 2);
  
  tsw = ceil(t_switch);
  dummy = 531;
  

  x[1] ~ normal(0, sigma_x2);
  for (i in 2:t) {
  if(i < dummy)
      x[i] ~ normal((x[i-1]-x1)*exp(-1/tau)+x1, sqrt(sigma_x2*(1-exp(-(2/tau)))));
    else
      x[i] ~ normal((x[i-1]-x2)*exp(-1/tau)+x2, sqrt(sigma_x2*(1-exp(-(2/tau)))));
  }
}

