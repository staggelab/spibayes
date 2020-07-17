functions{
}

data{
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
}

transformed data{
}

parameters{
  vector[N] mu;  
  real<lower=0> shape ;
  real<lower=0> rate ;
  real<lower=0> sigma_error;
//  real<lower=0> sigma_transition; 
}

transformed parameters{
  vector[N] spi_hat ;

  for(i in 1:N) {
  	spi_hat[i] = inv_Phi(gamma_cdf(y[i], shape, rate));
  }
}

model{
shape ~ gamma(2, 0.1);
rate ~ gamma(2, 0.1);
sigma_error ~ gamma(2, 0.1);
//sigma_transition ~ gamma(2, 0.1);

  for(i in 2:N) {
     mu[i] ~ normal(mu[i-1], 0.1474401);
  }

  y ~ gamma(shape, rate);

  mu ~ normal(spi_hat, sigma_error);

}
