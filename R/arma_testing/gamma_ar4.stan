functions{
}

data{
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
  cov_matrix[N] Sigma;	

}

transformed data{
  vector[N] zeros; 
  zeros = rep_vector(0, N) ;
}

parameters{
  vector[N] mu;  
//  vector<lower=0>[N] precip_guess;   
  real<lower=0> shape ;
  real<lower=0> rate ;
  real<lower=0> sigma_error;
//  real<lower=0> sigma_transition; 
}

transformed parameters{
  vector[N] spi_hat ;

  for(i in 1:N) {
//  	spi_hat[i] = inv_Phi(gamma_cdf(precip_guess[i], shape, rate));
  	spi_hat[i] = inv_Phi(gamma_cdf(y[i], shape, rate));
  }
}

model{
shape ~ gamma(2, 0.1);
rate ~ gamma(2, 0.1);
 sigma_error ~ gamma(2, 0.1);
//sigma_transition ~ gamma(2, 0.1);

  mu ~ multi_normal(zeros, Sigma);


  y ~ gamma(shape, rate);
  mu ~ normal( spi_hat, sigma_error);

}





