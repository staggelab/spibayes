functions{
}

data{
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
  int<lower=2> n_roll ;   //accumulation period
}

transformed data{
  int N_withna;
  int n_roll_count ;

  n_roll_count = n_roll - 1 ;
  N_withna = N + n_roll_count;
}

parameters{
  vector[N] innov;  
  real<lower=0> shape ;
  real<lower=0> rate ;
  real<lower=0> sigma_error;
}

transformed parameters{
  vector[N] spi_hat ;
  vector[N] mu ;

  for(i in 1:N) {
  	spi_hat[i] = inv_Phi(gamma_cdf(y[i], shape, rate));
  }

 for(i in n_roll:3652) {
	mu[i] = sum(innov[(i-n_roll_count):i]);
  }


}

model{
shape ~ gamma(2, 0.1);
rate ~ gamma(2, 0.1);
sigma_error ~ gamma(2, 0.1);
innov ~ normal(0, 0.1042572);

  y ~ gamma(shape, rate);

  for(i in n_roll:3652) {
    mu[i] ~ normal(spi_hat[i-n_roll_count], sigma_error);
  }

 // for(i in 2:N) {
 //    mu[i] ~ normal(spi_hat[i-1], 0.104);
 // }
}
