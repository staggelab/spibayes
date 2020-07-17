functions{
}

data{
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
}

transformed data{
}

parameters{
  real<lower=0> shape ;
  real<lower=0> mean_param ;
}

transformed parameters{
}

model{
shape ~ gamma(2, 0.1);
mean_param ~ gamma(2, 0.1);

  y ~ gamma(shape, shape/mean_param);

}
