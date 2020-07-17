functions {
  int fib(int n);
}
data {
  int<lower=0> N;  //number of points
  int<lower=0> n_day;  //number of days
  vector[N] y;     // list of y values from data
  int day[N];   // list of months from data
}
transformed data {

}
parameters {
  vector[n_day] scale_param ;
  vector[n_day] shape_param ;
  vector<lower=0, upper=1>[n_day] theta_param ;
}
transformed parameters {

}
model {
 shape_param ~ gamma(0.1, 0.02);
 scale_param ~ gamma(0.1, 0.02);
 theta_param ~ beta(2,2);

  for(i in 1:N){
   if (y[i] == 0)
      1 ~ bernoulli(theta_param[day[i]]);
    else {
      0 ~ bernoulli(theta_param[day[i]]);
      y[i] ~ gamma(shape_param[day[i]], 1/scale_param[day[i]]);
    }	
  }   
}
generated quantities {
int fib_result;

fib_result = fib(35);

}

