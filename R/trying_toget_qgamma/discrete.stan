data {
  int<lower=0> N;  //number of points
  int<lower=0> n_day;  //number of days
  vector[N] y;     // list of y values from data
  int day[N];   // list of months from data
}
parameters {
  vector<lower=0>[n_day] scale_param ;
  vector<lower=0>[n_day] shape_param ;
  vector<lower=0, upper=1>[n_day] theta ;
}
model {
 shape_param ~ gamma(0.1, 0.02);
 scale_param ~ gamma(0.1, 0.02);
 theta ~ beta(2,2);

  for(i in 1:N){
   if (y[i] == 0)
      1 ~ bernoulli(theta[day[i]]);
    else {
      0 ~ bernoulli(theta[day[i]]);
      y[i] ~ gamma(shape_param[day[i]], 1/scale_param[day[i]]);
    }	
  }   
}
generated quantities {
  
}






