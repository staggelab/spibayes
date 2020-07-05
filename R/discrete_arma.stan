functions { 

//vector qgamma(vector y,        // unknowns
 //             vector theta,    // parameters
 //             real[] x_r,      // data (real)
 //             int[] x_i) {     // data (integer)
//  vector[1] z;
//  z[1] = gamma_cdf(y, theta[1], 1/theta[2]) - x_r[1];
//  return z;
}

}
data {
  int<lower=0> N;  //number of points
  int<lower=0> n_day;  //number of days
  vector[N] y;     // list of y values from data
  int day[N];   // list of months from data
}
parameters {
  vector<lower=0>[n_day] scale[n_day];
  vector<lower=0>[n_day] shape[n_day];
  vector<lower=0, upper=1>[n_day] theta;
}
model {
  shape ~ gamma(0.1, 0.02);
  scale ~ gamma(0.1, 0.02);
  for(i in 1:N){
   if (y[i] == 0)
      1 ~ bernoulli(theta[day[i]]);
    else {
      0 ~ bernoulli(theta[day[i]]);
      y[i] ~ gamma(shape[day[i]], 1/scale[day[i]]);
    }	
  }   
}
generated quantities {
  
	


}

