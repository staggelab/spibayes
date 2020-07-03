functions {
vector qgamma(vector y,        // unknowns
              vector theta,    // parameters
              real[] x_r,      // data (real)
              int[] x_i) {     // data (integer)
  vector[1] z;
  z[1] = gamma_cdf(y, theta[1], 1/theta[2]) - x_r[1];
  return z;
}
}
data {
  int<lower=0> N;  //number of points
  int<lower=0> n_day;  //number of days
  vector[N] y;     // list of y values from data
  int day[N];   // list of months from data
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
vector[1] qgamma_result;
vector[2] theta = [shape_param[1], scale_param[1]]';
vector[1] y_guess = [0.5]';

int x_i[0];
real x_r[1]; // Return the 40th percentile

x_r[1] = 0.4;

qgamma_result = algebra_solver(qgamma, y_guess, theta, x_r, x_i);
}

