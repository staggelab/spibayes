functions {

vector qgamma(vector y,        // unknowns
              vector theta,    // parameters
              real[] x_r,      // data (real)
              int[] x_i) {     // data (integer)
  vector[1] z;
  z[1] = gamma_cdf(y, theta[1], 1/theta[2]) - theta[3];
  return z;
}

real qgamma_approx(real p, real shape, real scale) {
    real skew;
    real true_mean;
    real true_var;
    real w_variate; 
    real q_val;

    skew = 2/((shape^0.5));

    w_variate = (2/skew)*((1-(skew/6)^2 + (skew/6)*inv_Phi(p))^3-1);

    true_mean = shape .* scale ; 
    true_var = shape .* (scale^2) ; 

    q_val = true_mean + (true_var ^ 0.5) * w_variate;
	if (q_val < 0)
	    q_val = 0;

    return q_val;
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
  vector[N] prob_previous;    
  vector[N] spi_previous;   

prob_previous[1] = 0.5;
 for(i in 1:(N-1)){
   if (y[i+1] == 0) 
	  prob_previous[i] = theta_param[day[i]] / 2;
   else {
	  prob_previous[i] = (1-theta_param[day[i]])*(gamma_cdf(y[i], shape_param[day[i]], 1/scale_param[day[i]])) + theta_param[day[i]];
    }
  }

spi_previous = inv_Phi(prob_previous);

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


}
