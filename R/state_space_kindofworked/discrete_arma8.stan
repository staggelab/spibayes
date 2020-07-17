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
    real skew_six;
    real norm_dev;
    real true_mean;
    real true_var;
    real w_variate; 
    real q_val;

    skew = 2/((shape^0.5));
	skew_six = skew/6 ;
    norm_dev = inv_Phi(p);
    w_variate = (2/skew)*((1-(skew_six)^2 + (skew_six)*norm_dev)^3-1);

    true_mean = shape * scale ; 
    true_var = shape * (scale^2) ; 

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

  int x_i[0] ;
  real x_r[0]; 
}
parameters {
  vector[n_day] scale_param ;
  vector[n_day] shape_param ;
  vector<lower=0, upper=1>[n_day] theta_param ;
  real<lower=0> sigma;
	

}
transformed parameters {

  vector<lower=0, upper=1>[N] prob_prev_total;    
  vector[N] spi_previous;  
  vector<lower=0>[N] mu;   
  vector[N] prob_current_total;
  vector[N] prob_current_pos;
  vector[3] opt_param ;
  real y_guess_real ; 
  vector[1]  y_guess ;
  vector[N]  y_guess_save ;
  vector[1]  shaba;

  vector[n_day] scale_trans ;
  vector[n_day] shape_trans ;

	scale_trans = exp(scale_param) ;
	shape_trans = exp(shape_param) ;

 y_guess[1] = 5;

prob_prev_total[1] = 0.5;

 for(i in 2:N){
   if (y[i-1] == 0) {
	  prob_prev_total[i] = theta_param[day[i-1]] / 2;
   }
   else {
	  prob_prev_total[i] = (1-theta_param[day[i-1]])*(gamma_cdf(y[i-1], shape_trans[day[i-1]], 1/scale_trans[day[i-1]])) + theta_param[day[i-1]];
    }
  }

 spi_previous = inv_Phi(prob_prev_total);
// This is where you would add some movement to the SPI For now assume its random so stationary over long term
 prob_current_total = Phi(spi_previous);

 for(i in 1:N){
    prob_current_pos[i] = (1/(1-theta_param[day[i]])) * (prob_current_total[i] - theta_param[day[i]]); 
   
    if (prob_current_pos[i] > 0) { 

      y_guess_real = qgamma_approx(prob_current_pos[i], shape_trans[day[i]], scale_trans[day[i]]);
	  y_guess[1] = y_guess_real;
	  y_guess_save[i] = y_guess_real;

	  opt_param = [shape_trans[day[i]], scale_trans[day[i]], prob_current_pos[i]]';
	  shaba = algebra_solver(qgamma, y_guess, opt_param, x_r, x_i) ;
	   if (shaba[1] < 0) {
			shaba[1] = 0;
		}
	  mu[i] = shaba[1];
   }
   else {
     prob_current_pos[i] = theta_param[day[i]]/2;
     mu[i] = 0;
  }
}

}
model {
 shape_param ~ normal(1, 0.4);
 scale_param ~ normal(0.4, 0.4);
 theta_param ~ beta(2,2);
 sigma ~ exponential(1);

  for(i in 1:N){
   if (y[i] == 0)
      1 ~ bernoulli(theta_param[day[i]]);
    else {
      0 ~ bernoulli(theta_param[day[i]]);
  // y[i] ~ gamma(shape_param[day[i]], 1/scale_param[day[i]]); 
     y[i] ~ normal(mu[i], sigma);

    }	
  }   
}
generated quantities {


}
