functions {

  real qgamma(real p, real shape, real scale) {
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
  vector[N] y;     // list of y values from data
	
  int<lower=0> basis_dim;
  
  matrix[N,basis_dim] X; 
  matrix[basis_dim, basis_dim] S; 
    
  vector[2] b_0_mean_prior; 
  vector[2] b_0_scale_prior; 
  
  vector[basis_dim] b_mean_prior; 
  vector[basis_dim] b_scale_prior; 

  vector[2] lambda_mean_prior; 
  vector[2] lambda_scale_prior; 
 }
transformed data {  
 // vector[basis_dim] zero; 
 // zero = rep_vector(0, basis_dim) ;
}
parameters {
  real b_0_mean;
  real b_0_scale;  
  
  vector[basis_dim] b_mean;  
  vector[basis_dim] b_scale;   
  
  real<lower=0> lambda_mean ;
  real<lower=0> lambda_scale ;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_mean; 
  matrix[basis_dim, basis_dim] K_scale; 
  
  vector<lower=0>[N] mean_param;  
  vector<lower=0>[N] scale_param;
  
  K_mean = S * lambda_mean ;
  K_scale = S * lambda_scale ;
   
  mean_param = to_vector(X * b_mean) + b_0_mean;
  scale_param = to_vector(X * b_scale) + b_0_scale;
} 
model {
 
  lambda_mean ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
  lambda_scale ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
  //lambda_mean ~ gamma(0.05,0.005);
  //lambda_scale ~ gamma(0.05,0.005);
	
   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(b_mean_prior,K_mean); 
 
   b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
   b_scale ~ multi_normal_prec(b_scale_prior,K_scale); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
  y ~ gamma(mean_param ./ scale_param, rep_vector(1, N) ./ scale_param);  // shape is mean over scale, rate is 1 over scale
 
}
generated quantities {
 // vector[N] mean_est;
 // vector[N] scale_est;
 // vector[N] shape_est;   
 // vector[N] rate_est;  
 // real rho_mean;
 // real rho_scale;
  
 // mean_est = exp(mean_param);
 // scale_est = exp(scale_param);
  
 // shape_est = mean_param ./ exp(scale_param);
 // rate_est = rep_vector(1, N) ./ exp(scale_param);
  
 // rho_mean = log(lambda_mean);
 // rho_scale = log(lambda_scale);
}
