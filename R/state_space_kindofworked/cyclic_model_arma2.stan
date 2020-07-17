data {
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
	
  int<lower=0> basis_dim;
  
  matrix[N,basis_dim] X; 
  matrix[basis_dim, basis_dim] S; 
    
  vector[2] b_0_mean_prior; 
  vector[2] b_0_scale_prior; 
   vector[2] b_0_theta_prior; 
  
  vector[basis_dim] b_mean_prior; 
  vector[basis_dim] b_scale_prior; 
  vector[basis_dim] b_theta_prior; 

  vector[2] lambda_mean_prior; 
  vector[2] lambda_scale_prior; 
   vector[2] lambda_theta_prior;
 }
transformed data {  
  vector[basis_dim] zero; 
  zero = rep_vector(0, basis_dim) ;
}
parameters {
  real b_0_mean;
  real b_0_scale;  
  real b_0_theta;
  
  vector[basis_dim] b_mean;  
  vector[basis_dim] b_scale;   
 vector[basis_dim] b_theta; 
  
  real<lower=0> lambda_mean ;
  real<lower=0> lambda_scale ;
  real<lower=0> lambda_theta ;

//  real spi_origin;
  real<lower=0>  sigma_spi;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_mean; 
  matrix[basis_dim, basis_dim] K_scale; 
  matrix[basis_dim, basis_dim] K_theta; 
  
  vector[N] mean_param;  
  vector[N] scale_param;
  vector[N] theta_log;
  vector<lower=0, upper =1>[N] theta_param;
  
  K_mean = S * lambda_mean ;
  K_scale = S * lambda_scale ;
  K_theta = S * lambda_theta ;
   
  mean_param = to_vector(X * b_mean) + b_0_mean;
  scale_param = to_vector(X * b_scale) + b_0_scale;
  theta_log = to_vector(X * b_theta) + b_0_theta;
  theta_param = inv_logit(theta_log) ;
} 
model {
//   vector[N] p_est; 
 // vector[N] spi_est; 
//  vector[N] spi; 

   //Guess for Day 1
//   spi_origin ~ normal(0, 1); 

   //Gamma centered on 92 to produce a standard dev of 1 over sqrt 92
  // sigma_spi ~ gamma(50, 480);
 
  lambda_mean ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
  lambda_scale ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
   lambda_theta ~ gamma(lambda_theta_prior[1], lambda_theta_prior[2]);
  //lambda_mean ~ gamma(0.05,0.005);
  //lambda_scale ~ gamma(0.05,0.005);
	
   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(zero,K_mean); 
 
   b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
   b_scale ~ multi_normal_prec(zero,K_scale); 
  
   b_0_theta ~ normal(b_0_theta_prior[1],b_0_theta_prior[2]);   
   b_theta ~ multi_normal_prec(zero, K_theta); 

  // Set up day 1   
//	spi[1] = spi_origin;
//	p_est[1] = gamma_cdf(y[1], mean_param[1]/scale_param[1], 1 / scale_param[1]);
//	spi[1] = 	0 ;
//	p_est[1] = 0.5;
//	spi_est[1] = inv_Phi(p_est[1]);


  for(i in 2:N){
   if (y[i] == 0) {
      1 ~ bernoulli(theta_param[i]);
    //  p_est[i] = 0.5 * theta_param[i] ;
	}
    else {
      0 ~ bernoulli(theta_param[i]);
      y[i] ~ gamma(exp(mean_param[i]) / exp(scale_param[i]), inv(exp(scale_param[i])) );
    //	p_est[i] = gamma_cdf(y[i], exp(mean_param[i]) / exp(scale_param[i]), inv(exp(scale_param[i])) );
    }	
	//spi_est[i] = inv_Phi(p_est[i-1]);
	//spi[i] ~ normal(spi_est[i-1], sigma_spi);
  }  


//   for(i in 2:N) {
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
//  y[i] ~ gamma(mean_param[i] / scale_param[i], 1 / scale_param[i]);  // shape is mean over scale, rate is 1 over scale
//  p_est[i] = gamma_cdf(y[i], mean_param[i]/scale_param[i], 1 / scale_param[i]);
//  spi_est[i] = inv_Phi(p_est[i-1]);
 // spi[i] ~ normal(spi_est[i-1], sigma_spi);

//	}
 
}
generated quantities {
 // vector[N] mean_est;
 // vector[N] scale_est;
 // vector[N] shape_est;   
 // vector[N] rate_est;  
 // real rho_mean;
 // real rho_scale;
 //   real sd_r ;  
   vector[N] p_est; 
   vector[N] spi_est; 
//  vector[N] spi; 

  p_est[1] = 0.5;

  for(i in 2:N){
	if (y[i] == 0) {
  	  p_est[i] = 0.5 * (theta_param[i]) ;
	}
    else {
	  p_est[i] = gamma_cdf(y[i], exp(mean_param[i]) / exp(scale_param[i]), inv(exp(scale_param[i])) );
    }
   spi_est[i] = inv_Phi(p_est[i-1]);
  }

}
