
#' Cyclic Stan model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @export
cyclic_model <- "
data {
  int<lower=0> N;  //number of points
  int<lower=0> N_pos;  //number of days with positive precip

  int<lower=0> basis_dim;

  vector[N_pos] y_pos;     // list of y values with positive precip
  int<lower=0,upper=1> y_zero[N];     // list of one or zero for positive precip
 
  matrix[N_pos,basis_dim] X_pos;  // Basis for positive precip days
  matrix[N,basis_dim] X;   // Basis for all days

  matrix[basis_dim, basis_dim] S; 
    
  vector[2] b_0_shape_prior; 
  vector[2] b_0_rate_prior; 
  vector[2] b_0_theta_prior; 

  real lambda_shape_prior; 
  real lambda_rate_prior; 
  real lambda_theta_prior; 
 }
transformed data {  
  vector[basis_dim] zeros; 
  vector[N_pos] ones_pos; 

  zeros = rep_vector(0, basis_dim) ;
  ones_pos = rep_vector(1, N_pos) ;
}
parameters {
  real b_0_shape;
  real b_0_rate;  
  real b_0_theta;  
 
  vector[basis_dim] b_shape;  
  vector[basis_dim] b_rate;   
  vector[basis_dim] b_theta; 
 
  real<lower=0> lambda_shape ;
  real<lower=0> lambda_rate ;
  real<lower=0> lambda_theta ;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_shape; 
  matrix[basis_dim, basis_dim] K_rate; 
  matrix[basis_dim, basis_dim] K_theta; 
 
  vector[N_pos] shape_param;  
  vector[N_pos] rate_param;


  K_shape = S * lambda_shape ;
  K_rate = S * lambda_rate ;
  K_theta = S * lambda_theta ;
  
  shape_param = exp(X_pos * b_shape + b_0_shape);
  rate_param = exp(X_pos * b_rate + b_0_rate);

} 
model {
 // lambda_mean ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
 // lambda_scale ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
 // lambda_theta ~ gamma(lambda_theta_prior[1], lambda_theta_prior[2]);
	
  lambda_shape ~ cauchy(0, lambda_shape_prior);
  lambda_rate ~ cauchy(0, lambda_rate_prior);
  lambda_theta ~ cauchy(0, lambda_theta_prior);

   b_0_shape  ~ normal(b_0_shape_prior[1],b_0_shape_prior[2]);   
   b_shape ~ multi_normal_prec(zeros,K_shape); 
 
   b_0_rate  ~ normal(b_0_rate_prior[1],b_0_rate_prior[2]);  
   b_rate ~ multi_normal_prec(zeros,K_rate); 

   b_0_theta ~ normal(b_0_theta_prior[1],b_0_theta_prior[2]);   
   b_theta ~ multi_normal_prec(zeros, K_theta); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
   y_pos ~ gamma(shape_param, rate_param);  // shape is mean over scale, rate is 1 over scale

  // Estimate zeros using logit model
	y_zero ~ bernoulli_logit(b_0_theta + X * b_theta);

}
generated quantities {
  real rho_shape ;
  real rho_rate ;
  real rho_theta;

  rho_shape = exp(lambda_shape);
  rho_rate = exp(lambda_rate);
  rho_theta = exp(lambda_theta);
}
"

## //  y_pos ~ gamma((X_pos * b_mean + b_0_mean) ./ (X_pos * b_scale + b_0_scale), ones_pos ./ (X_pos * b_scale + b_0_scale));  // shape is mean over scale, rate is 1 over scale
## // y_pos ~ gamma(shape_param, rate_param);  // shape is mean over scale, rate is 1 over scale

#rho_mean = log(lambda_mean);
#rho_scale = log(lambda_scale);
#rho_theta = log(lambda_theta);

	
###lambda_shape_log ~ uniform(-5,0);#
###  lambda_rate_log ~ uniform(-5,0);
###  lambda_theta_log ~ uniform(-5, 0);
##secderiv_mean = b_mean' * S * b_mean; 
##secderiv_scale = b_scale' * S * b_scale; 
##




#' Cyclic Stan model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @export

cyclic_model_diag <- "
data {
  int<lower=0> N;  //number of points
  int<lower=0> N_pos;  //number of days with positive precip

  int<lower=0> basis_dim;

  vector[N_pos] y_pos;     // list of y values with positive precip
  int<lower=0,upper=1> y_zero[N];     // list of one or zero for positive precip
 
  matrix[N_pos,basis_dim] X_pos;  // Basis for positive precip days
  matrix[N,basis_dim] X;   // Basis for all days

  vector[2] b_0_shape_prior; 
  vector[2] b_0_rate_prior; 
  vector[2] b_0_theta_prior; 

  real lambda_shape_prior; 
  real lambda_rate_prior; 
  real lambda_theta_prior; 
 }
transformed data {  
  vector[basis_dim] zeros; 
  vector[N_pos] ones_pos; 

  zeros = rep_vector(0, basis_dim) ;
  ones_pos = rep_vector(1, N_pos) ;
}
parameters {
  real b_0_shape;
  real b_0_rate;  
  real b_0_theta;  
 
  vector[basis_dim] b_shape;  
  vector[basis_dim] b_rate;   
  vector[basis_dim] b_theta; 
 
  real lambda_shape ;
  real lambda_rate ;
  real lambda_theta ;

//  real rho_shape ;
//  real rho_rate ;
//  real rho_theta;
}
transformed parameters { 
  vector[N_pos] shape_param;  
  vector[N_pos] rate_param;

//  real lambda_shape ;
//  real lambda_rate ;
//  real lambda_theta ;

  shape_param = exp(X_pos * b_shape + b_0_shape);
  rate_param = exp(X_pos * b_rate + b_0_rate);

// lambda_shape = exp(rho_shape);
// lambda_rate = exp(rho_rate);
// lambda_theta = exp(rho_theta);
} 
model {
 // lambda_mean ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
 // lambda_scale ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
 // lambda_theta ~ gamma(lambda_theta_prior[1], lambda_theta_prior[2]);

  lambda_shape ~ cauchy(0, lambda_shape_prior);
  lambda_rate ~ cauchy(0, lambda_rate_prior);
  lambda_theta ~ cauchy(0, lambda_theta_prior);

//  rho_shape ~ uniform(-12, 12);
//  rho_rate ~ uniform(-12, 12);
//  rho_theta ~ uniform(-12, 12);


   b_0_shape  ~ normal(b_0_shape_prior[1],b_0_shape_prior[2]);   
   b_shape ~  normal(zeros,lambda_shape); 
 
   b_0_rate  ~ normal(b_0_rate_prior[1],b_0_rate_prior[2]);  
   b_rate ~   normal(zeros,lambda_rate); 

   b_0_theta ~ normal(b_0_theta_prior[1],b_0_theta_prior[2]);   
   b_theta ~  normal(zeros,lambda_theta); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
   y_pos ~ gamma(shape_param, rate_param);  // shape is mean over scale, rate is 1 over scale

  // Estimate zeros using logit model
	y_zero ~ bernoulli_logit(b_0_theta + X * b_theta);

}
generated quantities {
  real rho_shape ;
  real rho_rate ;
  real rho_theta;

  rho_shape = exp(lambda_shape);
  rho_rate = exp(lambda_rate);
  rho_theta = exp(lambda_theta);

}
"






#' Tensor Product Stan model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @export
tensor_model <- '
data {
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
	
  int<lower=0> basis_dim;
  
  matrix[N,basis_dim] X; 
  // Could declare cov_matrix but this involved a repeated check and so slows things down
  matrix[basis_dim, basis_dim] S_1; 
  matrix[basis_dim, basis_dim] S_2; 
  
  vector[2] b_0_mean_prior; 
  vector[2] b_0_scale_prior; 

  vector[4] lambda_mean_prior;
  vector[4] lambda_scale_prior;
 }
transformed data {  
 vector[basis_dim] zero; 
 vector[N] ones;

 zero = rep_vector(0, basis_dim) ;

 ones = rep_vector(1, N); 
}
parameters {
  real b_0_mean;
  real b_0_scale;  
  
  vector[basis_dim] b_mean;  
  vector[basis_dim] b_scale;   

  real<lower=0> lambda_mean_first ;
  real<lower=0> lambda_mean_second ;
  real<lower=0> lambda_scale_first ;
  real<lower=0> lambda_scale_second ;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_mean; 
  matrix[basis_dim, basis_dim] K_scale; 
  
  vector<lower=0>[N] mean_param;  
  vector<lower=0>[N] scale_param;
  
  K_mean = S_1 * lambda_mean_first  + S_2 * lambda_mean_second ;
  K_scale = S_1 * lambda_scale_first  + S_2 * lambda_scale_second ;
   
  mean_param = to_vector(X * b_mean) + b_0_mean;
  scale_param = to_vector(X * b_scale) + b_0_scale;
} 
model {
    lambda_mean_first ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
    lambda_mean_second ~ gamma(lambda_mean_prior[3], lambda_mean_prior[4]);

    lambda_scale_first ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
    lambda_scale_second ~ gamma(lambda_scale_prior[3], lambda_scale_prior[4]);
	
    b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
    b_mean ~ multi_normal_prec(zero,K_mean);   

    b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
    b_scale ~ multi_normal_prec(zero,K_scale); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
  // shape is mean over scale, rate is 1 over scale
  y ~ gamma(mean_param ./ scale_param, ones ./ scale_param);  
 
}
generated quantities {
}
'



#' Tensor Product Stan model with fixed lambda in the year axis
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @export
tensor_model_lambdafixed <- '
data {
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
	
  int<lower=0> basis_dim;
  
  matrix[N,basis_dim] X; 
  cov_matrix[basis_dim] S_1; 
  cov_matrix[basis_dim] S_2; 
  
  vector[2] b_0_mean_prior; 
  vector[2] b_0_scale_prior; 

  vector[3] lambda_mean_prior;
  vector[3] lambda_scale_prior;

 }
transformed data {  
 vector[basis_dim] zero; 
 vector[N] ones;

 zero = rep_vector(0, basis_dim) ;

 ones = rep_vector(1, N); 
}
parameters {
  real b_0_mean;
  real b_0_scale;  
  
  vector[basis_dim] b_mean;  
  vector[basis_dim] b_scale;   

  real<lower=0> lambda_mean_first ;
  real<lower=0> lambda_scale_first ;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_mean; 
  matrix[basis_dim, basis_dim] K_scale; 
  
  vector<lower=0>[N] mean_param;  
  vector<lower=0>[N] scale_param;
  
  K_mean = S_1 * lambda_mean_first  + S_2 * lambda_mean_prior[3] ;
  K_scale = S_1 * lambda_scale_first  + S_2 * lambda_scale_prior[3] ;
   
  mean_param = to_vector(X * b_mean) + b_0_mean;
  scale_param = to_vector(X * b_scale) + b_0_scale;
} 
model {
  lambda_mean_first ~  uniform(lambda_mean_prior[1], lambda_mean_prior[2]);
  lambda_scale_first ~ uniform(lambda_scale_prior[1], lambda_scale_prior[2]);
	
   //lambda_mean_first ~  gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
   //lambda_scale_first ~  gamma(lambda_mean_prior[1], lambda_mean_prior[2]);

   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(zero,K_mean);   

   b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
   b_scale ~ multi_normal_prec(zero,K_scale); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
  y ~ gamma(mean_param ./ scale_param, ones ./ scale_param);  
 
}
generated quantities {
}
'





