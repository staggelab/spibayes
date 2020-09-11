
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

  vector[2]  rho_shape_prior; 
  vector[2]  rho_rate_prior; 
  vector[2]  rho_theta_prior; 
 }
transformed data {  
  vector[basis_dim] zeros; 
//  vector[N_pos] ones_pos; 

  zeros = rep_vector(0, basis_dim) ;
//  ones_pos = rep_vector(1, N_pos) ;
}
parameters {
  real b_0_shape;
  real b_0_rate;  
  real b_0_theta;  
 
  vector[basis_dim] b_shape;  
  vector[basis_dim] b_rate;   
  vector[basis_dim] b_theta; 
 
  real rho_shape ;
  real rho_rate ;
  real rho_theta;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_shape; 
  matrix[basis_dim, basis_dim] K_rate; 
  matrix[basis_dim, basis_dim] K_theta; 
 
  vector[N_pos] shape_param;  
  vector[N_pos] rate_param;

  real lambda_shape ;
  real lambda_rate ;
  real lambda_theta ;

  lambda_shape = exp(rho_shape);
  lambda_rate = exp(rho_rate);
  lambda_theta = exp(rho_theta);

  K_shape = S * lambda_shape ;
  K_rate = S * lambda_rate ;
  K_theta = S * lambda_theta ;
  
  shape_param = exp(X_pos * b_shape + b_0_shape);
  rate_param = exp(X_pos * b_rate + b_0_rate);

} 
model {
  rho_shape ~ uniform(rho_shape_prior[1], rho_shape_prior[2]);
  rho_rate ~ uniform(rho_rate_prior[1], rho_rate_prior[2]);
  rho_theta ~ uniform(rho_theta_prior[1], rho_theta_prior[2]);

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
tensor_model <- "
data {
  int<lower=0> N;  //number of points
  int<lower=0> N_pos;  //number of days with positive precip

  vector[N_pos] y_pos;     // list of y values with positive precip
  int<lower=0,upper=1> y_zero[N];     // list of one or zero for positive precip

  // basis dimension for the 3 ti pieces
  int<lower=0> basis_dim_mean[3];
  int<lower=0> basis_dim_disp[3];
  int<lower=0> basis_dim_theta[3];

  // Read in the model matrices
  matrix[N_pos,basis_dim_mean[1]] X_mean_jdate ;
  matrix[N_pos,basis_dim_mean[2]] X_mean_year ;
  matrix[N_pos,basis_dim_mean[3]] X_mean_tensor ;

  matrix[N_pos,basis_dim_disp[1]] X_disp_jdate ;
  matrix[N_pos,basis_dim_disp[2]] X_disp_year ;
  matrix[N_pos,basis_dim_disp[3]] X_disp_tensor ;

  matrix[N,basis_dim_theta[1]] X_theta_jdate ;
  matrix[N,basis_dim_theta[2]] X_theta_year ;
  matrix[N,basis_dim_theta[3]] X_theta_tensor ;

  // Read in penalty matrices
  matrix[basis_dim_mean[1], basis_dim_mean[1]] s_mean_jdate; 
  matrix[basis_dim_mean[2], basis_dim_mean[2]] s_mean_year; 
  matrix[basis_dim_mean[2], basis_dim_mean[2]] s_mean_year_double;
  matrix[basis_dim_mean[3], basis_dim_mean[3]] s_mean_tensor_jdate; 
  matrix[basis_dim_mean[3], basis_dim_mean[3]] s_mean_tensor_year; 

  matrix[basis_dim_disp[1], basis_dim_disp[1]] s_disp_jdate; 
  matrix[basis_dim_disp[2], basis_dim_disp[2]] s_disp_year; 
  matrix[basis_dim_disp[2], basis_dim_disp[2]] s_disp_year_double;
  matrix[basis_dim_disp[3], basis_dim_disp[3]] s_disp_tensor_jdate; 
  matrix[basis_dim_disp[3], basis_dim_disp[3]] s_disp_tensor_year; 

  matrix[basis_dim_theta[1], basis_dim_theta[1]] s_theta_jdate; 
  matrix[basis_dim_theta[2], basis_dim_theta[2]] s_theta_year; 
  matrix[basis_dim_theta[2], basis_dim_theta[2]] s_theta_year_double;
  matrix[basis_dim_theta[3], basis_dim_theta[3]] s_theta_tensor_jdate; 
  matrix[basis_dim_theta[3], basis_dim_theta[3]] s_theta_tensor_year; 

  matrix[5,2] lambda_mean_prior;
  matrix[5,2] lambda_disp_prior;
  matrix[5,2] lambda_theta_prior;

  matrix[3,2] b_0_prior;
 }
transformed data {  
  vector[N_pos] one_vec;
  vector[N_pos] seven_vec;

  one_vec = rep_vector(1, N_pos) ;
  seven_vec = rep_vector(-7, N_pos) ;
}
parameters {
  // Come back you might want to put a zero bound
  real b_0_mean; 
  real b_0_disp; 
  real b_0_theta;  
 
  vector[basis_dim_mean[1]] b_mean_jdate;  
  vector[basis_dim_mean[2]] b_mean_year;  
  vector[basis_dim_mean[3]] b_mean_tensor;  

  vector[basis_dim_disp[1]] b_disp_jdate;  
  vector[basis_dim_disp[2]] b_disp_year;  
  vector[basis_dim_disp[3]] b_disp_tensor; 

  vector[basis_dim_theta[1]] b_theta_jdate;  
  vector[basis_dim_theta[2]] b_theta_year;  
  vector[basis_dim_theta[3]] b_theta_tensor; 

  vector<lower = 0>[5]  lambda_mean ;
  vector<lower = 0>[5]  lambda_disp ;
  vector<lower = 0>[5]  lambda_theta ;
}
transformed parameters { 
  vector[N_pos] mean_param;  
  // vector[N_pos] second_param; //<upper = 500>[
  vector[N_pos] disp_param;

  matrix[basis_dim_mean[1], basis_dim_mean[1]] K_mean_jdate; 
  matrix[basis_dim_mean[2], basis_dim_mean[2]] K_mean_year; 
  matrix[basis_dim_mean[3], basis_dim_mean[3]] K_mean_tensor; 

  matrix[basis_dim_disp[1], basis_dim_disp[1]] K_disp_jdate; 
  matrix[basis_dim_disp[2], basis_dim_disp[2]] K_disp_year; 
  matrix[basis_dim_disp[3], basis_dim_disp[3]] K_disp_tensor; 

  matrix[basis_dim_theta[1], basis_dim_theta[1]] K_theta_jdate; 
  matrix[basis_dim_theta[2], basis_dim_theta[2]] K_theta_year; 
  matrix[basis_dim_theta[3], basis_dim_theta[3]] K_theta_tensor; 

  // Penalize for mean
  // S3 and S5 are the double penalty for null space
  K_mean_jdate = s_mean_jdate * lambda_mean[1] ;
  K_mean_year = s_mean_year * lambda_mean[2] + s_mean_year_double * lambda_mean[3] ; 
  K_mean_tensor = s_mean_tensor_jdate * lambda_mean[4] + s_mean_tensor_year * lambda_mean[5]    ;

  // Penalize for disp
  K_disp_jdate = s_disp_jdate * lambda_disp[1] ;
  K_disp_year = s_disp_year * lambda_disp[2] + s_disp_year_double * lambda_disp[3] ; 
  K_disp_tensor = s_disp_tensor_jdate * lambda_disp[4] + s_disp_tensor_year * lambda_disp[5]    ;

  // Penalize for theta
  K_theta_jdate = s_theta_jdate * lambda_theta[1] ;
  K_theta_year = s_theta_year * lambda_theta[2] + s_theta_year_double * lambda_theta[3] ; 
  K_theta_tensor = s_theta_tensor_jdate * lambda_theta[4] + s_theta_tensor_year * lambda_theta[5]    ;


  // Calculate mean and dispersion
  mean_param = exp(b_0_mean + X_mean_jdate * b_mean_jdate + X_mean_year * b_mean_year + X_mean_tensor * b_mean_tensor) ;

  disp_param = exp(seven_vec + log(one_vec + exp(b_0_disp + X_disp_jdate * b_disp_jdate + X_disp_year * b_disp_year + X_disp_tensor * b_disp_tensor)));

} 
model {
  lambda_mean[1] ~ gamma(lambda_mean_prior[1,1], lambda_mean_prior[1,2]);
  lambda_mean[2] ~ gamma(lambda_mean_prior[2,1], lambda_mean_prior[2,2]);
  lambda_mean[3] ~ gamma(lambda_mean_prior[3,1], lambda_mean_prior[3,2]);
  lambda_mean[4] ~ gamma(lambda_mean_prior[4,1], lambda_mean_prior[4,2]);
  lambda_mean[5] ~ gamma(lambda_mean_prior[5,1], lambda_mean_prior[5,2]);

  lambda_disp[1] ~ gamma(lambda_disp_prior[1,1], lambda_disp_prior[1,2]);
  lambda_disp[2] ~ gamma(lambda_disp_prior[2,1], lambda_disp_prior[2,2]);
  lambda_disp[3] ~ gamma(lambda_disp_prior[3,1], lambda_disp_prior[3,2]);
  lambda_disp[4] ~ gamma(lambda_disp_prior[4,1], lambda_disp_prior[4,2]);
  lambda_disp[5] ~ gamma(lambda_disp_prior[5,1], lambda_disp_prior[5,2]);

  lambda_theta[1] ~ gamma(lambda_theta_prior[1,1], lambda_theta_prior[1,2]);
  lambda_theta[2] ~ gamma(lambda_theta_prior[2,1], lambda_theta_prior[2,2]);
  lambda_theta[3] ~ gamma(lambda_theta_prior[3,1], lambda_theta_prior[3,2]);
  lambda_theta[4] ~ gamma(lambda_theta_prior[4,1], lambda_theta_prior[4,2]);
  lambda_theta[5] ~ gamma(lambda_theta_prior[5,1], lambda_theta_prior[5,2]);

   b_0_mean ~ normal(b_0_prior[1,1], b_0_prior[1,2]);   
   b_mean_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_mean[1]), K_mean_jdate); 
   b_mean_year ~ multi_normal_prec(rep_vector(0, basis_dim_mean[2]), K_mean_year); 
   b_mean_tensor ~ multi_normal_prec(rep_vector(0, basis_dim_mean[3]), K_mean_tensor); 

   b_0_disp  ~ normal(b_0_prior[2,1], b_0_prior[2,2]);   
   b_disp_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_disp[1]), K_disp_jdate); 
   b_disp_year ~ multi_normal_prec(rep_vector(0, basis_dim_disp[2]), K_disp_year); 
   b_disp_tensor ~ multi_normal_prec(rep_vector(0, basis_dim_disp[3]), K_disp_tensor); 

   b_0_theta ~ normal(b_0_prior[3,1], b_0_prior[3,2]);   
   b_theta_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_theta[1]), K_theta_jdate); 
   b_theta_year ~ multi_normal_prec(rep_vector(0, basis_dim_theta[2]), K_theta_year); 
   b_theta_tensor ~ multi_normal_prec(rep_vector(0, basis_dim_theta[3]), K_theta_tensor); 
  
  // Estimate y values using a gamma distribution, Stan uses shape and inverse scale parameters
   y_pos ~ gamma(one_vec ./ disp_param, one_vec ./ (mean_param .* disp_param) );   // 

  // Estimate zeros using logit model
	y_zero ~ bernoulli_logit(b_0_theta + X_theta_jdate * b_theta_jdate + X_theta_year * b_theta_year + X_theta_tensor * b_theta_tensor);

}
generated quantities {

}
"






#' Tensor Product Nozero
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
tensor_model_nozero <- "
data {
  int<lower=0> N;  //number of points
  int<lower=0> N_pos;  //number of days with positive precip

  vector[N_pos] y_pos;     // list of y values with positive precip
  int<lower=0,upper=1> y_zero[N];     // list of one or zero for positive precip

  // basis dimension for the 3 ti pieces
  int<lower=0> basis_dim_mean[3];
  int<lower=0> basis_dim_disp[3];

  // Read in the model matrices
  matrix[N_pos,basis_dim_mean[1]] X_mean_jdate ;
  matrix[N_pos,basis_dim_mean[2]] X_mean_year ;
  matrix[N_pos,basis_dim_mean[3]] X_mean_tensor ;

  matrix[N_pos,basis_dim_disp[1]] X_disp_jdate ;
  matrix[N_pos,basis_dim_disp[2]] X_disp_year ;
  matrix[N_pos,basis_dim_disp[3]] X_disp_tensor ;

  // Read in penalty matrices
  matrix[basis_dim_mean[1], basis_dim_mean[1]] s_mean_jdate; 
  matrix[basis_dim_mean[2], basis_dim_mean[2]] s_mean_year; 
  matrix[basis_dim_mean[2], basis_dim_mean[2]] s_mean_year_double;
  matrix[basis_dim_mean[3], basis_dim_mean[3]] s_mean_tensor_jdate; 
  matrix[basis_dim_mean[3], basis_dim_mean[3]] s_mean_tensor_year; 

  matrix[basis_dim_disp[1], basis_dim_disp[1]] s_disp_jdate; 
  matrix[basis_dim_disp[2], basis_dim_disp[2]] s_disp_year; 
  matrix[basis_dim_disp[2], basis_dim_disp[2]] s_disp_year_double;
  matrix[basis_dim_disp[3], basis_dim_disp[3]] s_disp_tensor_jdate; 
  matrix[basis_dim_disp[3], basis_dim_disp[3]] s_disp_tensor_year; 

  matrix[5,2] lambda_mean_prior;
  matrix[5,2] lambda_disp_prior;

  matrix[3,2] b_0_prior;
 }
transformed data {  
  vector[N_pos] one_vec;
  vector[N_pos] seven_vec;

  one_vec = rep_vector(1, N_pos) ;
  seven_vec = rep_vector(-7, N_pos) ;
}
parameters {
  // Come back you might want to put a zero bound
  real b_0_mean; 
  real b_0_disp; 
  real b_0_theta;  
 
  vector[basis_dim_mean[1]] b_mean_jdate;  
  vector[basis_dim_mean[2]] b_mean_year;  
  vector[basis_dim_mean[3]] b_mean_tensor;  

  vector[basis_dim_disp[1]] b_disp_jdate;  
  vector[basis_dim_disp[2]] b_disp_year;  
  vector[basis_dim_disp[3]] b_disp_tensor; 

  vector<lower = 0>[5]  lambda_mean ;
  vector<lower = 0>[5]  lambda_disp ;
}
transformed parameters { 
  vector[N_pos] mean_param;  
  // vector[N_pos] second_param; //<upper = 500>[
  vector[N_pos] disp_param;

  matrix[basis_dim_mean[1], basis_dim_mean[1]] K_mean_jdate; 
  matrix[basis_dim_mean[2], basis_dim_mean[2]] K_mean_year; 
  matrix[basis_dim_mean[3], basis_dim_mean[3]] K_mean_tensor; 

  matrix[basis_dim_disp[1], basis_dim_disp[1]] K_disp_jdate; 
  matrix[basis_dim_disp[2], basis_dim_disp[2]] K_disp_year; 
  matrix[basis_dim_disp[3], basis_dim_disp[3]] K_disp_tensor; 

  // Penalize for mean
  // S3 and S5 are the double penalty for null space
  K_mean_jdate = s_mean_jdate * lambda_mean[1] ;
  K_mean_year = s_mean_year * lambda_mean[2] + s_mean_year_double * lambda_mean[3] ; 
  K_mean_tensor = s_mean_tensor_jdate * lambda_mean[4] + s_mean_tensor_year * lambda_mean[5]    ;

  // Penalize for disp
  K_disp_jdate = s_disp_jdate * lambda_disp[1] ;
  K_disp_year = s_disp_year * lambda_disp[2] + s_disp_year_double * lambda_disp[3] ; 
  K_disp_tensor = s_disp_tensor_jdate * lambda_disp[4] + s_disp_tensor_year * lambda_disp[5]    ;

  // Calculate mean and dispersion
  mean_param = exp(b_0_mean + X_mean_jdate * b_mean_jdate + X_mean_year * b_mean_year + X_mean_tensor * b_mean_tensor) ;

  disp_param = exp(seven_vec + log(one_vec + exp(b_0_disp + X_disp_jdate * b_disp_jdate + X_disp_year * b_disp_year + X_disp_tensor * b_disp_tensor)));

} 
model {
  lambda_mean[1] ~ gamma(lambda_mean_prior[1,1], lambda_mean_prior[1,2]);
  lambda_mean[2] ~ gamma(lambda_mean_prior[2,1], lambda_mean_prior[2,2]);
  lambda_mean[3] ~ gamma(lambda_mean_prior[3,1], lambda_mean_prior[3,2]);
  lambda_mean[4] ~ gamma(lambda_mean_prior[4,1], lambda_mean_prior[4,2]);
  lambda_mean[5] ~ gamma(lambda_mean_prior[5,1], lambda_mean_prior[5,2]);

  lambda_disp[1] ~ gamma(lambda_disp_prior[1,1], lambda_disp_prior[1,2]);
  lambda_disp[2] ~ gamma(lambda_disp_prior[2,1], lambda_disp_prior[2,2]);
  lambda_disp[3] ~ gamma(lambda_disp_prior[3,1], lambda_disp_prior[3,2]);
  lambda_disp[4] ~ gamma(lambda_disp_prior[4,1], lambda_disp_prior[4,2]);
  lambda_disp[5] ~ gamma(lambda_disp_prior[5,1], lambda_disp_prior[5,2]);

   b_0_mean ~ normal(b_0_prior[1,1], b_0_prior[1,2]);   
   b_mean_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_mean[1]), K_mean_jdate); 
   b_mean_year ~ multi_normal_prec(rep_vector(0, basis_dim_mean[2]), K_mean_year); 
   b_mean_tensor ~ multi_normal_prec(rep_vector(0, basis_dim_mean[3]), K_mean_tensor); 

   b_0_disp  ~ normal(b_0_prior[2,1], b_0_prior[2,2]);   
   b_disp_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_disp[1]), K_disp_jdate); 
   b_disp_year ~ multi_normal_prec(rep_vector(0, basis_dim_disp[2]), K_disp_year); 
   b_disp_tensor ~ multi_normal_prec(rep_vector(0, basis_dim_disp[3]), K_disp_tensor); 

   b_0_theta ~ normal(b_0_prior[3,1], b_0_prior[3,2]);   
  
  // Estimate y values using a gamma distribution, Stan uses shape and inverse scale parameters
   y_pos ~ gamma(one_vec ./ disp_param, one_vec ./ (mean_param .* disp_param) );   // 

  // Estimate zeros using logit model
	y_zero ~ bernoulli_logit(b_0_theta);

}
generated quantities {

}
"

