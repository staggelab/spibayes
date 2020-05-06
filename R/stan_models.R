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
cyclic_model <- '
data {
  int<lower=0> N;  //number of points
  vector[N] y;     // list of y values from data
	
  int<lower=0> basis_dim;
  
  matrix[N,basis_dim] X; 
  matrix[basis_dim, basis_dim] S; 
    
  vector[2] b_0_mean_prior; 
  vector[2] b_0_scale_prior; 
  
 }
transformed data {  
 vector[basis_dim] zero; 
 zero = rep_vector(0, basis_dim) ;
 
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
  
  vector[N] mean_param;  
  vector[N] scale_param;
  
  K_mean = S * lambda_mean ;
  K_scale = S * lambda_scale ;
   
  mean_param = to_vector(X * b_mean) + b_0_mean;
  scale_param = to_vector(X * b_scale) + b_0_scale;
} 
model {
 
  //lambda_mean ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
  //lambda_scale ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
  lambda_mean ~ gamma(0.05,0.005);
  lambda_scale ~ gamma(0.05,0.005);
	
   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(zero,K_mean); 
 
   b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
   b_scale ~ multi_normal_prec(zero,K_scale); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
  y ~ gamma(exp(mean_param) ./ exp(scale_param), rep_vector(1, N) ./ exp(scale_param));  // shape is mean over scale, rate is 1 over scale
 
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
'





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
  matrix[basis_dim,basis_dim] S_1; 
  matrix[basis_dim,basis_dim] S_2; 
  matrix[basis_dim,basis_dim] S_3; 
  
  vector[2] b_0_mean_prior; 
  vector[2] b_0_disp_prior; 
 }
transformed data {  
 vector[basis_dim] zero; 
 zero = rep_vector(0, basis_dim) ;
 
}
parameters {
  real b_0_mean;
  real b_0_disp;  
  
  vector[basis_dim] b_mean;  
  vector[basis_dim] b_disp;   
  
  vector<lower=0>[3] lambda_mean ;
  vector<lower=0>[3] lambda_disp ;
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_mean; 
  matrix[basis_dim, basis_dim] K_disp; 
  
  vector<lower=0> [N] mean_param;  
  vector<lower=0> [N] disp_param;
  
  K_mean = S_1 * lambda_mean[1]  + S_2 * lambda_mean[2] + S_3 * lambda_mean[3] ;
  K_disp = S_1 * lambda_disp[1]  + S_2 * lambda_disp[2] + S_3 * lambda_disp[3] ;
   
  mean_param = to_vector(X * b_mean) + b_0_mean;
  disp_param = to_vector(X * b_disp) + b_0_disp;
} 
model {
 
  lambda_mean ~ gamma(0.05,0.005);
  lambda_disp ~ gamma(0.05,0.005);
	
   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(zero,K_mean); 

   b_0_disp  ~ normal(b_0_disp_prior[1],b_0_disp_prior[2]);  
   b_disp ~ multi_normal_prec(zero,K_disp); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
  y ~ gamma(rep_vector(1, N) ./ exp(disp_param), rep_vector(1, N) ./ ( mean_param /* disp_param)) ;  // shape is 1 over disp, scale is mean over shape or 1 over mean times disp
 
}
generated quantities { }
'



