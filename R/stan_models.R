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
  
  vector[basis_dim] b_mean_prior; 
  vector[basis_dim] b_scale_prior; 

  vector[2] lambda_mean_prior; 
  vector[2] lambda_scale_prior; 
 }
transformed data {  
  vector[basis_dim] b_theta_prior; 
  b_theta_prior = rep_vector(0.5, basis_dim);
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
}
transformed parameters { 
  matrix[basis_dim, basis_dim] K_mean; 
  matrix[basis_dim, basis_dim] K_scale; 
  matrix[basis_dim, basis_dim] K_theta; 
 
  vector<lower=0>[N] mean_param;  
  vector<lower=0>[N] scale_param;
  vector[N] theta_param;
 
  K_mean = S * lambda_mean ;
  K_scale = S * lambda_scale ;
  K_theta = S * lambda_theta ;
  
  mean_param = to_vector(X * b_mean) + b_0_mean;
  scale_param = to_vector(X * b_scale) + b_0_scale;
  theta_param = to_vector(X * b_theta) + b_0_theta;
} 
model {
  lambda_mean ~ gamma(lambda_mean_prior[1], lambda_mean_prior[2]);
  lambda_scale ~ gamma(lambda_scale_prior[1], lambda_scale_prior[2]);
  lambda_theta ~ gamma(5, 5/50);
	
   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(b_mean_prior,K_mean); 
 
   b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
   b_scale ~ multi_normal_prec(b_scale_prior,K_scale); 

   b_0_theta ~ normal(-1,2);   
   b_theta ~ multi_normal_prec(b_theta_prior, K_theta); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
 // y ~ gamma(mean_param ./ scale_param, rep_vector(1, N) ./ scale_param);  // shape is mean over scale, rate is 1 over scale

  for(n in 1:N){
   if (y[n] == 0)
      1 ~ bernoulli_logit(theta_param[n]);
    else {
      0 ~ bernoulli_logit(theta_param[n]);
      y[n] ~ gamma(mean_param[n] / scale_param[n], 1/scale_param[n] );
    }	
  }  
 
}
generated quantities {
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
  cov_matrix[basis_dim] S_1; 
  cov_matrix[basis_dim] S_2; 
  
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





