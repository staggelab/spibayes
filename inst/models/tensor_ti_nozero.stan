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
