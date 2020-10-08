data {
  int<lower=0> N;  //number of points
  int<lower=0> N_pos;  //number of days with positive precip

  vector[N_pos] y_pos;     // list of y values with positive precip
  int<lower=0,upper=1> y_zero[N];     // list of one or zero for positive precip

  // basis dimension for the 1 ti piece
  int<lower=0> basis_dim_mean;
  int<lower=0> basis_dim_disp;
  int<lower=0> basis_dim_theta;

  // Read in the model matrices
  matrix[N_pos,basis_dim_mean] X_mean_jdate ;
  matrix[N_pos,basis_dim_disp] X_disp_jdate ;
  matrix[N,basis_dim_theta] X_theta_jdate ;

  // Read in penalty matrices
  matrix[basis_dim_mean, basis_dim_mean] s_mean_jdate; 
  matrix[basis_dim_disp, basis_dim_disp] s_disp_jdate; 
  matrix[basis_dim_theta, basis_dim_theta] s_theta_jdate; 

  matrix[1,2] lambda_mean_prior;
  matrix[1,2] lambda_disp_prior;
  matrix[1,2] lambda_theta_prior;

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
 
  vector[basis_dim_mean] b_mean_jdate;  
  vector[basis_dim_disp] b_disp_jdate;  
  vector[basis_dim_theta] b_theta_jdate;  

  vector<lower = 0>[1]  lambda_mean ;
  vector<lower = 0>[1]  lambda_disp ;
  vector<lower = 0>[1]  lambda_theta ;
}
transformed parameters { 
  vector[N_pos] mean_param;  
  vector[N_pos] disp_param;

  matrix[basis_dim_mean, basis_dim_mean] K_mean_jdate; 
  matrix[basis_dim_disp, basis_dim_disp] K_disp_jdate; 
  matrix[basis_dim_theta, basis_dim_theta] K_theta_jdate; 

  // Penalize for mean
  // S3 and S5 are the double penalty for null space
  K_mean_jdate = s_mean_jdate * lambda_mean[1] ;

  // Penalize for disp
  K_disp_jdate = s_disp_jdate * lambda_disp[1] ;

  // Penalize for theta
  K_theta_jdate = s_theta_jdate * lambda_theta[1] ;

  // Calculate mean and dispersion
  mean_param = exp(b_0_mean + X_mean_jdate * b_mean_jdate) ;

  disp_param = exp(seven_vec + log(one_vec + exp(b_0_disp + X_disp_jdate * b_disp_jdate)));

} 
model {
  lambda_mean[1] ~ gamma(lambda_mean_prior[1,1], lambda_mean_prior[1,2]);
  lambda_disp[1] ~ gamma(lambda_disp_prior[1,1], lambda_disp_prior[1,2]);
  lambda_theta[1] ~ gamma(lambda_theta_prior[1,1], lambda_theta_prior[1,2]);

   b_0_mean ~ normal(b_0_prior[1,1], b_0_prior[1,2]);   
   b_mean_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_mean), K_mean_jdate); 

   b_0_disp  ~ normal(b_0_prior[2,1], b_0_prior[2,2]);   
   b_disp_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_disp), K_disp_jdate); 


   b_0_theta ~ normal(b_0_prior[3,1], b_0_prior[3,2]);   
   b_theta_jdate ~ multi_normal_prec(rep_vector(0, basis_dim_theta), K_theta_jdate); 
  
  // Estimate y values using a gamma distribution, Stan uses shape and inverse scale parameters
   y_pos ~ gamma(one_vec ./ disp_param, one_vec ./ (mean_param .* disp_param) );   // 

  // Estimate zeros using logit model
	y_zero ~ bernoulli_logit(b_0_theta + X_theta_jdate * b_theta_jdate);

}
generated quantities {

}

