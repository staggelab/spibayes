#' Run cyclic model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @import rstan
#' @return A matrix of the infile
#' @export
spi_fit<- function(spi_input, n_chains=1, iter=1000, cores = 1){
	require(rstan)

	### Extract important values from the input object
	y <-  spi_input$data$precip
	X <- spi_input$x_reparam

	### Calculate some necesary dimensions
	N <- length(y)
	basis_dim <- dim(X)[2]
	
	### Error check that the basis and y values are of the same length
	assertthat::assert_that(N == dim(X)[1])

	### Create the data to send to stan model
	data_fitting <- list(N = N, 
		basis_dim = basis_dim, 
		y= y, 
		X = X,  
		b_0_mean_prior=spi_input$b_0$mean, 
		b_0_scale_prior=spi_input$b_0$scale)
	
	### Use the estimated lambda or 100 - whichever is smaller
	lambda_mean_init <- sapply(spi_input$lambda_init$mean, function(x){min(x, 100)})
	lambda_scale_init <- sapply(spi_input$lambda_init$scale, function(x){min(x, 100)})

	### Create the initial values
	init_vals <- list(list(
		b_0_mean = spi_input$b_0$mean[1], 
		b_0_scale = spi_input$b_0$scale[1], 
		b_mean = spi_input$b_init$mean, 
		b_scale = spi_input$b_init$scale, 
		lambda_mean = lambda_mean_init, 
		lambda_scale = lambda_scale_init)
	)

	### Loop through the initial values if we have more than one chain
	### Vary lambda values 2 orders of magnitude smaller or larger
	if(n_chains > 1) {
		for (j in seq(2, n_chains)) {
			init_vals[[j]] <- init_vals[[1]]
			init_vals[[j]]$lambda_mean <- lambda_mean_init*(10^runif(1, -2,2))
			init_vals[[j]]$lambda_scale <- lambda_scale_init*(10^runif(1, -2,2))
		}	
	}

	### Run the cyclic model if cyclic, tensor if tensor
	if (spi_input$type == "cyclic"){

		### Insert the S penaly matrix, which is different shape for cyclic than tensor
		data_fitting[["S"]] <- spi_input$s_reparam

		### Run model
		#model_fit <- spi_cyclic(data = data_fitting, init_vals = init_vals, n_chains = n_chains, iter = iter, cores = cores)

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
  lambda_mean ~ gamma(0.05,0.005);
  lambda_scale ~ gamma(0.05,0.005);
	
   b_0_mean  ~ normal(b_0_mean_prior[1],b_0_mean_prior[2]);   
   b_mean ~ multi_normal_prec(zero,K_mean); 
 
   b_0_scale  ~ normal(b_0_scale_prior[1],b_0_scale_prior[2]);  
   b_scale ~ multi_normal_prec(zero,K_scale); 
  
  // Estimate y values using a gamma distribution, Stan uses rate, rather than scale parameter
  y ~ gamma(exp(mean_param) ./ exp(scale_param), rep_vector(1, N) ./ exp(scale_param));  // shape is mean over scale, rate is 1 over scale
 
}
'
		### Fit the model
		model_fit <- rstan::stan(model_code = cyclic_model, 
			data = data_fitting, 
			init = init_vals,
			iter = iter, 
			chains = n_chains,
			cores = cores, 
			verbose = FALSE)

	} else if (spi_input$type == "tensor"){
		
		### Insert the 3 S penaly matrices
		data_fitting[["S_1"]] <- spi_input$s_reparam[[1]]
		data_fitting[["S_2"]] <- spi_input$s_reparam[[2]]
		data_fitting[["S_3"]] <- spi_input$s_reparam[[3]]

		### Run model
		model_fit <- spi_tensor(data = data_fitting, init_vals = init_vals, n_chains = n_chains, iter = iter, cores = cores)

	}

	return(model_fit)
}


#' Run cyclic model
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
spi_cyclic <- function(data, init_vals, n_chains, iter, cores){

	### Fit the model
	model_fit <- rstan::stan(model_code = cyclic_model, 
		data = data, 
		init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)

	return(model_fit)
}

#' Run tensor model
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
spi_tensor <- function(data, init_vals, n_chains, iter, cores){

	### Fit the model
	model_fit <- rstan::stan(model_code = tensor_model, 
		data = data, 
		init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)

	return(model_fit)
}




