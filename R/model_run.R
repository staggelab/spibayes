#' Run model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @param engine  options are sample (full MCMC), variational (variational approximation of posterior),  optimize (Stan's LBFGS algorithm), 
#' @import cmdstanr
#' @return A matrix of the infile
#' @export
model_fit <- function(spi_input, n_chains=1, iter=1000, cores = 1, engine = "sample", output_dir = getwd()){
	require(cmdstanr)

	spi_input <- spi_input$input

	### Run the correct model
	if(spi_input$type == "cyclic"){
		model_fit <- cyclic_fit(spi_input = spi_input, n_chains = n_chains, iter = iter, cores = cores, engine = engine, output_dir = output_dir)
	} else if (spi_input$type == "tensor"){
		model_fit <- tensor_fit(spi_input = spi_input, n_chains = n_chains, iter = iter, cores = cores, engine = engine, output_dir = output_dir)
	}

	fit_params <- list(n_chains = n_chains, iter = iter, engine = engine)

	### Add in all the input data 
	output_fit <- list(input = spi_input, model_fit = model_fit, fit_params = fit_params)

	return(output_fit)
}


#' Run Cyclic model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @param engine  options are sample (full MCMC), variational (variational approximation of posterior),  optimize (Stan's LBFGS algorithm), 
#' @import cmdstanr
#' @return A matrix of the infile
#' @export
cyclic_fit <- function(spi_input, n_chains=1, iter=1000, cores = 1, engine = "sample", output_dir = getwd()){
	require(cmdstanr)

	### Extract precipitation and zeros
	y_pos <-  spi_input$data_pos$precip
	y_zero <- spi_input$data_all$zero

	### Calculate some required dimensions
	N <- length(y_zero)
	N_pos <- length(y_pos)

	### Create the data to send to stan model
	data_fitting <- list(
		N = N, 
		N_pos = N_pos,

		basis_dim_mean = spi_input$basis_dim$mean, 
		basis_dim_disp = spi_input$basis_dim$disp, 
		basis_dim_theta = spi_input$basis_dim$theta, 

		y_pos= y_pos, 
		y_zero = y_zero,

		X_mean_jdate = spi_input$x_matrix$mean$jdate, 
		X_disp_jdate = spi_input$x_matrix$disp$jdate, 
		X_theta_jdate = spi_input$x_matrix$theta$jdate,  

		b_0_prior = spi_input$b_0_prior,

		s_mean_jdate = spi_input$s_matrix$mean$jdate, 
		s_disp_jdate = spi_input$s_matrix$disp$jdate, 
		s_theta_jdate = spi_input$s_matrix$theta$jdate, 

		lambda_mean_prior = spi_input$lambda_prior$mean,
		lambda_disp_prior = spi_input$lambda_prior$disp,
		lambda_theta_prior = spi_input$lambda_prior$theta, 

		sigma_mean_prior = spi_input$sigma_prior$mean,
		sigma_disp_prior = spi_input$sigma_prior$disp,
		sigma_theta_prior = spi_input$sigma_prior$theta
	)

	### Ensure the penalties positive definite
	penalty_index <- which(startsWith(names(data_fitting), "s_"))
	for (j in penalty_index){
		data_fitting[[j]] <- as.matrix(Matrix::nearPD(data_fitting[[j]])$mat)
	}

	### Reparse for init vals
	init_vals <- list(
		b_0_mean = spi_input$b_0_init$mean,
		b_0_disp = spi_input$b_0_init$disp,
		b_0_theta = spi_input$b_0_init$theta,
		b_mean_jdate = spi_input$b_init$mean$jdate,
		b_disp_jdate = spi_input$b_init$disp$jdate,
		b_theta_jdate = spi_input$b_init$theta$jdate,
		lambda_mean = unlist(spi_input$lambda_init$mean),
		lambda_disp = unlist(spi_input$lambda_init$disp),
		lambda_theta = unlist(spi_input$lambda_init$theta),
		sigma_mean = unlist(spi_input$sigma_init$mean),
		sigma_disp = unlist(spi_input$sigma_init$disp),
		sigma_theta = unlist(spi_input$sigma_init$theta)
	)
	###Initial values must be one list for each chain including the first
	init_vals <- list(init_vals)

	### Set up for multiple chains
	if(n_chains > 1){
	for (j in seq(2,n_chains)){
		init_vals[[j]] <- init_vals[[1]]
	
		### Only tweak the year penalties
		init_vals[[j]]$lambda_mean[1] <- init_vals[[j]]$lambda_mean[1] * exp(runif(1, -3,3))
		init_vals[[j]]$lambda_disp[1] <- init_vals[[j]]$lambda_disp[1] * exp(runif(1, -3,3))
		init_vals[[j]]$lambda_theta[1] <- init_vals[[j]]$lambda_theta[1] * exp(runif(1, -3,3))
	}
	}

	### Compile the model
	if (engine == "optimize"){
		cyclic_mod <- cmdstan_model(system.file("models/cyclic_ti_optimize.stan", package = "spibayes"))
	} else {
		cyclic_mod <- cmdstan_model(system.file("models/cyclic_ti.stan", package = "spibayes"))
	}

	#	cyclic_mod <- cmdstan_model(system.file("models/cyclic_ti.stan", package = "spibayes"))

	### Run the model
	if(engine == "sample"){
		#avail_cores <- parallel::detectCores()
		#options(mc.cores = avail_cores)
			
		model_fit <- cyclic_mod$sample(
  			data = data_fitting,
			#iter = iter,
			chains = n_chains,
			#parallel_chains = cores,			
			iter_warmup = iter*0.3,
			iter_sampling = iter *0.7, 
  			init = init_vals,
  			refresh = 10, 
			save_warmup=TRUE,
			output_dir = output_dir
		)
	} else if (engine == "variational"){
		model_fit <- cyclic_mod$variational(
  			data = data_fitting,
			iter = iter,
  			init = init_vals,
  			refresh = 10,
			output_dir = output_dir
		)
	} else if (engine == "optimize"){
		model_fit <- cyclic_mod$optimize(
  			data = data_fitting,
			iter = iter,
  			init = init_vals,
  			refresh = 10,
			output_dir = output_dir
		)
	}

	return(model_fit)
}







#' Run Tensor model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @param engine  options are sample (full MCMC), variational (variational approximation of posterior),  optimize (Stan's LBFGS algorithm), 
#' @import cmdstanr
#' @return A matrix of the infile
#' @export
tensor_fit <- function(spi_input, n_chains=1, iter=1000, cores = 1, engine = "sample", output_dir = getwd()){
	require(cmdstanr)

	### Extract precipitation and zeros
	y_pos <-  spi_input$data_pos$precip
	y_zero <- spi_input$data_all$zero

	### Calculate some required dimensions
	N <- length(y_zero)
	N_pos <- length(y_pos)

	### Create the data to send to stan model
	data_fitting <- list(
		N = N, 
		N_pos = N_pos,

		basis_dim_mean = spi_input$basis_dim$mean, 
		basis_dim_disp = spi_input$basis_dim$disp, 
		basis_dim_theta = spi_input$basis_dim$theta, 

		y_pos= y_pos, 
		y_zero = y_zero,

		X_mean_jdate = spi_input$x_matrix$mean$jdate, 
		X_mean_year = spi_input$x_matrix$mean$year, 
		X_mean_tensor = spi_input$x_matrix$mean$tensor, 

		X_disp_jdate = spi_input$x_matrix$disp$jdate, 
		X_disp_year = spi_input$x_matrix$disp$year, 
		X_disp_tensor = spi_input$x_matrix$disp$tensor, 

		X_theta_jdate = spi_input$x_matrix$theta$jdate,  
		X_theta_year = spi_input$x_matrix$theta$year, 
		X_theta_tensor = spi_input$x_matrix$theta$tensor, 

		b_0_prior = spi_input$b_0_prior,

		s_mean_jdate = spi_input$s_matrix$mean$jdate, 
		s_mean_year = spi_input$s_matrix$mean$year, 
		s_mean_year_double = spi_input$s_matrix$mean$year_double, 
		s_mean_tensor_jdate = spi_input$s_matrix$mean$tensor_jdate, 
		s_mean_tensor_year = spi_input$s_matrix$mean$tensor_year, 

		s_disp_jdate = spi_input$s_matrix$disp$jdate, 
		s_disp_year = spi_input$s_matrix$disp$year, 
		s_disp_year_double = spi_input$s_matrix$disp$year_double, 
		s_disp_tensor_jdate = spi_input$s_matrix$disp$tensor_jdate, 
		s_disp_tensor_year = spi_input$s_matrix$disp$tensor_year, 

		s_theta_jdate = spi_input$s_matrix$theta$jdate, 
		s_theta_year = spi_input$s_matrix$theta$year, 
		s_theta_year_double = spi_input$s_matrix$theta$year_double, 
		s_theta_tensor_jdate = spi_input$s_matrix$theta$tensor_jdate, 
		s_theta_tensor_year = spi_input$s_matrix$theta$tensor_year, 

		lambda_mean_prior = spi_input$lambda_prior$mean,
		lambda_disp_prior = spi_input$lambda_prior$disp,
		lambda_theta_prior = spi_input$lambda_prior$theta, 

		sigma_mean_prior = spi_input$sigma_prior$mean,
		sigma_disp_prior = spi_input$sigma_prior$disp,
		sigma_theta_prior = spi_input$sigma_prior$theta
	)

	### Ensure the penalties positive definite
	penalty_index <- which(startsWith(names(data_fitting), "s_"))
	for (j in penalty_index){
		data_fitting[[j]] <- as.matrix(Matrix::nearPD(data_fitting[[j]])$mat)
	}

	### Reparse for init vals
	init_vals <- list(
		b_0_mean = spi_input$b_0_init$mean,
		b_0_disp = spi_input$b_0_init$disp,
		b_0_theta = spi_input$b_0_init$theta,
		b_mean_jdate = spi_input$b_init$mean$jdate,
		b_mean_year = spi_input$b_init$mean$year,
		b_mean_tensor = spi_input$b_init$mean$tensor,
		b_disp_jdate = spi_input$b_init$disp$jdate,
		b_disp_year = spi_input$b_init$disp$year,
		b_disp_tensor = spi_input$b_init$disp$tensor,
		b_theta_jdate = spi_input$b_init$theta$jdate,
		b_theta_year = spi_input$b_init$theta$year,
		b_theta_tensor = spi_input$b_init$theta$tensor,
		lambda_mean = unlist(spi_input$lambda_init$mean),
		lambda_disp = unlist(spi_input$lambda_init$disp),
		lambda_theta = unlist(spi_input$lambda_init$theta),
		sigma_mean = unlist(spi_input$sigma_init$mean),
		sigma_disp = unlist(spi_input$sigma_init$disp),
		sigma_theta = unlist(spi_input$sigma_init$theta)
	)
	###Initial values must be one list for each chain including the first
	init_vals <- list(init_vals)

	### Set up for multiple chains
	if(n_chains > 1){
	for (j in seq(2,n_chains)){
		init_vals[[j]] <- init_vals[[1]]
	
		### Only tweak the year penalties
		init_vals[[j]]$lambda_mean[c(2,3,5)] <- init_vals[[j]]$lambda_mean[c(2,3,5)] * exp(runif(3, -3,3))
		init_vals[[j]]$lambda_disp[c(2,3,5)] <- init_vals[[j]]$lambda_disp[c(2,3,5)] * exp(runif(3, -3,3))
		init_vals[[j]]$lambda_theta[c(2,3,5)] <- init_vals[[j]]$lambda_theta[c(2,3,5)] * exp(runif(3, -3,3))
	}
	}

	### Compile the model
	if (engine == "optimize"){
		tensor_mod <- cmdstan_model(system.file("models/tensor_ti_optimize.stan", package = "spibayes"))
	} else {
		tensor_mod <- cmdstan_model(system.file("models/tensor_ti.stan", package = "spibayes"))
	}


	### Run the model
	if(engine == "sample"){
		#avail_cores <- parallel::detectCores()
		#options(mc.cores = avail_cores)
			
		model_fit <- tensor_mod$sample(
  			data = data_fitting,
			#iter = iter,
			chains = n_chains,
			#parallel_chains = cores,			
			iter_warmup = iter*0.3,
			iter_sampling = iter *0.7, 
  			init = init_vals,
  			refresh = 10, 
			save_warmup=TRUE,
			output_dir = output_dir
		)
	} else if (engine == "variational"){
		model_fit <- tensor_mod$variational(
  			data = data_fitting,
			iter = iter,
  			init = init_vals,
  			refresh = 10,
			output_dir = output_dir
		)
	} else if (engine == "optimize"){
		model_fit <- tensor_mod$optimize(
  			data = data_fitting,
			iter = iter,
  			init = init_vals,
  			refresh = 10,
			output_dir = output_dir
		)
	}

	return(model_fit)
}




