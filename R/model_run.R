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
spi_fit<- function(spi_input, n_chains=1, iter=1000, cores = 1, lambda_year = "free", diagonalize = FALSE){
	require(rstan)

	### Prepare for splitting
	data_extract <- spi_input$data

	data_extract <- data_extract %>%
		mutate(id = seq(1, dim(data_extract)[1])) %>%
		mutate(zero = c(precip == 0))

	id_pos <- data_extract$id[data_extract$zero == FALSE]

	### Extract precipitation
	y <-  spi_input$data$precip
	y_zero <- spi_input$data$zero

	if (spi_input$type == "cyclic" & diagonalize == TRUE){
		X <- spi_input$x_diag
	} else {
		X <- spi_input$x_reparam
	}

	### Extract important values from the input object
	y_pos <- y[id_pos]
	X_pos <- X[id_pos,]

	### Calculate some necesary dimensions
	N <- length(y)
	N_pos <- dim(X_pos)[1]
	basis_dim <- dim(X)[2]
	
	### Error check that the basis and y values are of the same length
	assertthat::assert_that(N == dim(X)[1])

	### Create the data to send to stan model
	data_fitting <- list(N = N, 
		N_pos = N_pos,
		basis_dim = basis_dim, 
		y_pos= y_pos, 
		y_zero = y_zero,
		X = X, 
 		X_pos = X_pos,
		b_0_shape_prior=spi_input$b_0$shape, 
		b_0_rate_prior=spi_input$b_0$rate,
		b_0_theta_prior=spi_input$b_0$theta,
		rho_shape_prior = spi_input$rho_prior$shape,
		rho_rate_prior = spi_input$rho_prior$rate, 
		rho_theta_prior = spi_input$rho_prior$theta
	)

#	if (spi_input$type == "cyclic"){
#		data_fitting$lambda_mean_prior <- c(0.5, data_fitting$lambda_mean_prior)
#		data_fitting$lambda_scale_prior <- c(0.5, data_fitting$lambda_scale_prior)
#	} else if (spi_input$type == "tensor"){
#		data_fitting$lambda_mean_prior <- c(rep(0.5,2), data_fitting$lambda_mean_prior)
#		data_fitting$lambda_scale_prior <- c(rep(0.5,2), data_fitting$lambda_scale_prior)
#	}

	### Add some variance around the b_0
	data_fitting$b_0_shape_prior[[2]] <- 100 * 	data_fitting$b_0_shape_prior[[2]] 
	data_fitting$b_0_rate_prior[[2]] <- 100 * 	data_fitting$b_0_rate_prior[[2]] 
	data_fitting$b_0_theta_prior[[2]] <- 100 * 	data_fitting$b_0_theta_prior[[2]] 

	### Create the initial values
	init_vals <- list(list(
		b_0_shape = spi_input$b_0$shape[1], 
		b_0_rate = spi_input$b_0$rate[1], 
		b_0_theta = spi_input$b_0$theta[1], 
		b_shape = spi_input$b_init$shape, 
		b_rate = spi_input$b_init$rate, 
		b_theta = spi_input$b_init$theta
	))

	if (spi_input$type == "cyclic"){
		init_vals[[1]]$rho_shape <- spi_input$rho_init$shape[1]
		init_vals[[1]]$rho_rate <-  spi_input$rho_init$rate[1]
		init_vals[[1]]$rho_theta <-  spi_input$rho_init$theta[1]
	}

	### Convert initial betas to diagonalized reparameterization
	if (spi_input$type == "cyclic" & diagonalize == TRUE){
		p_inv <- ginv(spi_input$p)
		init_vals[[1]]$b_shape <- c(p_inv %*% init_vals[[1]]$b_shape)
		init_vals[[1]]$b_rate <- c(p_inv %*% init_vals[[1]]$b_rate)
		init_vals[[1]]$b_theta <- c(p_inv %*% init_vals[[1]]$b_theta)
	}

	if (spi_input$type == "tensor"){
		init_vals[[1]]$rho_shape <- spi_input$rho_init$shape
		init_vals[[1]]$rho_rate <-  spi_input$rho_init$rate
		init_vals[[1]]$rho_theta <-  spi_input$rho_init$theta
	}


### Come back here
#	if (spi_input$type == "tensor"){
#		init_vals[[1]]$rho_shape[1] <- spi_input$rho_init$shapemean[1]
#		init_vals[[1]]$rho_shape[1] <- spi_input$lambda_init$mean[1]

#		init_vals[[1]]$lambda_scale_first <-  spi_input$lambda_init$scale[1]
#	}

#	if (spi_input$type == "tensor" & lambda_year == "free"){
#		init_vals[[1]]$lambda_mean_second <- spi_input$lambda_init$mean[2]
#		init_vals[[1]]$lambda_scale_second <- spi_input$lambda_init$scale[2]
#	}

	### Loop through the initial values if we have more than one chain
	### Vary lambda values 2 orders of magnitude smaller or larger
	if(n_chains > 1) {
		for (j in seq(2, n_chains)) {
			init_vals[[j]] <- init_vals[[1]]
			init_vals[[j]]$lambda_shape <- lambda_shape_init*(10^runif(1, -2,2))
			init_vals[[j]]$lambda_rate <- lambda_rate_init*(10^runif(1, -2,2))
		}	
	}

	### Run the cyclic model if cyclic, tensor if tensor
	if (spi_input$type == "cyclic"){

		### Insert the S penaly matrix, which is different shape for cyclic than tensor
		data_fitting[["S"]] <- as.matrix(Matrix::nearPD(spi_input$s_reparam)$mat)

		### Run model
		model_fit <- spi_cyclic(data = data_fitting, init_vals = init_vals, n_chains = n_chains, iter = iter, cores = cores, diagonalize = diagonalize)

	} else if (spi_input$type == "tensor"){
		
		### Insert the 2 S penaly matrices
		### Use the nearPD command to make sure they are positive definite
		data_fitting[["S_1"]] <- as.matrix(Matrix::nearPD(spi_input$s_reparam[[1]])$mat)
		data_fitting[["S_2"]] <- as.matrix(Matrix::nearPD(spi_input$s_reparam[[2]])$mat)

		### Run model
		model_fit <- spi_tensor(data = data_fitting, init_vals = init_vals, n_chains = n_chains, iter = iter, cores = cores, lambda_year = lambda_year)
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
spi_cyclic <- function(data, init_vals, n_chains, iter, cores, diagonalize){


	if (diagonalize == TRUE){
	model_fit <- rstan::stan(model_code = cyclic_model_diag, 
		data = data, 
		init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)

	} else {
	model_fit <- rstan::stan(model_code = cyclic_model, 
		data = data, 
		init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)

	}

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
spi_tensor <- function(data, init_vals, n_chains, iter, cores,  lambda_year){

	if(lambda_year == "free"){
	### Fit the model
	model_fit <- rstan::stan(model_code = tensor_model, 
		data = data, 
		#init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)
	} else if (lambda_year == "fixed"){
	model_fit <- rstan::stan(model_code = tensor_model_lambdafixed, 
		data = data, 
	#	init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)
	}

	return(model_fit)
}




