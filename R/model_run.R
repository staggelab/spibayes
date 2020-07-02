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
spi_fit<- function(spi_input, n_chains=1, iter=1000, cores = 1, lambda_year = "free"){
	require(cmdstanr)

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
		b_0_scale_prior=spi_input$b_0$scale,
		lambda_mean_prior = spi_input$lambda_prior$mean,
		lambda_scale_prior = spi_input$lambda_prior$scale
	)
	
	if (spi_input$type == "cyclic"){
		data_fitting$b_mean_prior = spi_input$b_init$mean
		data_fitting$b_scale_prior = spi_input$b_init$scale

		data_fitting$lambda_mean_prior <- c(5, 5/spi_input$lambda_init$mean)
 		data_fitting$lambda_scale_prior <- c(5, 5/spi_input$lambda_init$scale)
	}

#	if (spi_input$type == "cyclic"){
#		data_fitting$lambda_mean_prior <- c(0.5, data_fitting$lambda_mean_prior)
#		data_fitting$lambda_scale_prior <- c(0.5, data_fitting$lambda_scale_prior)
#	} else if (spi_input$type == "tensor"){
#		data_fitting$lambda_mean_prior <- c(rep(0.5,2), data_fitting$lambda_mean_prior)
#		data_fitting$lambda_scale_prior <- c(rep(0.5,2), data_fitting$lambda_scale_prior)
#	}

	### Create the initial values
	init_vals <- list(list(
		b_0_mean = spi_input$b_0$mean[1], 
		b_0_scale = spi_input$b_0$scale[1], 
		b_mean = spi_input$b_init$mean, 
		b_scale = spi_input$b_init$scale, 
		lambda_mean_first = spi_input$lambda_init$mean[1], 
		lambda_scale_first = spi_input$lambda_init$scale[1])
	)

	if (spi_input$type == "tensor" & lambda_year == "free"){
		init_vals[[1]]$lambda_mean_second <- spi_input$lambda_init$mean[2]
		init_vals[[1]]$lambda_scale_second <- spi_input$lambda_init$scale[2]
	}

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
		data_fitting[["S"]] <- as.matrix(Matrix::nearPD(spi_input$s_reparam)$mat)

		### Compile model
		model_fit <- spi_cyclic(data = data_fitting, init_vals = init_vals, n_chains = n_chains, iter = iter, cores = cores)

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
spi_cyclic <- function(data, init_vals, n_chains, iter, cores){

	### Write model to file
	stan_dir <- file.path(cmdstan_path(), "spibayes")
	dir.create(stan_dir, recursive=TRUE, showWarnings = FALSE)
	stan_program <- file.path(stan_dir, "cyclic_model.stan")
	writeLines(cyclic_model, stan_program)

	### Compile model
	mod <- cmdstan_model(stan_program)

	### Fit the model
	model_fit <- mod$sample(
  		data = data,
		  init = init_vals,
		  seed = 123,
		  chains = n_chains,
		  #output_dir = "/media/data/Documents/work_folder/projects_research/code/spi_paper/output/", 
		  iter_sampling = round(iter * 0.6), 
		  iter_warmup = round(iter * 0.4), 
		  save_warmup = TRUE
	)

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
		init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)
	} else if (lambda_year == "fixed"){
	model_fit <- rstan::stan(model_code = tensor_model_lambdafixed, 
		data = data, 
		init = init_vals,
		iter = iter, 
		chains = n_chains,
		cores = cores, 
		verbose = FALSE)
	}

	return(model_fit)
}




