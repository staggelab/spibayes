
#' Extract model output parameters
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param model_fit Dataframe with the underlying data. Columns must include variable names
#' @param newdata List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param type Simple gives back only the parameters, Marginal gives each marginal
#' @return A matrix of the infile
#' @export
predict_vals <- function(model_fit, newdata = NULL){

	if (is.null(newdata)){
		newdata_pos <- model_fit$input$data_pos
		newdata_all <- model_fit$input$data_all
	} else {
		newdata_pos <- newdata
		newdata_all <- newdata
	}

	### Whether cyclic or tensor
	if(model_fit$input$type == "cyclic") {

		var_list <- c("b_0_mean", "b_0_disp", "b_0_theta", "b_mean_jdate", "b_disp_jdate", "b_theta_jdate", "lambda_mean", "lambda_disp", "lambda_theta")

	x_mean_jdate <- PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_pos)
	x_disp_jdate <- PredictMat(model_fit$input$model$gamma$smooth[[2]],newdata_pos)
	x_theta_jdate <- PredictMat(model_fit$input$model$theta$smooth[[1]],newdata_all)

	} else if (model_fit$input$type == "tensor") {

		var_list <- c("b_0_mean", "b_0_disp", "b_0_theta", "b_mean_jdate", "b_mean_year", "b_mean_tensor", "b_disp_jdate", "b_disp_year", "b_disp_tensor","b_theta_jdate", "b_theta_year", "b_theta_tensor", "lambda_mean", "lambda_disp", "lambda_theta")

	x_mean_jdate <- PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_pos)
	x_mean_year <- PredictMat(model_fit$input$model$gamma$smooth[[2]],newdata_pos)
	x_mean_tensor <- PredictMat(model_fit$input$model$gamma$smooth[[3]],newdata_pos)

	x_disp_jdate <- PredictMat(model_fit$input$model$gamma$smooth[[4]],newdata_pos)
	x_disp_year <- PredictMat(model_fit$input$model$gamma$smooth[[5]],newdata_pos)
	x_disp_tensor <- PredictMat(model_fit$input$model$gamma$smooth[[6]],newdata_pos)

	x_theta_jdate <- PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_all)
	x_theta_year <- PredictMat(model_fit$input$model$gamma$smooth[[2]],newdata_all)
	x_theta_tensor <- PredictMat(model_fit$input$model$gamma$smooth[[3]],newdata_all)

	#x_mean <- cbind(1, 
	#	PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_pos),
	#	PredictMat(model_fit$input$model$gamma$smooth[[2]],newdata_pos),
	#	PredictMat(model_fit$input$model$gamma$smooth[[3]],newdata_pos)
	#	)

	#x_disp <- cbind(1, 
	#	PredictMat(model_fit$input$model$gamma$smooth[[4]],newdata_pos),
	#	PredictMat(model_fit$input$model$gamma$smooth[[5]],newdata_pos),
	#	PredictMat(model_fit$input$model$gamma$smooth[[6]],newdata_pos)
	#)

	#x_theta <- cbind(1, 
	#	PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_all),
	#	PredictMat(model_fit$input$model$theta$smooth[[2]],newdata_all),
	#	PredictMat(model_fit$input$model$theta$smooth[[3]],newdata_all)
	#)

	}

	### Separate first by init vs fit
	### Then by type of model first

	if("fit_params" %in% names(model_fit) == TRUE){
	if(model_fit$fit_params$engine == "optimize") {
		model_read <- cmdstanr::read_cmdstan_csv(model_fit$model_fit$output_files(), variables = var_list)

		param_est <- as.data.frame(model_read$point_estimates) %>%
			mutate(.chain = 0, .iteration = 1, .draw = 1)

		b_0_mean <- param_est %>% select(starts_with("b_0_mean"))
		b_mean_jdate <- param_est %>% select(starts_with("b_mean_jdate"))
		mean_0 <- t(b_0_mean)
		mean_jdate <- x_mean_jdate %*% t(b_mean_jdate)
		sigma_mean <- matrix(rep(NA, length(mean_jdate)), dim(mean_jdate))

		b_0_disp <- param_est %>% select(starts_with("b_0_disp"))
		b_disp_jdate <- param_est %>% select(starts_with("b_disp_"))
		disp_0 <- t(b_0_disp)
		disp_jdate <- x_disp_jdate %*% t(b_disp_jdate)
		sigma_disp <- matrix(rep(NA, length(disp_jdate)), dim(disp_jdate))

		b_0_theta <- param_est %>% select(starts_with("b_0_theta"))
		b_theta_jdate <- param_est %>% select(starts_with("b_theta_"))
		theta_0 <- t(b_0_theta)
		theta_jdate <- x_theta_jdate %*% t(b_theta_jdate)
		sigma_theta <- matrix(rep(NA, length(theta_jdate)), dim(theta_jdate))

		if (model_fit$input$type == "tensor") {
			b_mean_year <- param_est %>% select(starts_with("b_mean_year"))
			b_mean_tensor <- param_est %>% select(starts_with("b_mean_tensor"))	
			mean_year <- x_mean_year %*% t(b_mean_year)
			mean_tensor <- x_mean_tensor %*% t(b_mean_tensor)	

			b_disp_year <- param_est %>% select(starts_with("b_disp_year"))
			b_disp_tensor <- param_est %>% select(starts_with("b_disp_tensor"))	
			disp_year <- x_disp_year %*% t(b_disp_year)
			disp_tensor <- x_disp_tensor %*% t(b_disp_tensor)	

			b_theta_year <- param_est %>% select(starts_with("b_theta_year"))
			b_theta_tensor <- param_est %>% select(starts_with("b_theta_tensor"))	
			theta_year <- x_theta_year %*% t(b_theta_year)
			theta_tensor <- x_theta_tensor %*% t(b_theta_tensor)	
		}

	} else if (model_fit$fit_params$engine == "variational") {
	
	} else if (model_fit$fit_params$engine == "sample") {
		var_list <- c(var_list, c("sigma_mean", "sigma_disp", "sigma_theta"))

		param_est <- model_fit$model_fit$draws(var_list) %>%
			posterior::as_draws_df() %>%
			as.data.frame()

		b_0_mean <- param_est %>% select(starts_with("b_0_mean"))
		b_mean_jdate <- param_est %>% select(starts_with("b_mean_jdate"))
		mean_0 <- t(b_0_mean)
		mean_jdate <- x_mean_jdate %*% t(b_mean_jdate)
		sigma_mean <- t(param_est %>% select("sigma_mean"))

		b_0_disp <- param_est %>% select(starts_with("b_0_disp"))
		b_disp_jdate <- param_est %>% select(starts_with("b_disp_"))
		disp_0 <- t(b_0_disp)
		disp_jdate <- x_disp_jdate %*% t(b_disp_jdate)
		sigma_disp <- t(param_est %>% select("sigma_disp"))

		b_0_theta <- param_est %>% select(starts_with("b_0_theta"))
		b_theta_jdate <- param_est %>% select(starts_with("b_theta_"))
		theta_0 <- t(b_0_theta)
		theta_jdate <- x_theta_jdate %*% t(b_theta_jdate)
		sigma_theta <- t(param_est %>% select("sigma_theta"))

		if (model_fit$input$type == "tensor") {
			b_mean_year <- param_est %>% select(starts_with("b_mean_year"))
			b_mean_tensor <- param_est %>% select(starts_with("b_mean_tensor"))	
			mean_year <- x_mean_year %*% t(b_mean_year)
			mean_tensor <- x_mean_tensor %*% t(b_mean_tensor)	

			b_disp_year <- param_est %>% select(starts_with("b_disp_year"))
			b_disp_tensor <- param_est %>% select(starts_with("b_disp_tensor"))	
			disp_year <- x_disp_year %*% t(b_disp_year)
			disp_tensor <- x_disp_tensor %*% t(b_disp_tensor)	

			b_theta_year <- param_est %>% select(starts_with("b_theta_year"))
			b_theta_tensor <- param_est %>% select(starts_with("b_theta_tensor"))	
			theta_year <- x_theta_year %*% t(b_theta_year)
			theta_tensor <- x_theta_tensor %*% t(b_theta_tensor)	
		}

	}
	### This uses the initial estimates from preprocessing
	} else {
		var_list <- c(var_list, c("sigma_mean", "sigma_disp", "sigma_theta"))
		param_est <- data.frame(.chain = 0, .iteration = 1, .draw = 1)

		b_0_mean <- model_fit$input$b_0_init$mean
		b_mean_jdate <- matrix(c(unlist(model_fit$input$b_init$mean$jdate)), 1, dim(x_mean_jdate)[2])
		mean_0 <- t(b_0_mean)
		mean_jdate <- x_mean_jdate %*% t(b_mean_jdate)
		sigma_mean <- matrix(rep(model_fit$input$sigma_init$mean, length(mean_jdate)), dim(mean_jdate))

		b_0_disp <- model_fit$input$b_0_init$disp
		b_disp_jdate <- matrix(c(unlist(model_fit$input$b_init$disp$jdate)), 1, dim(x_disp_jdate)[2])
		disp_0 <- t(b_0_disp)
		disp_jdate <- x_disp_jdate %*% t(b_disp_jdate)
		sigma_disp <- matrix(rep(model_fit$input$sigma_init$disp, length(disp_jdate)), dim(disp_jdate))

		b_0_theta <- model_fit$input$b_0_init$theta
		b_theta_jdate <- matrix(c(unlist(model_fit$input$b_init$theta$jdate)), 1, dim(x_theta_jdate)[2])
		theta_0 <- t(b_0_theta)
		theta_jdate <- x_theta_jdate %*% t(b_theta_jdate)
		sigma_theta <- matrix(rep(model_fit$input$sigma_init$theta, length(theta_jdate)), dim(theta_jdate))


		if (model_fit$input$type == "tensor") {
			b_mean_year <-  matrix(c(unlist(model_fit$input$b_init$mean$year)), 1, dim(x_mean_year)[2])
			b_mean_tensor <-  matrix(c(unlist(model_fit$input$b_init$mean$tensor)), 1, dim(x_mean_tensor)[2])
			mean_year <- x_mean_year %*% t(b_mean_year)
			mean_tensor <- x_mean_tensor %*% t(b_mean_tensor)	

			b_disp_year <-  matrix(c(unlist(model_fit$input$b_init$disp$year)), 1, dim(x_disp_year)[2])
			b_disp_tensor <-  matrix(c(unlist(model_fit$input$b_init$disp$tensor)), 1, dim(x_disp_tensor)[2])
			disp_year <- x_disp_year %*% t(b_disp_year)
			disp_tensor <- x_disp_tensor %*% t(b_disp_tensor)	

			b_theta_year <-  matrix(c(unlist(model_fit$input$b_init$theta$year)), 1, dim(x_theta_year)[2])
			b_theta_tensor <-  matrix(c(unlist(model_fit$input$b_init$theta$tensor)), 1, dim(x_theta_tensor)[2])
			theta_year <- x_theta_year %*% t(b_theta_year)
			theta_tensor <- x_theta_tensor %*% t(b_theta_tensor)	
		}


	}

	for(j in seq(1, dim(param_est)[1])){
		marginal_j <- newdata_pos %>%
			mutate(chain = param_est$.chain[j], 
				iteration = param_est$.iteration[j],
				draw = param_est$.draw[j]) 

		### Add marginal columns for jdate
		marginal_mean_j <- marginal_j %>%
			mutate(mean_0 = mean_0[,j]) %>%
			mutate(mean_jdate = mean_jdate[,j], mean_year = NA, mean_tensor = NA) %>%
			mutate(sigma_mean = sigma_mean[,j])
		
		marginal_disp_j <- marginal_j %>%		
			mutate(disp_0 = disp_0[,j]) %>%
			mutate(disp_jdate = disp_jdate[,j], disp_year = NA, disp_tensor = NA) %>%
			mutate(sigma_disp = sigma_disp[,j])
	
		marginal_theta_j <- newdata_all %>%
			mutate(chain = param_est$.chain[j], 
				iteration = param_est$.iteration[j],
				draw = param_est$.draw[j]) %>%		
			mutate(theta_0 = theta_0[,j]) %>%
			mutate(theta_jdate = theta_jdate[,j], theta_year = NA, theta_tensor = NA) %>%
			mutate(sigma_theta = sigma_theta[,j])

		### If tensor add these columns	
		if(model_fit$input$type == "tensor") {
			marginal_mean_j <- marginal_mean_j %>%
				mutate(mean_year = mean_year[,j], mean_tensor = mean_tensor[,j])

			marginal_disp_j <- marginal_disp_j %>%
				mutate(disp_year = disp_year[,j], disp_tensor = disp_tensor[,j])

			marginal_theta_j <- marginal_theta_j %>%
				mutate(theta_year = theta_year[,j], theta_tensor = theta_tensor[,j])
		}

	marginal_mean_j <- marginal_mean_j %>%
			mutate(mean = rowSums(select(., starts_with("mean_")), na.rm=TRUE)) %>% 
			relocate(sigma_mean, .after = last_col())

	marginal_disp_j <- marginal_disp_j %>%
			mutate(disp = rowSums(select(., starts_with("disp_")), na.rm=TRUE)) %>% 
			relocate(sigma_disp, .after = last_col())

	marginal_theta_j <- marginal_theta_j %>%
			mutate(theta = rowSums(select(., starts_with("theta_")), na.rm=TRUE)) %>% 
			relocate(sigma_theta, .after = last_col())

		if(j == 1){
			marginal_mean <- marginal_mean_j
			marginal_disp <- marginal_disp_j
			marginal_theta <- marginal_theta_j
		} else {
			marginal_mean <- bind_rows(marginal_mean, marginal_mean_j)
			marginal_disp <- bind_rows(marginal_disp, marginal_disp_j)
			marginal_theta <- bind_rows(marginal_theta, marginal_theta_j)
		}
	}

	param_gamma <-	marginal_mean %>%
		select(-c(mean_0,mean_jdate, mean_year, mean_tensor,mean, sigma_mean)) %>%
		mutate(mean = exp(marginal_mean$mean)) %>%
		mutate(shape = 1/ exp(-7 + log(1 + exp(marginal_disp$disp)))) %>%
		mutate(disp = 1/shape) %>%
		mutate(scale = mean/shape) %>%
		mutate(rate = 1/scale) %>%
		relocate(chain, iteration, draw, .after = rate)

	param_theta <- marginal_theta %>%
		select(-c(theta_0,theta_jdate, theta_year, theta_tensor, sigma_theta)) %>%	
		mutate(theta = 	exp(theta)/(1+ exp(theta))) %>%
		relocate(chain, iteration, draw, .after = theta)


	output_list <- list(estimate = list(gamma = param_gamma, theta = param_theta), marginal = list(mean = marginal_mean, disp = marginal_disp, theta = marginal_theta))

	return(output_list)

}
	






#' Extract model output parameters
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
extract_params <- function(model_fit, var_list = NULL){

	### Whether cyclic or tensor
	if(model_fit$input$type == "cyclic" & is.null(var_list) ) {

		var_list <- c("b_0_mean", "b_0_disp", "b_0_theta", "b_mean_jdate", "b_disp_jdate", "b_theta_jdate", "lambda_mean", "lambda_disp", "lambda_theta")

	} else if (model_fit$input$type == "tensor" & is.null(var_list)) {

		var_list <- c("b_0_mean", "b_0_disp", "b_0_theta", "b_mean_jdate", "b_mean_year", "b_mean_tensor", "b_disp_jdate", "b_disp_year", "b_disp_tensor","b_theta_jdate", "b_theta_year", "b_theta_tensor", "lambda_mean", "lambda_disp", "lambda_theta")

	} else {
		var_list <- var_list
	}

	### Separate by type of model first
	if(model_fit$fit_params$engine == "optimize") {
		model_read <- cmdstanr::read_cmdstan_csv(model_fit$model_fit$output_files(), variables = var_list)

		param_est <- as.data.frame(model_read$point_estimates)

	} else if (model_fit$fit_params$engine == "variational") {
	
	
	} else if (model_fit$fit_params$engine == "sample") {
	

	}

	return(param_est)

}
	

