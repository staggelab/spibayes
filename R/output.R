
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

	x_mean <- cbind(1, 
		PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_pos)
		)
	x_disp <- cbind(1, 
		PredictMat(model_fit$input$model$gamma$smooth[[2]],newdata_pos)
	)

	x_theta <- cbind(1, 
		PredictMat(model_fit$input$model$theta$smooth[[1]],newdata_all)
	)

	} else if (model_fit$input$type == "tensor") {

		var_list <- c("b_0_mean", "b_0_disp", "b_0_theta", "b_mean_jdate", "b_mean_year", "b_mean_tensor", "b_disp_jdate", "b_disp_year", "b_disp_tensor","b_theta_jdate", "b_theta_year", "b_theta_tensor", "lambda_mean", "lambda_disp", "lambda_theta")

	x_mean <- cbind(1, 
		PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_pos),
		PredictMat(model_fit$input$model$gamma$smooth[[2]],newdata_pos),
		PredictMat(model_fit$input$model$gamma$smooth[[3]],newdata_pos)
		)
	x_disp <- cbind(1, 
		PredictMat(model_fit$input$model$gamma$smooth[[4]],newdata_pos),
		PredictMat(model_fit$input$model$gamma$smooth[[5]],newdata_pos),
		PredictMat(model_fit$input$model$gamma$smooth[[6]],newdata_pos)
	)

	x_theta <- cbind(1, 
		PredictMat(model_fit$input$model$gamma$smooth[[1]],newdata_all),
		PredictMat(model_fit$input$model$theta$smooth[[2]],newdata_all),
		PredictMat(model_fit$input$model$theta$smooth[[3]],newdata_all)
	)


	}

	### Separate first by init vs fit
	### Then by type of model first

	if("fit_params" %in% names(model_fit) == TRUE){
	if(model_fit$fit_params$engine == "optimize") {
		model_read <- cmdstanr::read_cmdstan_csv(model_fit$model_fit$output_files(), variables = var_list)

		param_est <- as.data.frame(model_read$point_estimates) %>%
			mutate(.chain == 0, .iteration = 1, .draw = 1)

		b_0_mean <- param_est %>% select(starts_with("b_0_mean"))
		b_mean <- param_est %>% select(starts_with("b_mean_"))
		b_mean <- as.matrix(cbind(b_0_mean, b_mean))
		mean_est <- x_mean %*% t(b_mean)
		mean_est <- exp(mean_est)

		b_0_disp <- param_est %>% select(starts_with("b_0_disp"))
		b_disp <- param_est %>% select(starts_with("b_disp_"))
		b_disp <- as.matrix(cbind(b_0_disp, b_disp))
		disp_est <- x_disp %*% t(b_disp)
		shape_est <- 1/ exp(-7 + log(1 + exp(disp_est)))

		b_0_theta <- param_est %>% select(starts_with("b_0_theta"))
		b_theta <- param_est %>% select(starts_with("b_theta_"))
		b_theta <- as.matrix(cbind(b_0_theta, b_theta))
		theta_est <- x_theta %*% t(b_theta)
		theta_est <- exp(theta_est)/(1+ exp(theta_est))
	} else if (model_fit$fit_params$engine == "variational") {
	
	} else if (model_fit$fit_params$engine == "sample") {
	
		param_est <- model_fit$model_fit$draws(var_list) %>%
			posterior::as_draws_df() %>%
			as.data.frame()

		b_0_mean <- param_est %>% select(starts_with("b_0_mean"))
		b_mean <- param_est %>% select(starts_with("b_mean_"))
		b_mean <- as.matrix(cbind(b_0_mean, b_mean))
		mean_est <- x_mean %*% t(b_mean)
		mean_est <- exp(mean_est)

		b_0_disp <- param_est %>% select(starts_with("b_0_disp"))
		b_disp <- param_est %>% select(starts_with("b_disp_"))
		b_disp <- as.matrix(cbind(b_0_disp, b_disp))
		disp_est <- x_disp %*% t(b_disp)
		shape_est <- 1/ exp(-7 + log(1 + exp(disp_est)))

		b_0_theta <- param_est %>% select(starts_with("b_0_theta"))
		b_theta <- param_est %>% select(starts_with("b_theta_"))
		b_theta <- as.matrix(cbind(b_0_theta, b_theta))
		theta_est <- x_theta %*% t(b_theta)
		theta_est <- exp(theta_est)/(1+ exp(theta_est))

	}
	} else {

		b_0_mean <- model_fit$input$b_0_init$mean
		b_mean <- c(unlist(model_fit$input$b_init$mean))
		b_mean <- t(as.matrix(c(b_0_mean, b_mean)))
		mean_est <- x_mean %*% t(b_mean)
		mean_est <- exp(mean_est)

		b_0_disp <- model_fit$input$b_0_init$disp
		b_disp <- c(unlist(model_fit$input$b_init$disp))
		b_disp <- t(as.matrix(c(b_0_disp, b_disp)))
		disp_est <- x_disp %*% t(b_disp)
		shape_est <- 1/ exp(-7 + log(1 + exp(disp_est)))

		b_0_theta <- model_fit$input$b_0_init$theta
		b_theta <- c(unlist(model_fit$input$b_init$theta))
		b_theta <- t(as.matrix(c(b_0_theta, b_theta)))
		theta_est <- x_theta %*% t(b_theta)
		theta_est <- exp(theta_est)/(1+ exp(theta_est))
	}

	for(j in seq(1, dim(param_est)[1])){
		param_gamma_j <- newdata_pos %>%
			mutate(chain = param_est$.chain[j], 
				iteration = param_est$.iteration[j],
				draw = param_est$.draw[j]) %>%
			mutate(mean = mean_est[,j], shape = shape_est[,j])

		param_theta_j <- newdata_all %>%
			mutate(chain = param_est$.chain[j], 
				iteration = param_est$.iteration[j],
				draw = param_est$.draw[j]) %>%
			mutate(theta = theta_est[,j])

		if(j == 1){
			param_gamma <- param_gamma_j
			param_theta <- param_theta_j
		} else {
			param_gamma <- bind_rows(param_gamma, param_gamma_j)
			param_theta <- bind_rows(param_theta, param_theta_j)
		}
	}

	param_gamma <- param_gamma %>%
		mutate(disp = 1/shape) %>%
		mutate(scale = mean/shape) %>%
		mutate(rate = 1/scale) %>%
		relocate(chain, iteration, draw, .after = rate)

	output_list <- list(gamma = param_gamma, theta = param_theta)

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
	

