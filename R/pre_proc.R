#' Create Spline Basis
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
pre_proc <- function(basis){
	if (basis$type == "cyclic"){
		preproc_output <- prior_cyclic(data = basis$data, knot_loc = basis$knot_loc)
	} else if (basis$type == "tensor"){
		
		cyclic_output <- prior_cyclic(data = basis$data, knot_loc = list(jdate = basis$knot_loc$jdate))


		preproc_output <- prior_tensor_fromcyclic(basis = basis, cyclic_output=cyclic_output)
	}

	output_list <- append(basis, preproc_output)
	return(output_list)
}

#' Estimate priors for cyclic spline
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @import dplyr
#' @export
prior_cyclic <- function(data, knot_loc ){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)
	
	### How many years of data to sample from MLE estimate
	if("year" %in% colnames(data)) {
		n_years <- length(unique(data$year))
	} else {
		n_years <- length(unique(lubridate::year(data$date)))
	}

	### Create MLE estimate
	### Fit each day with MLE
	mle_fit <- mle_gamma_fit(jdate = data$jdate, values = data$precip, n = n_years)
	

	#mle_fit <- data %>%
	#	dplyr::group_by(jdate) %>%
	#	dplyr::filter(precip > 0) %>%
	#	dplyr::summarise(shape = fitdistrplus::fitdist(precip, "gamma")[[1]][[1]], rate = fitdistrplus::fitdist(precip, "gamma")[[1]][[2]]) %>%
	#	dplyr::mutate(scale = 1/rate) %>%
	#	dplyr::mutate(mean = scale * shape) %>%
	#	dplyr::mutate(disp = 1/shape) %>%
	#	dplyr::ungroup()

	### Fit a cyclic spline model using mgcv
	shape_model <- mgcv::gam(log(shape) ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=mle_fit$draws, knots = knot_loc, select=FALSE)
	rate_model <- mgcv::gam(log(rate) ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=mle_fit$draws, knots = knot_loc, select=FALSE)
	theta_model <- mgcv::gam(zero ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=data, knots = knot_loc, select=FALSE,family=binomial)

	### Create the prior for the mean intercept
	b_0_shape_prior <- c(summary(shape_model)$p.table[1], summary(shape_model)$p.table[2])
	b_0_rate_prior <- c(summary(rate_model)$p.table[1], summary(rate_model)$p.table[2])
	b_0_theta_prior <- c(summary(theta_model)$p.table[1], summary(theta_model)$p.table[2])

	### Create a vector for intializing the mean
	b_shape_init <- c(coef(shape_model)[2:length(coef(shape_model))])
	b_rate_init <- c(coef(rate_model)[2:length(coef(rate_model))])
	b_theta_init <- c(coef(theta_model)[2:length(coef(theta_model))])

	### Extract smoothing penalty
	lambda_shape_init <- unname(c(shape_model$sp ))
	lambda_rate_init <- unname(c(rate_model$sp))
	lambda_theta_init <- unname(c(theta_model$sp))

	### Calculate alpha the rescaling factor and then rescale
	shape_alpha <- sapply(shape_model$smooth, "[[", "S.scale") / lambda_shape_init
	lambda_shape_init <- c(lambda_shape_init / shape_alpha)

	rate_alpha <- sapply(rate_model$smooth, "[[", "S.scale") / lambda_rate_init
	lambda_rate_init <- c(lambda_rate_init / rate_alpha)

	theta_alpha <- sapply(theta_model$smooth, "[[", "S.scale") / lambda_theta_init
	lambda_theta_init <- c(lambda_theta_init / theta_alpha)

	### Extract penalty matrix
	s_matrix <- shape_model$smooth[[1]]$S[[1]]
	#s_inv <- MASS::ginv(s_matrix)

	### Rescale to match standard devs
	### Might not need to do this
	#shape_samples <- mvrnorm(n=100, b_shape_init, s_inv)
	#rate_samples <- mvrnorm(n=100, b_rate_init, s_inv)
	#theta_samples <- mvrnorm(n=100, b_theta_init, s_inv)

	#lambda_shape <- 1/((sd(b_shape_init)/mean(apply(shape_samples, 1, sd)))^2)
	#lambda_rate <- 1/((sd(b_rate_init)/mean(apply(rate_samples, 1, sd)))^2)
	#lambda_theta <- 1/((sd(b_theta_init)/mean(apply(theta_samples, 1, sd)))^2)

	#mean_testing <- mvrnorm(n=100, b_mean_init, s_inv/lambda_mean)
	#scale_testing <- mvrnorm(n=100, b_scale_init, s_inv/lambda_scale)	
	#theta_testing <- mvrnorm(n=100, b_theta_init, s_inv/lambda_theta_init)
	
	rho_init_shape <- log(lambda_shape_init)
	rho_init_rate <- log(lambda_rate_init)
	rho_init_theta <- log(lambda_theta_init)

	rho_shape_prior <- c(rho_init_shape-3, rho_init_shape + 3)
	rho_rate_prior <- c(rho_init_rate-3, rho_init_rate + 3)
	rho_theta_prior <- c(rho_init_theta-3, rho_init_theta + 3)

	### Create output list and return
	output_list <- list(
		b_0 = list(shape = b_0_shape_prior, rate = b_0_rate_prior,  theta = b_0_theta_prior), 
		b_init = list(shape = b_shape_init, rate = b_rate_init, theta = b_theta_init), 
		lambda_init = list(shape = lambda_shape_init, rate = lambda_rate_init, theta = lambda_theta_init), 
		rho_init = list(shape = rho_init_shape, rate = rho_init_rate, theta = rho_init_theta),
		rho_prior = list(shape = rho_shape_prior, rate = rho_rate_prior, theta = rho_theta_prior)
	)
	return(output_list)

}





#' Estimate priors for tensor spline
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @import dplyr
#' @export
prior_tensor_fromcyclic <- function(basis = basis, cyclic_output=cyclic_output){
	### Extract the number of knots
	knot_loc <- basis$knot_loc
	n_knots <- sapply(knot_loc, FUN = length)
	
	### Add a year column if it doesn't exist
	### It should, but just in case
	data <- basis$data
	if(!("year" %in% colnames(data))) {
		data <- data %>%
			mutate(year = lubridate::year(date))		
	}

	### Remove any leap days for fitting purposes
	data <- data %>%
		filter(jdate <=365)

	### Add in a draw column with the unique years so this can later be joined
	### Base number of samples from MLE off of this
	data <- data %>% 
		mutate(draw = as.numeric(factor(year)))
	n_years <- max(data$draw)


	### Create cyclic basis internally for conversion
	cyclic_basis <- create_cyclic_basis(data = data, knot_loc = list(jdate = basis$knot_loc$jdate))

	### Process the beta init
	b_init_cyclic_shape <- cyclic_output$b_init$shape
	b_init_cyclic_rate <- cyclic_output$b_init$rate
	b_init_cyclic_theta <- cyclic_output$b_init$theta

	### Convert to original basis
	b_init_cyclic_shape_orig <- c(t(cyclic_basis$z %*% b_init_cyclic_shape ))
	b_init_cyclic_rate_orig <- c(t(cyclic_basis$z %*% b_init_cyclic_rate ))
	b_init_cyclic_theta_orig <- c(t(cyclic_basis$z %*% b_init_cyclic_theta ))

	### Repeat to create identical init values through time
	b_init_shape_orig <-  rep(b_init_cyclic_shape_orig, each = length(knot_loc$year))
	b_init_rate_orig <-  rep(b_init_cyclic_rate_orig, each = length(knot_loc$year))
	b_init_theta_orig <-  rep(b_init_cyclic_theta_orig, each = length(knot_loc$year))

	### Convert to reparameterized basis 
	b_shape_init <- c(t( ginv(basis$z) %*% b_init_shape_orig))
	b_rate_init <- c(t(ginv(basis$z) %*% b_init_rate_orig))
	b_theta_init <- c(t(ginv(basis$z) %*% b_init_theta_orig))

	### Estimate rho init and priors
	rho_init_shape <- c(cyclic_output$rho_init$shape, 10)
	rho_init_rate <- c(cyclic_output$rho_init$rate, 10)
	rho_init_theta <- c(cyclic_output$rho_init$theta, 10)

	rho_shape_prior <- c(rho_init_shape[1]-3, rho_init_shape[1] + 3, 0,  12)
	rho_rate_prior <- c(rho_init_rate[1]-3, rho_init_rate[1] + 3, 0,  12)
	rho_theta_prior <- c(rho_init_theta[1]-3, rho_init_theta[1] + 3, 0,  12)

	### Convert to lambda
	lambda_shape_init <- exp(rho_init_shape)
	lambda_rate_init <- exp(rho_init_rate)
	lambda_theta_init <- exp(rho_init_theta)

	### Create output list and return
	output_list <- list(
		b_0 = cyclic_output$b_0, 
		b_init = list(shape = b_shape_init, rate = b_rate_init, theta = b_theta_init), 
		lambda_init = list(shape = lambda_shape_init, rate = lambda_rate_init, theta = lambda_theta_init), 
		rho_init = list(shape = rho_init_shape, rate = rho_init_rate, theta = rho_init_theta),
		rho_prior = list(shape = rho_shape_prior, rate = rho_rate_prior, theta = rho_theta_prior)
	)

	return(output_list)
}





#' Estimate priors for tensor spline
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @return A matrix of the infile
#' @import dplyr
#' @export
prior_tensor <- function(data, knot_loc ){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)
	
	### Add a year column if it doesn't exist
	### It should, but just in case
	if(!("year" %in% colnames(data))) {
		data <- data %>%
			mutate(year = lubridate::year(date))		
	}

	### Remove any leap days for fitting purposes
	data <- data %>%
		filter(jdate <=365)

	### Add in a draw column with the unique years so this can later be joined
	### Base number of samples from MLE off of this
	data <- data %>% 
		mutate(draw = as.numeric(factor(year)))
	n_years <- max(data$draw)

	### Create MLE estimate. This will be used to set priors and initial values
	### Create preliminary periods based on knot locations
#	year_cuts <- unique(c(min(data$year, na.rm=TRUE),knot_loc$year, max(data$year, na.rm=TRUE)))
#	year_cuts <- c(0, year_cuts)

	### Check how many observations 
#	countbycut <- data %>% 
#		drop_na(precip) %>%
#		mutate(period = cut(year, year_cuts)) %>%
#		group_by(period, .drop = FALSE) %>% 
#		count() %>%
#		ungroup() 

	### Merge periods that have zero observations or less than 4 years
#	nonzero_test <- c(TRUE, countbycut$n != 0)
#	multiyear_test <- c(TRUE, countbycut$n > (365*4))
#	year_cuts <- year_cuts[nonzero_test & multiyear_test]

	### Cut data to periods and generate MLE estiamtes for the gamma distribution
#	data <- data %>%
#		mutate(period = cut(year, year_cuts)) 

	### Create MLE estimate
	### Fit each day with MLE
	mle_fit <- mle_gamma_fit(jdate = data$jdate, values = data$precip, n = n_years)
	

#	mle_fit <- data %>%
###		dplyr::group_by(jdate, period) %>%
#		dplyr::group_by(jdate) %>%
#		dplyr::filter(precip > 0) %>%
#		dplyr::summarise(shape = fitdistrplus::fitdist(precip, "gamma")[[1]][[1]], rate = fitdistrplus::fitdist(precip, "gamma")[[1]][[2]]) %>%
#		dplyr::mutate(scale = 1/rate) %>%
#		dplyr::mutate(mean = scale * shape) %>%
#		dplyr::mutate(disp = 1/shape) %>%
#		dplyr::ungroup()

	### Do some reorganization to join and then remove 30 year periods
#	mle_fit <- data %>% 
#		full_join(mle_fit) %>%
#		select(-period)

#	data <- data %>%
#		select(-period)

	data <- data %>%
		select(date, jdate, year) %>%
		inner_join(mle_fit$draws, by = c("jdate"))
	
	### Fit a tensor product spline model using mgcv
	shape_model <- mgcv::gam(log(shape) ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=data, knots = knot_loc, select=FALSE)
	rate_model <- mgcv::gam(log(rate) ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=data, knots = knot_loc, select=FALSE)
	theta_model <- mgcv::gam(zero ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=data, knots = knot_loc, select=FALSE,family=binomial)


	### Create the prior for the mean intercept
	b_0_shape_prior <- c(summary(shape_model)$p.table[1], summary(shape_model)$p.table[2])
	b_0_rate_prior <- c(summary(rate_model)$p.table[1], summary(rate_model)$p.table[2])
	b_0_theta_prior <- c(summary(theta_model)$p.table[1], summary(theta_model)$p.table[2])

	### Create a vector for intializing the mean
	b_shape_init <- c(coef(shape_model)[2:length(coef(shape_model))])
	b_rate_init <- c(coef(rate_model)[2:length(coef(rate_model))])
	b_theta_init <- c(coef(theta_model)[2:length(coef(theta_model))])

	### Extract smoothing penalty
	lambda_shape_init <- unname(c(shape_model$sp ))
	lambda_rate_init <- unname(c(rate_model$sp))
	lambda_theta_init <- unname(c(theta_model$sp))

	### Calculate alpha the rescaling factor and then rescale
	shape_alpha <- sapply(shape_model$smooth, "[[", "S.scale") / lambda_shape_init
	lambda_shape_init <- c(lambda_shape_init / shape_alpha)

	rate_alpha <- sapply(rate_model$smooth, "[[", "S.scale") / lambda_rate_init
	lambda_rate_init <- c(lambda_rate_init / rate_alpha)

	theta_alpha <- sapply(theta_model$smooth, "[[", "S.scale") / lambda_theta_init
	lambda_theta_init <- c(lambda_theta_init / theta_alpha)

	### Extract penalty matrix
	s_matrix <- sapply(shape_model$smooth, "[[", "S")

	s_mat_1 <- s_matrix[[1]]
	s_mat_2 <- s_matrix[[2]]
#	s_inv_1 <- ginv(s_mat_1)
#	s_inv_2 <- ginv(s_mat_2)

	#mean_samples <- mvrnorm(n=100, b_mean_init, ginv(s_mat_1*lambda_mean_init[1] + s_mat_2*lambda_mean_init[2]))
	#scale_samples <- mvrnorm(n=100, b_mean_init, ginv(s_mat_1*lambda_mean_init[1] + s_mat_2*lambda_mean_init[2]))

	#plot(b_mean_init, type="l", ylim=c(-2,2))
	#lines(mean_samples[1,], col="red")

	#plot(mean_samples[1,1:6], type="l", ylim=c(-1,1))
	#lines(mean_samples[2,1:6], col="red")
	#lines(mean_samples[3,1:6], col="red")
	#lines(mean_samples[4,1:6], col="red")
	#lines(mean_samples[5,1:6], col="red")
	#lines(mean_samples[1,7:13], col="blue")
	#lines(mean_samples[2,7:13], col="blue")
	#lines(mean_samples[3,7:13], col="blue")
	#lines(mean_samples[4,7:13], col="blue")
	#lines(mean_samples[5,7:13], col="blue")

	rho_init_shape <- log(lambda_shape_init)
	rho_init_rate <- log(lambda_rate_init)
	rho_init_theta <- log(lambda_theta_init)

	rho_shape_prior <- c(rho_init_shape-3, rho_init_shape + 3)
	rho_rate_prior <- c(rho_init_rate-3, rho_init_rate + 3)
	rho_theta_prior <- c(rho_init_theta-3, rho_init_theta + 3)

	### Create output list and return
	output_list <- list(
		b_0 = list(shape = b_0_shape_prior, rate = b_0_rate_prior,  theta = b_0_theta_prior), 
		b_init = list(shape = b_shape_init, rate = b_rate_init, theta = b_theta_init), 
		lambda_init = list(shape = lambda_shape_init, rate = lambda_rate_init, theta = lambda_theta_init), 
		rho_init = list(shape = rho_init_shape, rate = rho_init_rate, theta = rho_init_theta),
		rho_prior = list(shape = rho_shape_prior, rate = rho_rate_prior, theta = rho_theta_prior)
	)

	return(output_list)

}






#' Magic Lambda
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
magic_lambda <- function(pre_model, lambda_year = -1, lambda_ratio = 200){

	if(lambda_year == -1 & lambda_ratio == -1){
		lambda_year <- 5000
	}

	### Extract the penalty matrix for the first parameter
	s_1 <- pre_model$s_reparam[[1]]
	s_2 <- pre_model$s_reparam[[2]]

	match_extremes <- function(lambda_jdate, lambda_year, target){
		penalty <- lambda_jdate * s_1 + (lambda_year) * s_2
		random_draw <- mvrnorm(n=1000, rep(0,length(mean_beta)), ginv(penalty))
		row_max <- apply(random_draw, 1, function(x){max(abs(x))})
		return(abs(target - median(row_max)))
	}

	match_extremes_ratio <- function(lambda_jdate, lambda_ratio, target){
		penalty <- lambda_jdate * s_1 + (lambda_jdate * lambda_ratio) * s_2
		random_draw <- mvrnorm(n=1000, rep(0,length(mean_beta)), ginv(penalty))
		row_max <- apply(random_draw, 1, function(x){max(abs(x))})
		return(abs(target - median(row_max)))
	}

	mean_beta <- pre_model$b_init$mean
	scale_beta <- pre_model$b_init$scale

	if(lambda_year == -1 & lambda_ratio > 0){
	lambda_lower <- optimize(match_extremes_ratio, interval = c(1E-2,1E6), lambda_ratio = lambda_ratio, target = 0.8*max(abs(mean_beta)))
	lambda_lower <- lambda_lower$minimum

	lambda_upper <- optimize(match_extremes_ratio, interval = c(1E-2,1E6), lambda_ratio = lambda_ratio, target = 0.3*max(abs(mean_beta)))
	lambda_upper <- lambda_upper$minimum

	lambda_mean_prior <- c(lambda_lower, lambda_upper, lambda_lower* lambda_ratio, lambda_upper * lambda_ratio)

	lambda_lower <- optimize(match_extremes_ratio, interval = c(1E-2,1E6), lambda_ratio = lambda_ratio, target = 0.8*max(abs(scale_beta)))
	lambda_lower <- lambda_lower$minimum

	lambda_upper <- optimize(match_extremes_ratio, interval = c(1E-2,1E6), lambda_ratio = lambda_ratio, target = 0.3*max(abs(scale_beta)))
	lambda_upper <- lambda_upper$minimum

	lambda_scale_prior <- c(lambda_lower, lambda_upper, lambda_lower* lambda_ratio, lambda_upper * lambda_ratio)
	} else {
	lambda_lower <- optimize(match_extremes, interval = c(1E-2,1E6), lambda_year = lambda_year, target = 0.8*max(abs(mean_beta)))
	lambda_lower <- lambda_lower$minimum

	lambda_upper <- optimize(match_extremes, interval = c(1E-2,1E6), lambda_year = lambda_year, target = 0.3*max(abs(mean_beta)))
	lambda_upper <- lambda_upper$minimum

	lambda_mean_prior <- c(lambda_lower, lambda_upper, lambda_year)

	lambda_lower <- optimize(match_extremes, interval = c(1E-2,1E6), lambda_year = lambda_year, target = 0.8*max(abs(scale_beta)))
	lambda_lower <- lambda_lower$minimum

	lambda_upper <- optimize(match_extremes, interval = c(1E-2,1E6), lambda_year = lambda_year, target = 0.3*max(abs(scale_beta)))
	lambda_upper <- lambda_upper$minimum

	lambda_scale_prior <- c(lambda_lower, lambda_upper, lambda_year)
	}


	return(list(mean = lambda_mean_prior, scale = lambda_scale_prior))
}



