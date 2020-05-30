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
		preproc_output <- prior_tensor(data = basis$data, knot_loc = basis$knot_loc)
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
	mean_model <- mgcv::gam(mean ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=mle_fit$draws, knots = knot_loc, select=FALSE)
	scale_model <- mgcv::gam(scale ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=mle_fit$draws, knots = knot_loc, select=FALSE)
	
	### Create the prior for the mean intercept
	b_0_mean_prior <- c(summary(mean_model)$p.table[1], summary(mean_model)$p.table[2])
	b_0_scale_prior <- c(summary(scale_model)$p.table[1], summary(scale_model)$p.table[2])

	### Create a vector for intializing the mean
	b_mean_init <- c(coef(mean_model)[2:length(coef(mean_model))])
	b_scale_init <- c(coef(scale_model)[2:length(coef(scale_model))])

	### Extract smoothing penalty
	lambda_mean_init <- unname(c(mean_model$sp ))
	lambda_scale_init <- unname(c(scale_model$sp))

	### Calculate alpha the rescaling factor and then rescale
	mean_alpha <- sapply(mean_model$smooth, "[[", "S.scale") / lambda_mean_init
	lambda_mean_init <- lambda_mean_init / mean_alpha

	scale_alpha <- sapply(scale_model$smooth, "[[", "S.scale") / lambda_scale_init
	lambda_scale_init <- lambda_scale_init / scale_alpha

	### Extract penalty matrix
	s_matrix <- mean_model$smooth[[1]]$S[[1]]
	s_inv <- MASS::ginv(s_matrix)

	### Rescale to match standard devs
	### Might not need to do this
	mean_samples <- mvrnorm(n=100, b_mean_init, s_inv)
	scale_samples <- mvrnorm(n=100, b_scale_init, s_inv)

	lambda_mean <- 1/((sd(b_mean_init)/mean(apply(mean_samples, 1, sd)))^2)
	lambda_scale <- 1/((sd(b_scale_init)/mean(apply(scale_samples, 1, sd)))^2)

	#mean_testing <- mvrnorm(n=100, b_mean_init, s_inv/lambda_mean)
	#scale_testing <- mvrnorm(n=100, b_scale_init, s_inv/lambda_scale)	
	#mean_testing <- mvrnorm(n=100, b_mean_init, s_inv/lambda_mean_init)
	
	### Create output list and return
	output_list <- list(b_0 = list(mean = b_0_mean_prior, scale = b_0_scale_prior), b_init = list(mean = b_mean_init, scale = b_scale_init), lambda_init = list(mean = lambda_mean, scale = lambda_scale))
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
		inner_join(mle_fit$draws, by = c("jdate", "draw"))
	
	### Fit a tensor product spline model using mgcv
	mean_model <- mgcv::gam(mean ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=data, knots = knot_loc, select=FALSE)
	scale_model <- mgcv::gam(scale ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=data, knots = knot_loc, select=FALSE)

	### Create the prior for the mean intercept
	b_0_mean_prior <- c(summary(mean_model)$p.table[1], summary(mean_model)$p.table[2])
	b_0_scale_prior <- c(summary(scale_model)$p.table[1], summary(scale_model)$p.table[2])

	### Create a vector for intializing the mean
	b_mean_init <- c(coef(mean_model)[2:length(coef(mean_model))])
	b_scale_init <- c(coef(scale_model)[2:length(coef(scale_model))])

	### Extract smoothing penalty
	lambda_mean_init <- unname(c(mean_model$sp ))
	lambda_scale_init <- unname(c(scale_model$sp))

	### Calculate alpha the rescaling factor and then rescale
	mean_alpha <- sapply(mean_model$smooth, "[[", "S.scale") / lambda_mean_init
	lambda_mean_init <- lambda_mean_init / mean_alpha

	scale_alpha <- sapply(scale_model$smooth, "[[", "S.scale") / lambda_scale_init
	lambda_scale_init <- lambda_scale_init / scale_alpha

	### Extract penalty matrix
	s_matrix <- sapply(mean_model$smooth, "[[", "S")

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

	### Catch for excessively high values
	lambda_scale_init[lambda_scale_init > 1E11] <- 1E11
	lambda_scale_init[lambda_scale_init > 1E11] <- 1E11

	### Create output list and return
	output_list <- list(b_0 = list(mean = b_0_mean_prior, scale = b_0_scale_prior), b_init = list(mean = b_mean_init, scale = b_scale_init), lambda_init = list(mean = lambda_mean_init, scale = lambda_scale_init))
	return(output_list)

}

