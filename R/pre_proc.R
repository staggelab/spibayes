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
pre_proc <- function(data, type, knot_loc, lambda_shape = c(500, 5000, 0.5, 5000, 5000), year_pen_ratio = 100){

	### Catch - if there is no precip column
	if(!("precip" %in% colnames(data))){
 		stop('There must be a column titled precip (lowercase)')
	}

	### Create a column for zero precip
	data <- data %>%
		mutate(zero = precip == 0)

	### Remove any leap days for fitting purposes
	data <- data %>%
		filter(jdate <=365)

	if (type == "cyclic"){
		preproc_output <- prior_cyclic(data = data, knot_loc = knot_loc, lambda_shape = lambda_shape[1])
	} else if (type == "tensor"){
		preproc_output <- prior_tensor(data = data, knot_loc = knot_loc, lambda_shape = lambda_shape, year_pen_ratio = year_pen_ratio)
	}

	### Create the output and return
	preproc_output$data_all <- data
	preproc_output$data_pos <- data %>% filter(precip > 0)
	preproc_output$type <- type
	preproc_output$knot_loc <- knot_loc

	return(list(input=preproc_output))
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
prior_cyclic <- function(data, knot_loc, lambda_shape ){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)
	
	### Extract positive precipitation and assign all data
	data_pos <- data %>%
		filter(precip > 0)
	data_all <- data

	### How many years of data to sample from MLE estimate
	if("year" %in% colnames(data_all)) {
		n_years <- length(unique(data_all$year))
	} else {
		n_years <- length(unique(lubridate::year(data_all$date)))
	}

	### Fit a cyclic spline model for gamma distributed precip
	gamma_model <- gam(list(precip 
		~ ti(jdate, bs=c("cc"), k = n_knots[1]),
		~ ti(jdate, bs=c("cc"), k = n_knots[1])),
		data = data_pos, 
		family=gammals(link=list("identity","log")),
		fit = TRUE, 
		select = TRUE)

	### Fit a cyclic spline model using mgcv for zeros
	theta_model <- mgcv::gam(zero ~ ti(jdate, bs=c("cc"), k = c(n_knots[1])), 
		data=data_all, 
		family=binomial,
		fit = TRUE,
		select=TRUE)

	### Extract basis matrix for mean and dispersion
	X_mean_jdate <- PredictMat(gamma_model$smooth[[1]], data_pos)
	X_disp_jdate <- PredictMat(gamma_model$smooth[[2]], data_pos)

	### Extract basis matrix For zeros
	X_theta_jdate <- PredictMat(gamma_model$smooth[[1]], data)

	### Extract S Penalty matrices for mean and dispersion
	s_mean_jdate <- gamma_model$smooth[[1]]$S[[1]]
	s_disp_jdate <- gamma_model$smooth[[2]]$S[[1]]

	### Extract S Penalty matrices for zero precip
	s_theta_jdate <- theta_model$smooth[[1]]$S[[1]]

	### Extract S scaling function
	s_scale_pos <- sapply(gamma_model$smooth, function(x){x$S.scale})
	s_scale_zero <- sapply(gamma_model$smooth, function(x){x$S.scale})
	
	### Extract basis dimensions
	basis_dim_mean <- dim(X_mean_jdate)[2]
	basis_dim_disp <- dim(X_disp_jdate)[2]
	basis_dim_theta <- dim(X_theta_jdate)[2]
	
	### Extract model coefficients
	coef_init_df <- data.frame(t(coef(gamma_model)))

	### Process the beta intercept
	b_0_mean_init <- unlist(select(coef_init_df,contains("Intercept"))[1])
	b_0_disp_init <- unlist(select(coef_init_df,contains("Intercept"))[2])
	b_0_theta_init <- unlist(coef(theta_model)[1])

	### Calculate b_0 prior matrix
	b_0_prior <- as.matrix(summary(gamma_model)$p.table[,1:2])
	b_0_prior <- rbind(b_0_prior, matrix(summary(theta_model)$p.table[1:2], 1, 2))
	colnames(b_0_prior) <- c("est", "std_dev")
	rownames(b_0_prior) <- c("mean", "disp", "theta")
	### Expand the range by an order of magnitude
	b_0_prior[,2] <- b_0_prior[,2] * 10

	### Process the beta init
	b_init_mean_jdate <- c(unlist(select(coef_init_df,contains("ti.jdate"))))
	b_init_disp_jdate <- c(unlist(select(coef_init_df,contains("ti.1.jdate"))))
	b_init_theta_jdate <- coef(theta_model)[2:length(coef(theta_model))]

	### Pull the estimated penalty from a jdate only model
	lambda_mean_jdate_est <- gamma_model$sp[1]
	lambda_disp_jdate_est <- gamma_model$sp[2]
	lambda_theta_jdate_est <- theta_model$sp[1]

	### Create matrices for gamma distributed priors on lambda
	lambda_mean_prior <- matrix(NA, 1, 2)
	### Add the shape parameter for lambda. Tight near normal for everything except the double penalty
	lambda_mean_prior[,1] <- lambda_shape

	### Replicate for other lambdas
	lambda_disp_prior <-lambda_mean_prior
	lambda_theta_prior <-lambda_mean_prior

	### Calculate the rate parameter
	lambda_mean_prior[,2] <-lambda_mean_prior[,1]/lambda_mean_jdate_est
	lambda_disp_prior[,2] <-lambda_disp_prior[,1]/lambda_disp_jdate_est
	lambda_theta_prior[,2] <-lambda_theta_prior[,1]/lambda_theta_jdate_est

	### Sigma
	sigma_df <- data_pos %>%
		group_by(jdate) %>%
		group_modify(~ est_gamma_sigma(.x$precip, n = 100)) %>%
		ungroup()

	sigma_theta <- data_all %>%
		group_by(jdate) %>%
		group_modify(~ est_theta_sigma(.x$zero)) %>%
		ungroup()

	sigma_mean <- median(sigma_df$sigma_mean)
	sigma_disp <- median(sigma_df$sigma_disp)
	sigma_theta <- median(sigma_theta$sigma_theta)

	sigma_mean_prior <- matrix(NA, 1, 2)
	sigma_mean_prior[,1] <- 99
	sigma_disp_prior <- sigma_mean_prior
	sigma_theta_prior <- sigma_mean_prior
	sigma_mean_prior[,2] <-sigma_mean_prior[,1]/sigma_mean
	sigma_disp_prior[,2] <-sigma_disp_prior[,1]/sigma_disp
	sigma_theta_prior[,2] <-sigma_theta_prior[,1]/sigma_theta

	### Create output list and return
	init_vals_output <- list(
		x_matrix = list(mean = list(jdate = X_mean_jdate), 
						disp = list(jdate = X_disp_jdate), 
						theta = list(jdate = X_theta_jdate)
			),
		basis_dim = list(mean = basis_dim_mean, 
						disp = basis_dim_disp, 
						theta = basis_dim_theta
			),
		s_matrix = list(mean = list(jdate = s_mean_jdate), 
						disp = list(jdate = s_disp_jdate), 
						theta = list(jdate = s_theta_jdate)
			),
		b_0_init = list(mean = b_0_mean_init, 
						disp = b_0_disp_init, 
						theta = b_0_theta_init
			),
		b_0_prior = b_0_prior,
		b_init = list(mean = list(jdate = b_init_mean_jdate), 
						disp = list(jdate = b_init_disp_jdate), 
						theta = list(jdate = b_init_theta_jdate)
			),
		lambda_init = list(mean = lambda_mean_jdate_est, 
						disp = lambda_disp_jdate_est, 
						theta = lambda_theta_jdate_est
			),
		lambda_prior = list(mean = lambda_mean_prior, 
						disp = lambda_disp_prior, 
						theta = lambda_theta_prior
			),
		sigma_init = list(mean = sigma_mean, 
						disp = sigma_disp, 
						theta = sigma_theta
			),
		sigma_prior = list(mean = sigma_mean_prior, 
						disp = sigma_disp_prior, 
						theta = sigma_theta_prior
			),
		model = list(gamma = gamma_model, 
						theta = theta_model
			)
	)

	return(init_vals_output)
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
prior_tensor <- function(data, knot_loc, lambda_shape, year_pen_ratio){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)
	
	### Extract positive precipitation and assign all data
	data_pos <- data %>%
		filter(precip > 0)
	data_all <- data

	### How many years of data to sample from MLE estimate
	if("year" %in% colnames(data_all)) {
		n_years <- length(unique(data_all$year))
	} else {
		n_years <- length(unique(lubridate::year(data_all$date)))
	}

	### Create cyclic spline basis for precipitation
	gamma_model <- mgcv::gam(list(precip 
		~ ti(jdate, bs=c("cc"), k = n_knots[1]) + ti(year, bs=c("cr"), k = n_knots[2]) + ti(jdate,year, bs=c("cc", "cr"), k = n_knots),
		~ ti(jdate, bs=c("cc"), k = n_knots[1]) + ti(year, bs=c("cr"), k = n_knots[2]) + ti(jdate,year, bs=c("cc", "cr"), k = n_knots)),
		data = data_pos, 
		knots =  knot_loc, 
		family=gammals(link=list("identity","log")),
		fit = FALSE, 
		select = TRUE)

	### Create cyclic spline basis for zero precipitation
	theta_model <- mgcv::gam(zero 
		~ ti(jdate, bs=c("cc"), k = n_knots[1]) + ti(year, bs=c("cr"), k = n_knots[2]) + ti(jdate,year, bs=c("cc", "cr"), k = n_knots), 
		data = data, 
		knots =  knot_loc, 
		family=binomial,
		fit = FALSE, 
		select = TRUE)

	### Extract basis matrix for mean
	X_mean_jdate <- PredictMat(gamma_model$smooth[[1]], data_pos)
	X_mean_year <- PredictMat(gamma_model$smooth[[2]], data_pos)
	X_mean_tensor <- PredictMat(gamma_model$smooth[[3]], data_pos)

	### Extract basis matrix for dispersion
	X_disp_jdate <- PredictMat(gamma_model$smooth[[4]], data_pos)
	X_disp_year <- PredictMat(gamma_model$smooth[[5]], data_pos)
	X_disp_tensor <- PredictMat(gamma_model$smooth[[6]], data_pos)

	### Extract basis matrix For zeros
	X_theta_jdate <- PredictMat(theta_model$smooth[[1]], data)
	X_theta_year <- PredictMat(theta_model$smooth[[2]], data)
	X_theta_tensor <- PredictMat(theta_model$smooth[[3]], data)

	### Extract S Penalty matrices for mean and dispersion
	s_mean_jdate <- gamma_model$smooth[[1]]$S[[1]]
	s_mean_year <- gamma_model$smooth[[2]]$S[[1]]
	s_mean_year_double <- gamma_model$smooth[[2]]$S[[2]]
	s_mean_tensor_jdate <- gamma_model$smooth[[3]]$S[[1]]
	s_mean_tensor_year <- gamma_model$smooth[[3]]$S[[2]]

	s_disp_jdate <- gamma_model$smooth[[4]]$S[[1]]
	s_disp_year <- gamma_model$smooth[[5]]$S[[1]]
	s_disp_year_double <- gamma_model$smooth[[5]]$S[[2]]
	s_disp_tensor_jdate <- gamma_model$smooth[[6]]$S[[1]]
	s_disp_tensor_year <- gamma_model$smooth[[6]]$S[[2]]

	### Extract S Penalty matrices for zero precip
	s_theta_jdate <- theta_model$smooth[[1]]$S[[1]]
	s_theta_year <- theta_model$smooth[[2]]$S[[1]]
	s_theta_year_double <- theta_model$smooth[[2]]$S[[2]]
	s_theta_tensor_jdate <- theta_model$smooth[[3]]$S[[1]]
	s_theta_tensor_year <- theta_model$smooth[[3]]$S[[2]]

	### Extract S scaling function
	s_scale_pos <- sapply(gamma_model$smooth, function(x){x$S.scale})
	s_scale_zero <- sapply(theta_model$smooth, function(x){x$S.scale})
	
	### Extract basis dimensions
	basis_dim_mean <- c(dim(X_mean_jdate)[2], dim(X_mean_year)[2], dim(X_mean_tensor)[2])
	basis_dim_disp <- c(dim(X_disp_jdate)[2], dim(X_disp_year)[2], dim(X_disp_tensor)[2])
	basis_dim_theta <- c(dim(X_theta_jdate)[2], dim(X_theta_year)[2], dim(X_theta_tensor)[2])

	### Evaluate cyclic spline to get init values
	cyclic_init <- prior_cyclic(data = data, knot_loc = list(jdate = knot_loc$jdate), lambda_shape = lambda_shape[1])

	### Add on zeros for everything that is not the cyclic spline
	b_init_mean_jdate <- cyclic_init$b_init$mean$jdate
	b_init_mean_year <- rep(0,basis_dim_mean[2])
	b_init_mean_tensor <- rep(0,basis_dim_mean[3])

	b_init_disp_jdate <- cyclic_init$b_init$disp$jdate
	b_init_disp_year <- rep(0,basis_dim_disp[2])
	b_init_disp_tensor <- rep(0,basis_dim_disp[3])

	b_init_theta_jdate <- cyclic_init$b_init$theta$jdate
	b_init_theta_year <- rep(0,basis_dim_theta[2])
	b_init_theta_tensor <- rep(0,basis_dim_theta[3])

	### Extract the intercept from the cyclic run
	b_0_mean_init <- cyclic_init$b_0_init$mean
	b_0_disp_init <- cyclic_init$b_0_init$disp
	b_0_theta_init <- cyclic_init$b_0_init$theta

	### Extract jdate lambda from cyclic
	lambda_mean_jdate_est <- cyclic_init$lambda_init$mean
	lambda_disp_jdate_est <- cyclic_init$lambda_init$disp
	lambda_theta_jdate_est <- cyclic_init$lambda_init$theta

	### Apply penalty ratio to each lambda
	lambda_mean_init <- c(lambda_mean_jdate_est, lambda_mean_jdate_est*year_pen_ratio, lambda_mean_jdate_est*year_pen_ratio^2, lambda_mean_jdate_est, lambda_mean_jdate_est*year_pen_ratio)
	lambda_disp_init <- c(lambda_disp_jdate_est, lambda_disp_jdate_est*year_pen_ratio, lambda_disp_jdate_est*year_pen_ratio^2, lambda_disp_jdate_est, lambda_disp_jdate_est*year_pen_ratio)
	lambda_theta_init <- c(lambda_theta_jdate_est, lambda_theta_jdate_est*year_pen_ratio, lambda_theta_jdate_est*year_pen_ratio^2, lambda_theta_jdate_est, lambda_theta_jdate_est*year_pen_ratio)

	### Add correct names
	names(lambda_mean_init) <- c("jdate", "year", "year_double", "tensor_jdate", "tensor_year")
	names(lambda_disp_init) <- c("jdate", "year", "year_double", "tensor_jdate", "tensor_year")
	names(lambda_theta_init) <- c("jdate", "year", "year_double", "tensor_jdate", "tensor_year")

	### Create matrices for gamma distributed priors on lambda
	lambda_mean_prior <- matrix(NA, 5, 2)
	### Add the shape parameter for lambda. Tight near normal for everything except the double penalty
	lambda_mean_prior[,1] <- lambda_shape

	### Replicate for other lambdas
	lambda_disp_prior <-lambda_mean_prior
	lambda_theta_prior <-lambda_mean_prior

	### Calculate the rate parameter
	lambda_mean_prior[,2] <-lambda_mean_prior[,1]/lambda_mean_init
	lambda_disp_prior[,2] <-lambda_disp_prior[,1]/lambda_disp_init
	lambda_theta_prior[,2] <-lambda_theta_prior[,1]/lambda_theta_init


	### Sigma
	sigma_df <- data_pos %>%
		group_by(jdate) %>%
		group_modify(~ est_gamma_sigma(.x$precip, n = 100)) %>%
		ungroup()

	sigma_theta <- data_all %>%
		group_by(jdate) %>%
		group_modify(~ est_theta_sigma(.x$zero)) %>%
		ungroup()

	sigma_mean <- median(sigma_df$sigma_mean)
	sigma_disp <- median(sigma_df$sigma_disp)
	sigma_theta <- median(sigma_theta$sigma_theta)

	sigma_mean_prior <- matrix(NA, 1, 2)
	sigma_mean_prior[,1] <- 99
	sigma_disp_prior <- sigma_mean_prior
	sigma_theta_prior <- sigma_mean_prior
	sigma_mean_prior[,2] <-sigma_mean_prior[,1]/sigma_mean
	sigma_disp_prior[,2] <-sigma_disp_prior[,1]/sigma_disp
	sigma_theta_prior[,2] <-sigma_theta_prior[,1]/sigma_theta

	### Create output list and return
	init_vals_output <- list(
		x_matrix = list(mean = list(jdate = X_mean_jdate, year = X_mean_year, tensor = X_mean_tensor), 
						disp = list(jdate = X_disp_jdate, year = X_disp_year, tensor = X_disp_tensor), 
						theta = list(jdate = X_theta_jdate, year = X_theta_year, tensor = X_theta_tensor)
			),
		basis_dim = list(mean = basis_dim_mean, 
						disp = basis_dim_disp, 
						theta = basis_dim_theta
			),
		s_matrix = list(mean = list(jdate = s_mean_jdate, year = s_mean_year, year_double = s_mean_year_double, tensor_jdate = s_mean_tensor_jdate, tensor_year = s_mean_tensor_year), 
						disp = list(jdate = s_disp_jdate, year = s_disp_year, year_double = s_disp_year_double, tensor_jdate = s_disp_tensor_jdate, tensor_year = s_disp_tensor_year), 
						theta = list(jdate = s_theta_jdate, year = s_theta_year, year_double = s_theta_year_double, tensor_jdate = s_theta_tensor_jdate, tensor_year = s_theta_tensor_year)
			),
		b_0_init = list(mean = b_0_mean_init, 
						disp = b_0_disp_init, 
						theta = b_0_theta_init
			),
		b_0_prior = cyclic_init$b_0_prior,
		b_init = list(mean = list(jdate = b_init_mean_jdate, year = b_init_mean_year, tensor = b_init_mean_tensor), 
						disp = list(jdate = b_init_disp_jdate, year = b_init_disp_year, tensor = b_init_disp_tensor), 
						theta = list(jdate = b_init_theta_jdate, year = b_init_theta_year, tensor = b_init_theta_tensor)
			),
		lambda_init = list(mean = lambda_mean_init, 
						disp = lambda_disp_init, 
						theta = lambda_theta_init
			),
		lambda_prior = list(mean = lambda_mean_prior, 
						disp = lambda_disp_prior, 
						theta = lambda_theta_prior
			),
		sigma_init = list(mean = sigma_mean, 
						disp = sigma_disp, 
						theta = sigma_theta
			),
		sigma_prior = list(mean = sigma_mean_prior, 
						disp = sigma_disp_prior, 
						theta = sigma_theta_prior
			),
		model = list(gamma = gamma_model, 
						theta = theta_model
			)
	)

	return(init_vals_output)
}





#' Extract draws from a fitdistrplus distribution using standard error variance covariance matrix
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param n Number of draws to calculate
#' @export
est_gamma_sigma <- function(data, n){
	data <- data[data > 0 & !is.na(data)]
	fitted <- fitdist(data, "gamma")
	draws <- mvrnorm(n, mu = coef(fitted), Sigma = vcov(fitted)) 
	draws <- draws %>%
		as.data.frame() %>%
		mutate(scale = 1/rate) %>%
		mutate(mean = scale * shape) %>%
		mutate(disp = 1/shape) %>%
		mutate(draw = seq(1,n)) %>%
		mutate(type = "draw")

	estimate <- data.frame(shape = coef(fitted)[[1]], rate = coef(fitted)[[2]]) %>%
		mutate(scale = 1/rate) %>%
		mutate(mean = scale * shape) %>%
		mutate(disp = 1/shape) %>%
		mutate(draw = 0) %>%
		mutate(type = "estimate")

	output_df <- data.frame(sigma_mean = sd(log(draws$mean) - log(estimate$mean)), sigma_disp = sd(log(draws$disp) - log(estimate$disp)))

	return(output_df)
}


est_theta_sigma <- function(data){
	glm_pos <- glm(data ~ 1, family = "binomial")
	output_df <- data.frame(sigma_theta = summary(glm_pos)$coefficients[[2]])

	return(output_df)
}


