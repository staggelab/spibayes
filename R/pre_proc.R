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
	
	### Create MLE estimate
	### Fit each day with MLE
	mle_fit <- data %>%
		dplyr::group_by(jdate) %>%
		dplyr::filter(precip > 0) %>%
		dplyr::summarise(shape = fitdistrplus::fitdist(precip, "gamma")[[1]][[1]], rate = fitdistrplus::fitdist(precip, "gamma")[[1]][[2]]) %>%
		dplyr::mutate(scale = 1/rate) %>%
		dplyr::mutate(mean = scale * shape) %>%
		dplyr::mutate(disp = 1/shape) %>%
		dplyr::ungroup()

	### Fit a cyclic spline model using mgcv
	mean_model <- mgcv::gam(mean ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=mle_fit, knots = knot_loc, select=TRUE)
	disp_model <- mgcv::gam(disp ~ s(jdate, bs=c("cc"), k = c(n_knots)), data=mle_fit, knots = knot_loc, select=TRUE)
	
	### Create the prior for the mean intercept
	b_0_mean_prior <- c(summary(mean_model)$p.table[1], summary(mean_model)$p.table[2])
	b_0_disp_prior <- c(summary(disp_model)$p.table[1], summary(disp_model)$p.table[2])

	### Create a vector for intializing the mean
	b_mean_init <- c(coef(mean_model)[2:length(coef(mean_model))])
	b_disp_init <- c(coef(disp_model)[2:length(coef(disp_model))])

	lambda_mean_init <- unname(c(mean_model$sp ))
	lambda_disp_init <- unname(c(disp_model$sp))

	### Create output list and return
	output_list <- list(b_0 = list(mean = b_0_mean_prior, disp = b_0_disp_prior), b_init = list(mean = b_mean_init, disp = b_disp_init), lambda_init = list(mean = lambda_mean_init, disp = lambda_disp_init))
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
	
	### Create MLE estimate. This will be used to set priors and initial values
	### Do this in 30 year subsets
	data <- data %>%
		mutate(period = cut(year, seq(0,2200, 30))) 

	mle_fit <- data %>%
		dplyr::group_by(jdate, period) %>%
		dplyr::filter(precip > 0) %>%
		dplyr::summarise(shape = fitdistrplus::fitdist(precip, "gamma")[[1]][[1]], rate = fitdistrplus::fitdist(precip, "gamma")[[1]][[2]]) %>%
		dplyr::mutate(scale = 1/rate) %>%
		dplyr::mutate(mean = scale * shape) %>%
		dplyr::mutate(disp = 1/shape) %>%
		dplyr::ungroup()

	### Do some reorganization to join and then remove 30 year periods
	mle_fit <- data %>% 
		full_join(mle_fit) %>%
		select(-period)

	data <- data %>%
		select(-period)

	### Fit a tensor product spline model using mgcv
	mean_model <- mgcv::gam(mean ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=mle_fit, knots = knot_loc, select=TRUE)
	disp_model <- mgcv::gam(disp ~ te(jdate,year, bs=c("cc", "cr"), k = c(n_knots)), data=mle_fit, knots = knot_loc, select=TRUE)

	### Create the prior for the mean intercept
	b_0_mean_prior <- c(summary(mean_model)$p.table[1], summary(mean_model)$p.table[2])
	b_0_disp_prior <- c(summary(disp_model)$p.table[1], summary(disp_model)$p.table[2])

	### Create a vector for intializing the mean
	b_mean_init <- c(coef(mean_model)[2:length(coef(mean_model))])
	b_disp_init <- c(coef(disp_model)[2:length(coef(disp_model))])

	lambda_mean_init <- unname(c(mean_model$sp ))
	lambda_disp_init <- unname(c(disp_model$sp))

	### Create output list and return
	output_list <- list(b_0 = list(mean = b_0_mean_prior, disp = b_0_disp_prior), b_init = list(mean = b_mean_init, disp = b_disp_init), lambda_init = list(mean = lambda_mean_init, disp = lambda_disp_init))
	return(output_list)

}

