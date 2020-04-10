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
extract_params <- function(model_fit, basis){

	### Extract basis 
	basis_reparam <- cbind(rep(1,dim(basis$x_reparam)[1]), basis$x_reparam)
	data <- basis$data
	type <- basis$type

	### Extract the spline coefficients and intercept for mean
	b_mean <- extract(model_fit, "b_mean")$b_mean
	b_0_mean <- extract(model_fit, "b_0_mean")$b_0_mean
	### Combine the intercept and spline coefficients into a single matrix
	b_full_mean <- cbind(matrix(b_0_mean, dim(b_mean)[1], 1), b_mean)

	### Extract the spline coefficients and intercept for scale
	b_scale <- extract(model_fit, "b_scale")$b_scale
	b_0_scale <- extract(model_fit, "b_0_scale")$b_0_scale
	### Combine the intercept and spline coefficients into a single matrix
	b_full_scale <- cbind(matrix(b_0_scale, dim(b_scale)[1], 1), b_scale)

	### Calculate the estimate of mean and scale for the demo basis
	b_mean_est <- exp(basis_reparam %*% t(b_full_mean))
	b_scale_est <- exp(basis_reparam %*% t(b_full_scale))

	### Gather the results into a long dataframe
	mean_est <- data.frame(demo_basis$data, b_mean_est) %>%
		pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="mean") 
	scale_est <- data.frame(demo_basis$data, b_scale_est) %>%
		pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="scale") 

	### Add in the derived parameters
	param_est <- full_join(mean_est, scale_est) %>%
		mutate(rate = 1/scale) %>%
		mutate(shape = mean * rate) %>%
		mutate(disp = 1/shape) 

	return(param_est)
}


#' Convert parameter estimate to long format
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
param_long <- function(param_est){
	### Convert to long format to plot all together
	param_est_long <- param_est %>%
		pivot_longer(cols=c(-jdate, -contains("year"), -draw), names_to = "param", values_to="value") 	

	return(param_est_long)
}



#' Convert parameter estimate to long format
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
param_summary <- function(param_est){
	### Convert to long format
	param_est_long <- param_long(param_est)

	### Group by parameter and time
	param_est_summary <- param_est_long %>%
		group_by(param, jdate)

	### For tensor product add a group on year
	if("year" %in% names(param_est_long)) {
		param_est_summary <- param_est_summary %>%
			group_by(year, add=TRUE)
	}

 	### Calculate summary statistics
	param_est_summary <- param_est_summary %>%	
		summarise(median = median(value), 
			mean = mean(value, na.rm=TRUE),
			sd = sd(value, na.rm=TRUE),
			perc_95_lower = quantile(value, 0.025), 
			perc_95_upper = quantile(value, 0.975), 
			perc_50_lower = quantile(value, 0.25), 
			perc_50_upper = quantile(value, 0.75)) %>%
		ungroup()

	return(param_est_summary)
}

