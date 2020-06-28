
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
extract_params <- function(model_fit, basis, newdata = NULL){
	### Extract type
	type <- basis$type

	### Extract the spline coefficients and intercept for mean
	b_mean <- extract(model_fit, "b_mean")$b_mean
	b_0_mean <- extract(model_fit, "b_0_mean")$b_0_mean
	### Convert to original basis
	b_mean_orig <- t(apply(b_mean, 1, function(x){basis$z %*% x}))
	### Combine the intercept and spline coefficients into a single matrix
	b_full_mean <- cbind(matrix(b_0_mean, dim(b_mean_orig)[1], 1), b_mean_orig)

	### Extract the spline coefficients and intercept for scale
	b_scale <- extract(model_fit, "b_scale")$b_scale
	b_0_scale <- extract(model_fit, "b_0_scale")$b_0_scale
	### Convert to original basis
	b_scale_orig <- t(apply(b_scale, 1, function(x){basis$z %*% x}))
	### Combine the intercept and spline coefficients into a single matrix
	b_full_scale <- cbind(matrix(b_0_scale, dim(b_scale_orig)[1], 1), b_scale_orig)

	### Extract theta
	b_theta <- extract(model_fit, "b_theta")$b_theta
	b_0_theta <- extract(model_fit, "b_0_theta")$b_0_theta
	### Convert to original basis
	b_theta_orig <- t(apply(b_theta, 1, function(x){basis$z %*% x}))
	### Combine the intercept and spline coefficients into a single matrix
	b_full_theta <- cbind(matrix(b_0_theta, dim(b_theta_orig)[1], 1), b_theta_orig)

	### If new data isn't provided, assume it is the data used to fit the model
	if(is.null(newdata)){
		x_orig <- cbind(1, basis$x_orig)
		newdata <- basis$data%>%
			mutate(data = "original") 
	} else {
		newdata <- newdata %>%
			mutate(data = "newdata") %>%
			bind_rows(mutate(basis$data, data = "original"))

		basis_newdata <- create_basis(data = newdata, type = basis$type, knot_loc = basis$knot_loc)
		x_orig <- cbind(1, basis_newdata$x_orig)

		newdata_test <- newdata$data == "newdata"

		x_orig <- x_orig[newdata_test, ]
		newdata <- newdata[newdata_test,]
	}

	### Calculate the estimate of mean and scale for the demo basis
	b_mean_est <- x_orig %*% t(b_full_mean)
	b_scale_est <- x_orig %*% t(b_full_scale)
	b_scale_est <- x_orig %*% t(b_theta_scale)

	if (type == "cyclic"){
		mean_est <- data.frame(select(newdata, "jdate", "year"), b_mean_est) %>%
			pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="mean") 
		scale_est <- data.frame(select(newdata, "jdate", "year"), b_scale_est)  %>%
			pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="scale") 
		theta_est <- data.frame(select(newdata, "jdate", "year"), theta_est) %>%
			pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="theta_logodds") 
	} else {
	### Gather the results into a long dataframe
		mean_est <- data.frame(select(newdata, -"data", -"date", -"precip"), b_mean_est) %>%
			#select("jdate", "year", "mean") %>%
			pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="mean") 
		scale_est <- data.frame(select(newdata, -"data", -"date", -"precip"), b_scale_est)  %>%
			pivot_longer(cols=c(-jdate, -contains("year")),  names_to = "draw", values_to="scale") 
	}

	### Add in the derived parameters
	param_est <- full_join(mean_est, scale_est) %>%
		full_join(theta_est) %>%
		mutate(rate = 1/scale) %>%
		mutate(shape = mean * rate) %>%
		mutate(disp = 1/shape) %>%
		mutate(theta = exp(theta_logodds)/(1+exp(theta_logodds)))

	return(list(param_est = param_est, b_mean = b_full_mean, b_scale = b_full_scale,  b_theta = b_full_theta))
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

