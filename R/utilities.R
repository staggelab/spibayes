
#' Extract draws from a fitdistrplus distribution using standard error variance covariance matrix
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param n Number of draws to calculate
#' @export
draw_fitgamma <- function(data, n){
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

	return(bind_rows(estimate, draws))
}



#' Extract draws from a fitdistrplus distribution using bootstrap resampling
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
draw_fitgamma_boot <- function(data, n){
	data <- data[data > 0 & !is.na(data)]
	fitted <- fitdist(data, "gamma")
	draws <- bootdist(fitted, niter=n)
	draws <- draws$estim
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

	return(bind_rows(estimate, draws))
}


#' Confidence interval on fitdisrplus gamma
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
mle_gamma_fit <- function(jdate, values, n, ci = 0.95, method = "boot"){

	interval <- (1-ci)/2
	interval <- c(interval, 1-interval)

	data <- data.frame(jdate = jdate, values = values)

	### Performing the analysis this way is needed to generate the upper and lower estimates for other parameters
	gamma_draws <- data %>%
		group_by(jdate) %>%
		filter(values > 0) %>%
		filter(jdate <=366) 

	if(method == "boot"){
		gamma_draws <- gamma_draws %>%
		group_modify(~ draw_fitgamma_boot(.x$values, n = n)) %>%
		ungroup()
	} else {
		gamma_draws <- gamma_draws %>%
		group_modify(~ draw_fitgamma(.x$values, n = n)) %>%
		ungroup()
	}

	estimate_df <- gamma_draws %>%
		filter(type == "estimate") %>%
		select(-draw, -type) %>%
		pivot_longer(-jdate, names_to = "param") %>%
		mutate(ci = "estimate")
	
	lower_df <- gamma_draws %>%
		filter(type == "draw") %>%		
		group_by(jdate) %>%
		summarize(mean = quantile(mean, interval[1]), 
			scale = quantile(scale, interval[1]), 
			rate = quantile(rate, interval[1]), 
			shape = quantile(shape, interval[1]), 			
			disp = quantile(disp, interval[1])) %>%
		ungroup() %>%
		pivot_longer(-jdate, names_to = "param") %>%
		mutate(ci = "lower_ci")

	upper_df <- gamma_draws %>%
		filter(type == "draw") %>%	
		group_by(jdate) %>%
		summarize(mean = quantile(mean, interval[2]), 
			scale = quantile(scale, interval[2]), 
			rate = quantile(rate, interval[2]), 
			shape = quantile(shape, interval[2]), 			
			disp = quantile(disp, interval[2])) %>%
		ungroup() %>%
		pivot_longer(-jdate, names_to = "param") %>%
		mutate(ci = "upper_ci")


	ci_df <- estimate_df%>%
		bind_rows(lower_df) %>%
		bind_rows(upper_df)

	gamma_draws <- gamma_draws %>%
		filter(type != "estimate")

	return(list(estimate = ci_df, draws = gamma_draws))
}



