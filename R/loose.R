
#' Run model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @param engine  options are sample (full MCMC), variational (variational approximation of posterior),  optimize (Stan's LBFGS algorithm), 
#' @import cmdstanr
#' @return A matrix of the infile
#' @export
disp_to_shape <- function(disp){
	1/ exp(-7 + log(1 + exp(disp)))
}


#' Run model
#'
#' This is where you describe the function itself
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#'
#' @param data Dataframe with the underlying data. Columns must include variable names
#' @param spline_type List of spline types. Either cc or cr are accepted. Names must agree with knot_loc
#' @param knot_loc List of knot locations
#' @param engine  options are sample (full MCMC), variational (variational approximation of posterior),  optimize (Stan's LBFGS algorithm), 
#' @import cmdstanr
#' @return A matrix of the infile
#' @export
logodds_to_p <- function(logodds){
	exp(logodds)/(1+exp(logodds))
}


