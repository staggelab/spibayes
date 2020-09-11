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
create_basis <- function(data, type, knot_loc, newdata = FALSE){

	### Catch - if there is no precip column, add a placeholder
	if(!("precip" %in% colnames(data))){
		cat("No precip column in data")

		data <- data %>%
			mutate(precip = 1)
	}

	### Create a column for zero precip
	data <- data %>%
		mutate(zero = precip == 0)
	
	### Remove any leap days for fitting purposes
	data <- data %>%
		filter(jdate <=365)

	### Send to appropriate basis function
	if (type == "cyclic"){
		output_list <- create_cyclic_basis(data = data, knot_loc = knot_loc, newdata = newdata)
	} else if (type == "tensor"){
		output_list <- create_tensor_basis(data = data, knot_loc = knot_loc, newdata = newdata)
	}

	### Create the output and return
	output_list$data <- data
	output_list$type <- type
	output_list$knot_loc <- knot_loc
	return(output_list)
}


#' Create Cyclic Spline Basis
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
create_cyclic_basis <- function(data, knot_loc, newdata){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)

	### Extract positive precipitation
	data_pos <- data %>%
		filter(precip > 0)

	### Check for newdata and whether there is a precipitation column
	### If false, give back the data with gaps for precipitation of zero
	### Else, use the gaps for fitting, but use full dataset for basis
	if( newdata == FALSE) {
		newdata <- data
		newdata_pos <- data_pos
	} else {
		newdata <- newdata
		newdata_pos <- newdata
	}

	### Create cyclic spline basis for precipitation
	gamma_model <- mgcv::gam(list(precip 
		~ ti(jdate, bs=c("cc"), k = n_knots[1]),
		~ ti(jdate, bs=c("cc"), k = n_knots[1])),
		data = data_pos, 
		knots = list(jdate = knot_loc$jdate), 
		family=gammals(link=list("identity","log")),
		fit = FALSE, 
		select = TRUE)

	### Create cyclic spline basis for zero precipitation
	theta_model <- mgcv::gam(zero 
		~ ti(jdate, bs=c("cc"), k = c(n_knots[1])), 
		data = data, 
		knots = list(jdate = knot_loc$jdate), 
		family=binomial,
		fit = FALSE, 
		select = TRUE)

	### Extract basis matrix for mean
	X_mean_jdate <- PredictMat(gamma_model$smooth[[1]],newdata_pos)

	### Extract basis matrix for dispersion
	X_disp_jdate <- PredictMat(gamma_model$smooth[[2]],newdata_pos)

	### For zeros
	X_theta_jdate <- PredictMat(gamma_model$smooth[[1]],newdata)

	### Extract S Penalty matrices for positive
	s_mean_jdate <- gamma_model$S[[1]]
	s_disp_jdate <- gamma_model$S[[2]]

	### Extract S Penalty matrices for zero precip
	s_theta_jdate <- gamma_model$S[[1]]

	### Extract S scaling function
	s_scale_pos <- sapply(gamma_model$smooth, function(x){x$S.scale})
	s_scale_zero <- sapply(gamma_model$smooth, function(x){x$S.scale})
	
	### Extract basis dimensions
	basis_dim_mean <- dim(X_mean_jdate)[2]
	basis_dim_disp <- dim(X_disp_jdate)[2]
	basis_dim_theta <- dim(X_theta_jdate)[2]
	
	### Create the output and return
	output_list <- list(x_mean = list(jdate = X_mean_jdate), 
		x_disp = list(jdate = X_disp_jdate),
		x_theta = list(jdate = X_theta_jdate),
		s_mean = list(jdate = s_mean_jdate),
		s_disp = list(jdate = s_disp_jdate),
		s_theta = list(jdate = s_theta_jdate),
		s_scale = list(pos = s_scale_pos, zero = s_scale_zero),
		basis_dim = list(mean = basis_dim_mean, disp = basis_dim_disp, theta = basis_dim_theta),
		data_all = data,
		data_pos = data_pos,
		gamma_model = gamma_model,
		theta_model = theta_model
	)

	### Return the object
	return(output_list)
}

#' Create Tensor Spline Basis
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
create_tensor_basis <- function(data, knot_loc, newdata){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)

	### Extract positive precipitation
	data_pos <- data %>%
		filter(precip > 0)

	### Check for newdata
	### If false, give back the data with gaps for precipitation
	### Else, use the gaps for fitting, but use full dataset for basis
	if( newdata == FALSE & sum(names(data) %in% "precip") > 0) {
		newdata <- data
		newdata_pos <- data_pos
	} else {
		newdata <- newdata
		newdata_pos <- newdata
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
	X_mean_jdate <- PredictMat(gamma_model$smooth[[1]],newdata_pos)
	X_mean_year <- PredictMat(gamma_model$smooth[[2]],newdata_pos)
	X_mean_tensor <- PredictMat(gamma_model$smooth[[3]],newdata_pos)

	### Extract basis matrix for dispersion
	X_disp_jdate <- PredictMat(gamma_model$smooth[[4]],newdata_pos)
	X_disp_year <- PredictMat(gamma_model$smooth[[5]],newdata_pos)
	X_disp_tensor <- PredictMat(gamma_model$smooth[[6]],newdata_pos)

	### For zeros
	X_theta_jdate <- PredictMat(theta_model$smooth[[1]],newdata)
	X_theta_year <- PredictMat(theta_model$smooth[[2]],newdata)
	X_theta_tensor <- PredictMat(theta_model$smooth[[3]],newdata)

	### Extract S Penalty matrices for positive
	s_mean_jdate <- gamma_model$S[[1]]
	s_mean_year <- gamma_model$S[[2]]
	s_mean_year_double <- gamma_model$S[[3]]
	s_mean_tensor_jdate <- gamma_model$S[[4]]
	s_mean_tensor_year <- gamma_model$S[[5]]

	s_disp_jdate <- gamma_model$S[[6]]
	s_disp_year <- gamma_model$S[[7]]
	s_disp_year_double <- gamma_model$S[[8]]
	s_disp_tensor_jdate <- gamma_model$S[[9]]
	s_disp_tensor_year <- gamma_model$S[[10]]

	### Extract S Penalty matrices for zero precip
	s_theta_jdate <- theta_model$S[[1]]
	s_theta_year <- theta_model$S[[2]]
	s_theta_year_double <- theta_model$S[[3]]
	s_theta_tensor_jdate <- theta_model$S[[4]]
	s_theta_tensor_year <- theta_model$S[[5]]

	### Extract S scaling function
	s_scale_pos <- sapply(gamma_model$smooth, function(x){x$S.scale})
	s_scale_zero <- sapply(theta_model$smooth, function(x){x$S.scale})
	
	### Extract basis dimensions
	basis_dim_mean <- c(dim(X_mean_jdate)[2], dim(X_mean_year)[2], dim(X_mean_tensor)[2])
	basis_dim_disp <- c(dim(X_disp_jdate)[2], dim(X_disp_year)[2], dim(X_disp_tensor)[2])
	basis_dim_theta <- c(dim(X_theta_jdate)[2], dim(X_theta_year)[2], dim(X_theta_tensor)[2])
	
	### Create the output and return
	output_list <- list(
		x_mean = list(jdate = X_mean_jdate, year = X_mean_year, tensor = X_mean_tensor),
		x_disp = list(jdate = X_disp_jdate, year = X_disp_year, tensor = X_disp_tensor),
		x_theta = list(jdate = X_theta_jdate, year = X_theta_year, tensor = X_theta_tensor),
		s_mean = list(jdate = s_mean_jdate, year = s_mean_year, year_double = s_mean_year_double, tensor_jdate = s_mean_tensor_jdate , tensor_year =s_mean_tensor_year ),
		s_disp = list(jdate = s_disp_jdate, year = s_disp_year, year_double = s_disp_year_double, tensor_jdate = s_disp_tensor_jdate , tensor_year =s_disp_tensor_year ),
		s_theta = list(jdate = s_theta_jdate, year = s_theta_year, year_double = s_theta_year_double, tensor_jdate = s_theta_tensor_jdate , tensor_year =s_theta_tensor_year ),
		s_scale = list(pos = s_scale_pos, zero = s_scale_zero),
		basis_dim = list(mean = basis_dim_mean, disp = basis_dim_disp, theta = basis_dim_theta),
		data_all = data,
		data_pos = data_pos,
		gamma_model = gamma_model,
		theta_model = theta_model
	)

	### Return the object
	return(output_list)

}

