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
create_basis <- function(data, type, knot_loc){
	if (type == "cyclic"){
		output_list <- create_cyclic_basis(data = data, knot_loc = knot_loc)
	} else if (type == "tensor"){
		output_list <- create_tensor_basis(data = data, knot_loc = knot_loc)
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
create_cyclic_basis <- function(data, knot_loc){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)
	
	### Create cyclic spline basis
	spline_orig <- mgcv::smoothCon(s(jdate, bs="cc", k = n_knots), data=data, knots = knot_loc, null.space.penalty = FALSE)

	### Create cyclic spline basis with absorbed constraint
	spline_reparam <- mgcv::smoothCon(s(jdate, bs="cc", k = n_knots), data=data, knots = knot_loc, absorb.cons=TRUE, null.space.penalty = FALSE)

	### Create cyclic spline basis with absorbed constraint
	spline_diag <- mgcv::smoothCon(s(jdate, bs="cc", k = n_knots), data=data, knots = knot_loc, absorb.cons=TRUE, null.space.penalty = FALSE, diagonal.penalty = TRUE)

	### Extract the matrices for basis and penalty term
	X_orig <- spline_orig[[1]]$X

	### Reparameterize the penalty matrix
	s_reparam <- spline_reparam[[1]]$S[[1]]

	### Reparameterize both using the QR decomposition following Wood 
	### Where Z is the Q matrix without the first column, used to reparameterize betas from reparam to original
	C <- rep(1, nrow(X_orig)) %*% X_orig
	qrc <- qr(t(C))
	Z <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]

	### Calculate reparameterized matrices for basis and penalty
	X_reparam <- X_orig%*%Z

	### Extract the matrix for diagonalized basis
	### Incorporates multivariate normal into basis. Only works for single predictor, not used for tensor
	X_diag <- spline_diag[[1]]$X

	### Extract the matrix to transfrom from diagonalized to reparameterized betas
	p <- spline_diag[[1]]$diagRP

	### Create the output and return
	output_list <- list(x_orig = X_orig, x_reparam = X_reparam, x_diag = X_diag, z = Z, s_reparam = s_reparam, c = C, qrc = qrc, p = p)
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
create_tensor_basis <- function(data, knot_loc){
	### Extract the number of knots
	n_knots <- sapply(knot_loc, FUN = length)
	
	### Create cyclic spline basis
	spline_orig <- mgcv::smoothCon(te(jdate,year, bs=c("cc", "cr"), k = n_knots), data=data, knots = knot_loc, null.space.penalty = TRUE)

	### Create cyclic spline basis with absorbed constraint
	spline_reparam <- mgcv::smoothCon(te(jdate,year, bs=c("cc", "cr"), k = n_knots), data=data, knots = knot_loc, absorb.cons=TRUE, null.space.penalty = TRUE)

	### Extract the matrices for basis and penalty term
	X_orig <- spline_orig[[1]]$X

	### Reparameterize the penalty matrix
	s_reparam <- list()
	s_reparam[[1]] <- spline_reparam[[1]]$S[[1]]
	s_reparam[[2]]  <- spline_reparam[[1]]$S[[2]]
	s_reparam[[3]]  <- spline_reparam[[1]]$S[[3]]

	### Reparameterize both using the QR decomposition following Wood 
	### Where Z is the Q matrix without the first column, used to reparameterize
	C <- rep(1, nrow(X_orig)) %*% X_orig
	qrc <- qr(t(C))
	Z <- qr.Q(qrc,complete=TRUE)[,(nrow(C)+1):ncol(C)]

	### Calculate reparameterized matrices for basis and penalty
	X_reparam <- X_orig%*%Z

	### Create the output and return
	output_list <- list(x_orig = X_orig, x_reparam = X_reparam, z = Z, s_reparam = s_reparam, c = C, qrc = qrc)
	return(output_list)
}

