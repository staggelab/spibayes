% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{predict_vals}
\alias{predict_vals}
\title{Extract model output parameters}
\usage{
predict_vals(model_fit, newdata = NULL, saved_model = NULL)
}
\arguments{
\item{model_fit}{spline object created from fit_model}

\item{newdata}{dataframe where you would like new estimates. Defaults to NULL or}

\item{saved_model}{Simple gives back only the parameters, Marginal gives each marginal}
}
\value{
A matrix of the infile
}
\description{
This is where you describe the function itself
contains the rownames and the subsequent columns are the sample identifiers.
Any rows with duplicated row names will be dropped with the first one being
}
