% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_C}
\alias{get_C}
\title{Get the connectance of a network}
\usage{
get_C(mat)
}
\arguments{
\item{mat}{The matrix to be valuated}
}
\value{
The complexity of the matrix mat
}
\description{
Returns the connectance of a network (square matrix), defined as the 
proportion of non-zero off-diagonal elements.
}
\examples{
eg_mat  <- matrix(dat = rnorm(n = 16), nrow = 4);
mat_cmp <- get_C(mat = eg_mat);
}
