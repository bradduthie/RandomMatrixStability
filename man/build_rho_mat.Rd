% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{build_rho_mat}
\alias{build_rho_mat}
\title{Build a matrix with a specified correlation structure}
\usage{
build_rho_mat(S, sigma, rho, mn = 0, dval = 1)
}
\arguments{
\item{S}{The size of the network (number of components)}

\item{sigma}{The standard deviation of network interaction strengths}

\item{rho}{The correlation between component interaction strengths (i.e.,
between off-diagonal matrix elements A_{ij} and A_{ji}}

\item{mn}{The mean interaction strength}

\item{dval}{Self-regulation of network elements (1 by default)}
}
\value{
A matrix with a pre-specified correlation structure.
}
\description{
Builds a random network with a pre-specified correlation structure between 
off-diagonal elements A_{ij} and A_{ji}
}
\examples{
eg_mat <- build_rho_mat(S = 8, sigma = 0.4, rho = 0.2);
}
