% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ga.R
\name{get_top_evo_out}
\alias{get_top_evo_out}
\title{Find the component response rate vectors of highest size}
\usage{
get_top_evo_out(evo_out, size)
}
\arguments{
\item{evo_out}{Output of component response values}

\item{size}{The number of gamma vectors to be returned}
}
\value{
A list of vectors of component response rates (gammas)
}
\description{
Takes the output of the genetic algorithm showing the component response rate
(gamma) vectors found to be stabilising and returns the ones that are from 
the largest systems (i.e., highest S)
}
