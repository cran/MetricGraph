% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_spde.R
\name{graph_spde_basis}
\alias{graph_spde_basis}
\title{Deprecated - Observation/prediction matrices for 'SPDE' models}
\usage{
graph_spde_basis(graph_spde, repl = NULL, drop_na = FALSE, drop_all_na = TRUE)
}
\arguments{
\item{graph_spde}{An \code{inla_metric_graph_spde} object built with the
\code{graph_spde()} function.}

\item{repl}{Which replicates? If there is no replicates, or to
use all replicates, one can set to \code{NULL}.}

\item{drop_na}{Should the rows with at least one NA for one of the columns be removed? DEFAULT is \code{FALSE}.}

\item{drop_all_na}{Should the rows with all variables being NA be removed? DEFAULT is \code{TRUE}.}
}
\value{
The observation matrix.
}
\description{
Constructs observation/prediction weight matrices
for metric graph models.
}
