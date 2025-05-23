% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spde_covariance.R
\name{spde_covariance}
\alias{spde_covariance}
\title{Covariance function for Whittle-Matérn fields}
\usage{
spde_covariance(P, kappa, tau, range, sigma, alpha, graph, directional = F)
}
\arguments{
\item{P}{Location (edge number and normalized location on the edge) for the
location to evaluate the covariance function at.}

\item{kappa}{Parameter kappa from the SPDE.}

\item{tau}{Parameter tau from the SPDE.}

\item{range}{Range parameter.}

\item{sigma}{Standard deviation parameter.}

\item{alpha}{Smoothness parameter (1 or 2).}

\item{graph}{A \code{metric_graph} object.}

\item{directional}{bool is the model a directional or not. directional only works for alpha=1}
}
\value{
Vector with the covariance function evaluate at the mesh locations.
}
\description{
Computes the covariance function for a Whittle-Matérn field.
}
\details{
Compute the covariance function \eqn{\rho(P,s_i)}{\rho(P,s_i)} where
P is the provided location and \eqn{s_i}{s_i} are all locations in the mesh
of the graph.
}
