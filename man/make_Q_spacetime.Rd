% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/space_time.R
\name{make_Q_spacetime}
\alias{make_Q_spacetime}
\title{Space-time precision operator discretization}
\usage{
make_Q_spacetime(graph, t, kappa, rho, gamma, alpha, beta, sigma)
}
\arguments{
\item{graph}{A \code{metric_graph} object.}

\item{t}{Vector of time points.}

\item{kappa}{Spatial range parameter.}

\item{rho}{Drift parameter.}

\item{gamma}{Temporal range parameter.}

\item{alpha}{Smoothness parameter (integer) for spatial operator.}

\item{beta}{Smoothness parameter (integer) for Q-Wiener process.}

\item{sigma}{Variance parameter.}
}
\value{
Precision matrix.
}
\description{
The precision matrix for all vertices for space-time field.
}
