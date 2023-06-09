% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_inference.R
\name{posterior_crossvalidation}
\alias{posterior_crossvalidation}
\title{Leave-one-out crossvalidation for \code{graph_lme} models assuming observations at
the vertices of metric graphs}
\usage{
posterior_crossvalidation(object)
}
\arguments{
\item{object}{A fitted model using the \code{graph_lme()} function.}
}
\value{
Vector with the posterior expectations and variances as well as
mean absolute error (MAE), root mean squared errors (RMSE), and three
negatively oriented proper scoring rules: log-score, CRPS, and scaled
CRPS.
}
\description{
Leave-one-out crossvalidation for \code{graph_lme} models assuming observations at
the vertices of metric graphs
}
