% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L1_graph.R
\name{get_knn}
\alias{get_knn}
\title{function to find the k-nearest neighbor graph}
\usage{
get_knn(X, K = 5)
}
\arguments{
\item{X}{number of rows in the returned eye matrix (D * N)}

\item{K}{number of neighbors of each data point excluding itself}
}
\value{
a matrix
}
\description{
function to find the k-nearest neighbor graph
}
