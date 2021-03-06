% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solver_ffa.R
\name{ffa.heuristical_swap}
\alias{ffa.heuristical_swap}
\title{Heuristical Swap (step2)}
\usage{
ffa.heuristical_swap(p, w, v, satur)
}
\arguments{
\item{p}{vector of colors}

\item{w}{vector of weights}

\item{v}{the uncolored vertex which will be swapped}

\item{satur}{vector of saturation degrees}
}
\value{
a permutation of w
}
\description{
Swaps an uncolored vertex with the vertex that has highest saturation degree
if tie, choose randomly
}
\details{
Steps:
1. get the first uncolored vertex
2. sort the predecessors according to the saturation degree descending
3. swap the uncolored vertex with the vertex that has highest saturation degree, if tie, choose randomly (from the tie set)
}
