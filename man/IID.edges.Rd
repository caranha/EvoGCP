% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_culberson.R
\name{IID.edges}
\alias{IID.edges}
\title{IID - Independent Random Edge Assignment
Source: https://webdocs.cs.ualberta.ca/~joe/Coloring/Generators/manual.html}
\usage{
IID.edges(V1, V2, V3, probability)
}
\arguments{
\item{V1}{Vector - group of nodes with color 1.}

\item{V2}{Vector - group of nodes with color 2.}

\item{V3}{Vector - group of nodes with color 3.}

\item{probability}{0..1 - edge's probability.}
}
\value{
a list of edges.
}
\description{
Assigns edges between nodes from different groups with a given probability
}
