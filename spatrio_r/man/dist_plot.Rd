% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{dist_plot}
\alias{dist_plot}
\title{Visualization of distance}
\usage{
dist_plot(
  object,
  coord = object@images$image@coordinates[, c("x", "y")],
  root = "Cell type 1",
  level = c("Cell type 1", "Cell type 2", "Cell type 3")
)
}
\arguments{
\item{object}{Input data in SeuratObject format}

\item{coord}{Matrix of spatial coordinate}

\item{root}{The cluster id used as root when calculating distance}

\item{level}{The level of cell clusters}
}
\description{
Visualization of distance
}
