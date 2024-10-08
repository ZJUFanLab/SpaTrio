% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatrio.R
\name{spatrio}
\alias{spatrio}
\title{Spatrio analysis using optimal transport}
\usage{
spatrio(
  input_path = getwd(),
  ref_path = NULL,
  output_path = getwd(),
  py_path = NULL,
  spatrio_path = NULL,
  marker_use = TRUE,
  top_marker_num = 100,
  hvg_use = FALSE,
  alpha = 0.1,
  dissimilarity = "scaled_euc",
  k = 10,
  graph_mode = "connectivity",
  aware_spatial = TRUE,
  aware_multi = TRUE,
  aware_power = 2,
  top_num = 5,
  random = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{input_path}{Path to load input data}

\item{output_path}{Path to save output data}

\item{py_path}{Path to python interpreter}

\item{spatrio_path}{Path to spatrio python package}

\item{marker_use}{Determines whether to select differential genes of each cell/spot type for subsequent analysis}

\item{top_marker_num}{Number of differential genes in each cell/spot type}

\item{hvg_use}{Determines whether to only reserve highly variable genes}

\item{alpha}{Alignment tuning parameter. Note: 0 <= alpha <= 1}

\item{dissimilarity}{Expression dissimilarity measure}

\item{k}{Number of neighbors to be used when constructing kNN graph}

\item{graph_mode}{Determines whether to use a connectivity graph or a distance graph}

\item{aware_spatial}{Determines whether to adjust the distance between spots according to areas (or other meta info)}

\item{aware_multi}{Determines whether to adjust the distance between cells according to types (or other meta info)}

\item{aware_power}{Type aware parameter. The greater the parameter, the greater the distance between different areas/types of spots/cells in the graph}

\item{random}{Determines whether to randomly assign cell coordinates or assign coordinates based on pearson correlation coefficient}

\item{verbose}{}

\item{The}{maximum number of cells allocated in a spot}
}
\description{
Spatrio analysis using optimal transport
}
