% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{spatrio_heatmap}
\alias{spatrio_heatmap}
\title{Visualization of spatrio results using heatmap}
\usage{
spatrio_heatmap(
  input_list,
  group = NULL,
  main = NA,
  cols = c("#00578C", "#8ED6EC", "white", "#FF7B67", "#B00028"),
  ann_cols = NA,
  name = "Module",
  ...
)
}
\arguments{
\item{input_list}{Spatrio module analysis results}

\item{group}{Modules to be plotted}

\item{main}{The title of the plot}

\item{cols}{Vector of colors used in heatmap}

\item{ann_cols}{Group color}

\item{name}{Prefix of module serial number}
}
\description{
Visualization of spatrio results using heatmap
}
