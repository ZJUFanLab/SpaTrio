% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{spatrio_score}
\alias{spatrio_score}
\title{Calcualte the score of modules}
\usage{
spatrio_score(
  object,
  input_list,
  nbin = 10,
  ctrl = 50,
  seed = 42,
  clean = T,
  name = "Module"
)
}
\arguments{
\item{object}{Input data in SeuratObject format}

\item{input_list}{Modules to be plotted}

\item{nbin}{Number of bins of aggregate expression levels for all analyzed feature}

\item{ctrl}{Number of control features}

\item{seed}{Random seed}

\item{clean}{Whether to clear the previous score}

\item{name}{Prefix of module serial number}
}
\description{
Calcualte the score of modules
}
