% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{gene_motif_analysis}
\alias{gene_motif_analysis}
\title{Calcualte significant gene-motif combos for celltypes}
\usage{
gene_motif_analysis(
  object,
  motif_gene_df,
  aver_exp_all = NULL,
  aver_chromvar_all = NULL
)
}
\arguments{
\item{object}{Input data in SeuratObject format}

\item{motif_gene_df}{Dataframe of motifs and genes}

\item{aver_exp_all}{Dataframe of fold changes of gene expression}

\item{aver_chromvar_all}{Dataframe of fold changes of motif activity}
}
\description{
Calcualte significant gene-motif combos for celltypes
}
