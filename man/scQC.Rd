% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scQC.R
\name{scQC}
\alias{scQC}
\title{function definitions ##### Common single-cell quality control.
This step will only save the VlnPlot before and after QC but not show.
Users can do this step by themselves using Seurat.}
\usage{
scQC(
  mtx_dir,
  prefix = "./",
  label = "",
  species = "Hs",
  gene_frac = 0.01,
  nFeature = c(200, 5000),
  nCount = 1000,
  mt = 10,
  blank_NTC = FALSE
)
}
\description{
function definitions ##### Common single-cell quality control.
This step will only save the VlnPlot before and after QC but not show.
Users can do this step by themselves using Seurat.
}
