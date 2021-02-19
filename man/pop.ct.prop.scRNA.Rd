\name{pop.ct.prop.scRNA}
\alias{pop.ct.prop.scRNA}
\title{
calculate population-level cell type proportions from single-cell data
}
\description{
calculate population-level cell type proportions from single-cell data
}
\usage{
pop.ct.prop.scRNA(
scRNA,
cluster,
sample,
sep.group=F,
group=NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scRNA}{ExpressionSet object for single-cell data}
  \item{cluster}{character string specifying the variable name for cell types}
  \item{sample}{character string specifying the variable name for subject/samples}
  \item{sep.group}{logical. Whether to seperately calculate population-level cell type proportions by groups. The default is False.}
  \item{group}{character string specifying the variable name for group}
}
\value{population-level cell type proportions}
