\name{estimate.geneprofile}
\alias{estimate.geneprofile}
\title{
Estimate mean gene profiles
}
\description{
Mean gene profiles are estimated given estimated proportions
}
\usage{
estimate.geneprofile(
bulk.data,
gene.used,
celltype.unique,
cluster.identifier,
prop.est
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bulk.data}{ExpressionSet object for a given bulk data}
  \item{gene.used}{Pool of marker genes to be estimated}
  \item{celltype.unique}{A list of cell types. It should match the order in list.marker}
  \item{cluster.identifier}{Indentifiers to indicate which cell type a given marker gene belongs to. They are coded as 1,2,3,...,T to represent T cell types, respectively}
  \item{prop.est}{Estimated cell type proprotion matrix}
}
\value{
This function returns a matrix containing estimated mean gene profiles for all cell types, with rows indicating genes and columns indicating cell types.
}
