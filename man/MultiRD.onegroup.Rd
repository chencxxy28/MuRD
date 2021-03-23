\name{MultiRD.onegroup}
\alias{MultiRD.onegroup}
\title{MultiRD algorithm to deconvolute cell type proportions from a target bulk data}
\description{
This is the main algorithm that MultiRD is relied on to implement deconvolution.
}
\usage{
MultiRD.onegroup(
bulk.data,
list.marker,
celltype.unique,
subject.level.proportion,
population.level.proportion,
proportion.sd=1,
lambda.option=c(seq(from=0,to=0.075,length=15),10,50,100,500,1000),
tol=0.001,
iter.num=1000
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bulk.data}{ExpressionSet object for a target bulk data}
  \item{list.marker}{A list of pre-specified marker genes corresponding to each cell type}
  \item{celltype.unique}{A list of cell types. It should match the order in list.marker}
  \item{subject.level.proportion}{A pre-specified cell type proportions for the target bulk data, which could be obtained from reference-based deconvolution approach.}
  \item{population.level.proportion}{A pre-specified population-level cell type proportions, which could be obtained from single-cell RNA-seq and external expression data from different studies, species, or data types}
  \item{proportion.sd}{an optional adjustment based on pre-specified standard deviation of cell-type proportion estimation. The default is 1 for each cell type.}
  \item{lambda.option}{a sequence of values for the tunning parameter}
  \item{tol}{a small constant used for convergence criteria. The default is 0.001}
  \item{iter.num}{The maximum number of iteration. The default is 1000.}
}
\value{
\item{est.prop}{a list containing estimated cell type proportions corresponding to each tuning value.}
\item{metrics}{a sequence of goodness-of-fit values corresponding to  each tuning value. The smaller the better.}
}
