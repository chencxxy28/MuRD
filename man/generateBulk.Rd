\name{generateBulk}
\alias{generateBulk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{bulk data generation function
}
\description{
Construct Pseudo bulk samples by random sampled number of cells per subject, but based on the actual numbers.
}
\usage{
generateBulk(
eset,
ct.varname,
sample,
disease = NULL,
ct.sub,
prop_mat = NULL,
nbulk = 50,
samplewithRep = F,
low_s = 0.3,
upp_s = 0.7
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ExpressionSet object for single cells}
  \item{ct.varname}{variable name for 'cell types'}
  \item{sample}{variable name for subject/samples}
  \item{disease}{indicate the health condition of subjects}
  \item{ct.sub}{a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.}
  \item{prop_mat}{manually input the cell-type proportion for pseudo bulk samples}
  \item{nbulk}{number of pseudo bulk samples to be constructed}
  \item{samplewithRep}{logical, randomly sample single cells with replacement. Default is F.}
  \item{low_s}{lower supprot a for uniform distribution U[a,b]}
  \item{upp_s}{upper supprot b for uniform distribution U[a,b]}
}
\value{
pseudo bulk samples ExpressionSet, and actual cell-type proportions
}
