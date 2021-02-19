\name{poolcells}
\alias{poolcells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
pool cells
}
\description{
pool cells from all samples
}
\usage{
poolcells(
eset,
sample
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ExpressionSet object for single cells}
  \item{sample}{variable name for subject/samples}
}
\value{ExpressionSet object for constructed single cell data by pooling all cells}
