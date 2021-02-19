\name{testset}
\alias{testset}
\title{
 A test data set generation
}
\description{
Create the test data set based on the target bulk data and a list of marker genes
}
\usage{
testset(
eset,
list.marker
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{ExpressionSet object for the target bulk data}
  \item{list.marker}{A list of pre-specified marker genes corresponding to each cell type}
}
\value{A new ExpressionSet object served as a test data set}
