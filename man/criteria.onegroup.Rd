\name{criteria.onegroup}
\alias{criteria.onegroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Assess a goodness-of-fit metric
}
\description{
A goodness-of-fit metric is calculated to assess the performance of given proportion estimation
}
\usage{
criteria.onegroup(
bulk.data,
prop.used
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bulk.data}{ExpressionSet object for a test bulk data}
  \item{prop.used}{A given proportion estimation}
}
\value{An evaluation based on three metrics, that is, median, mean absolute deviance (MAD), and root of mean square deviance (RMSD). The latter two are more recommended to use in practice.}
