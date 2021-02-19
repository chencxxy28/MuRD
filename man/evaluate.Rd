\name{evaluate}
\alias{evaluate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
evaluation function for estimated cell type proportions
}
\description{
Several evaluation metrics are provided, such as mean absolute deviance (MAD), concordance correlation coefficient (CCC), Pearson correlation coefficient (Pearson), given true cell type proportions,
}
\usage{
evaluate(
est.prop,
true.prop
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est.prop}{estimated cell type proportions}
  \item{true.prop}{true cell type proportions}
}
\value{
\item{cell.type.eva}{cell-type level evaluations based on MAD, CCC, and Pearson. }
\item{all.eva}{Overall evaluations based on averaged MAD, CCC, and Pearson.}
}
