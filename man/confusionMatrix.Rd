% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{confusionMatrix}
\alias{confusionMatrix}
\title{Confusion matrix}
\usage{
confusionMatrix(i = NULL, j = NULL)
}
\arguments{
\item{i}{The first vector}

\item{j}{The second vector}
}
\value{
a matrix of the confusion matrix of the two vectors
}
\description{
This function will make a confusion matrix based on two vectors. I took this
function from the ArchR package.
}
\examples{
\dontrun{
cm <- confusionMatrix(i, j)
}
}
