% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_backward.R
\name{mov_backward}
\alias{mov_backward}
\title{Movements Survey Backward Compatibility}
\usage{
mov_backward(data)
}
\arguments{
\item{data}{A data frame containing mov survey data}
}
\value{
A modified data frame compatible with the latest version
}
\description{
Manages compatibility between different versions of mov survey data,
applying necessary transformations to maintain consistency with the latest version.
}
\details{
Handle backward compatibility for movement survey versions
}
\examples{
mov <- read.csv("mov.csv")
compatible_data <- mov_backward(mov)
}
