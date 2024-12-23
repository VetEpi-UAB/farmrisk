% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/region_code.R
\name{get_region_code}
\alias{get_region_code}
\title{Regional Code Management Functions}
\usage{
get_region_code(data, module = NULL, pathway = NULL)
}
\arguments{
\item{data}{dataframe with gid0, gid1, gid2 and gid3 data}

\item{module}{Optional module filter}

\item{pathway}{Optional pathway filter}
}
\value{
A dataframe with region_code column
}
\description{
A collection of functions for handling regional codes and hierarchical joins

Get the lowest level (more detailed) region code
}
\examples{
get_region_code(data = mov, module = "purchase")
}
