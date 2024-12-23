% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_model.R
\name{get_previous_nodes}
\alias{get_previous_nodes}
\title{Get Nodes from Previous Module}
\usage{
get_previous_nodes(mcmodule, get_nodes = NULL)
}
\arguments{
\item{mcmodule}{An mcmodule or mcnode_list object}

\item{get_nodes}{Optional vector of node names to retrieve}
}
\value{
A subset of the node list containing requested nodes
}
\description{
Retrieves nodes from a previous module and assigns them to the parent environment
}
\examples{
get_previous_nodes(mcmodule)
}
