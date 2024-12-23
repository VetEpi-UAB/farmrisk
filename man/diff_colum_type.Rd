% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidy_custom_table.R
\name{diff_colum_type}
\alias{diff_colum_type}
\title{Column type comparison checker}
\usage{
diff_colum_type(data.x, data.y)
}
\arguments{
\item{data.x}{First comparison dataframe}

\item{data.y}{Second comparison dataframe}
}
\value{
Detailed report of column type differences
}
\description{
Identifies column type mismatches between two dataframes that could cause binding issues.
}
\examples{
diff_colum_type(data.x = visit_veh_data, data.y = visit_people_data)
}
