% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_input_files.R
\name{check_file_modifications}
\alias{check_file_modifications}
\title{Check if files need to be processed from JSON}
\usage{
check_file_modifications(user_path, user_files_info)
}
\arguments{
\item{user_path}{Character string. Path to user directory}

\item{user_files_info}{File information from file.info()}
}
\value{
Logical. TRUE if JSON processing is needed
}
\description{
Check if files need to be processed from JSON
}
