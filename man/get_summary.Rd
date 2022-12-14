% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_data.R
\name{get_summary}
\alias{get_summary}
\title{Summarize data}
\usage{
get_summary(data, generation, trait)
}
\arguments{
\item{data}{A data frame that contains the generation column and trait column.}

\item{generation}{A vector that contains the generations.}

\item{trait}{A numeric vector that contains the trait value.}
}
\value{
A data frame that summarizes the number of observations, mean, mean
of variance, model.
}
\description{
Summarize data
}
\examples{
get_summary(data, generation, trait)
}
