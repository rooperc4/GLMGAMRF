% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DesignBased.R
\name{get_strata_area}
\alias{get_strata_area}
\title{Function to get the strata area for a stratum}
\usage{
get_strata_area(data, strata_column, region)
}
\arguments{
\item{data}{data set where the stratum area is needed}

\item{strata_column}{Strata designation for each haul}

\item{region}{Either "GOA" or "AI"}
}
\description{
This function reads in the GOA and AI stratum tables so that a stratum area can be designated 
for each bottom trawl haul where the stratum is known.
}
\examples{
Design.data<-get_strata_area(Juvenile_POP_Data,"STRATUM","GOA")
}
\keyword{random}
\keyword{stratified}
\keyword{survey}
