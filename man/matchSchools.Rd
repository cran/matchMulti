% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchMulti-internal.R
\name{matchSchools}
\alias{matchSchools}
\title{Match Schools on Student-based Distance}
\usage{
matchSchools(
  dmat,
  students,
  treatment,
  school.id,
  school.fb,
  penalty,
  verbose,
  tol,
  solver = "rlemon"
)
}
\arguments{
\item{dmat}{a distance matrix for schools, with a row for each treated
school and a column for each control school.}

\item{students}{a dataframe containing student and school covariates, with a
different row for each student.}

\item{treatment}{the column name of the binary treatment status indicator in
the \code{students} dataframe.}

\item{school.id}{the column name of the unique school ID in the
\code{students} dataframe.}

\item{school.fb}{an optional list of character vectors, each containing a
subset of the column names of \code{students}.  Each element of the list
should contain all the names in previous elements (producing a nested
structure).}

\item{penalty}{a numeric value, treated as the cost to the objective
function of excluding a treated school.  If it is set lower, more schools
will be excluded.}

\item{verbose}{a logical value indicating whether detailed output should be
printed.}

\item{tol}{a numeric tolerance value for comparing distances.  It may need
to be raised above the default when matching with many levels of refined
balance.}

\item{solver}{Name of package used to solve underlying network flow problem, 
one of 'rlemon' and 'rrelaxiv'.  rrelaxiv carries an 
academic license and is not hosted on CRAN so it must be installed 
separately.}
}
\value{
a dataframe with two columns, one containing treated school IDs and
the other containing matched control school IDs.
}
\description{
Takes in a school distance matrix created using information from the
first-stage student match and matches schools optimally, potentially
}
\details{
The \code{school.fb} argument encodes a refined covariate balance
constraint: the matching algorithm optimally balances the interaction of the
variables in the first list element, then attempts to further balance the
interaction in the second element, and so on.  As such variables should be
added in order of priority for balance.
}
\author{
Luke Keele, Penn State University, \email{ljk20@psu.edu}

Sam Pimentel, University of California, Berkeley, \email{spi@berkeley.edu}
}
\keyword{internal}
