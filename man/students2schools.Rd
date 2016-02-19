\name{students2schools}
\alias{students2schools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Aggregate Student Data into School Data
}
\description{
Takes a dataframe of student-level covariates and aggregates selected columns into a dataframe of school covariates.
}
\usage{
students2schools(students, school.cov, school.id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{students}{
	a dataframe of students.
}
  \item{school.cov}{
	a character vector of column names in \code{students} that should be aggregated by school.
}
  \item{school.id}{
	the name of the column in \code{students} containing the unique school identifier.
}
}
\details{
	Aggregation is either done by taking averages or by selecting the unique factor value when a school has only one value for a factor.  As a result, \code{school.covs} should only include variables that are numeric or do not vary within schools.	 
}
\value{
	a dataframe of aggregated data, with one row for each school and columns in \code{school.covs} and \code{school.id}.
}
\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

\keyword{internal}