\name{matchSchools}
\alias{matchSchools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Match Schools on Student-based Distance
}
\description{
Takes in a school distance matrix created using information from the first-stage student match and matches schools optimally, potentially 
}
\usage{
matchSchools(dmat, students, treatment, school.id, school.fb, penalty, verbose, tol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{
	a distance matrix for schools, with a row for each treated school and a column for each control school.
}
  \item{students}{
	a dataframe containing student and school covariates, with a different row for each student.
}
  \item{treatment}{
	the column name of the binary treatment status indicator in the \code{students} dataframe.
}
  \item{school.id}{
	the column name of the unique school ID in the \code{students} dataframe.
}
  \item{school.fb}{
	an optional list of character vectors, each containing a subset of the column names of \code{students}.  Each element of the list should contain all the names in previous elements (producing a nested structure).
}
  \item{penalty}{
	a numeric value, treated as the cost to the objective function of excluding a treated school.  If it is set lower, more schools will be excluded.
}
  \item{verbose}{
	a logical value indicating whether detailed output should be printed.
}

\item{tol}{
	a numeric tolerance value for comparing distances.  It may need to be raised above the default when matching with many levels of refined balance.
}

}
\details{
	The \code{school.fb} argument encodes a refined covariate balance constraint: the matching algorithm optimally balances the interaction of the variables in the first list element, then attempts to further balance the interaction in the second element, and so on.  As such variables should be added in order of priority for balance.
}
\value{
	a dataframe with two columns, one containing treated school IDs and the other containing matched control school IDs.
}

\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
