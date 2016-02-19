\name{match2distance}
\alias{match2distance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute School Distance from a Student Match
}
\description{
	Defines a distance between two schools whose students have been matched based on the size of the resulting matched sample and on the student-level covariate balance.
}
\usage{
match2distance(matchFrame, treatFrame, ctrlFrame, student.vars, treatment, largeval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matchFrame}{
	dataframe containing all matched students.
}
  \item{treatFrame}{
	dataframe containing all students from the treated school.
}
  \item{ctrlFrame}{
	dataframe containing all students from the control school.
}
  \item{student.vars}{
	names of variables on which to evaluate balance in the matched sample.  Must be present in the column names of each of \code{matchFrame}, \code{treatFrame} and \code{ctrlFrame}.
}
  \item{treatment}{
	name of the treatment variable. Must be present in the column names of each of \code{matchFrame}, \code{treatFrame} and \code{ctrlFrame}.
}
  \item{largeval}{
	a large penalty value to be added to the distance for each student-level imbalance. 
}
}
\details{
	The distance is computed by (1) subtracting the harmonic mean of the treated and control counts in the matched sample from \code{largeval}  (2) adding \code{largeval} for each covariate among \code{studentvars} that has an absolute standardized difference exceeding 0.2.  This encourages the school match to choose larger schools with better balance.
}
\value{
a numeric distance.
}

\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\keyword{internal}
