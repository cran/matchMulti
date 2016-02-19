\name{assembleMatch}
\alias{assembleMatch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Collect Matched Samples
}
\description{
After students and schools have both been matched separately, assembles the matched student samples corresponding to the school match into a single dataframe of student-level data.
}
\usage{
assembleMatch(student.matches, school.match, school.id, treatment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{student.matches}{
	a list of lists object produced by \code{matchStudents}, with each element of the second list containing a dataframe composed of a matched sample for a different treated-control school pairing.
}
  \item{school.match}{
	a dataframe, produced by \code{matchSchools}, with two columns, one containing treated school IDs and the other containing matched control school IDs.
}
  \item{school.id}{
	the name of the column storing the unique school identifier (in the dataframes stored in \code{student.matches})
}
  \item{treatment}{
	the name of the column storing the binary treatment status indicator (in the dataframes stored in \code{student.matches})
}
}
\value{
a dataframe containing the full set of matched samples for the multilevel match.
}

\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}


\keyword{internal}
