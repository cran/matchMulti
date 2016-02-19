\name{matchStudents}
\alias{matchStudents}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Compute Student Matches for all Pairs of Schools
}
\description{
Iterates over all possible treated-control school pairs, optionally computes and stores an optimal student match for each one, and generates a distance matrix for schools based on the quality of each student match.
}
\usage{
matchStudents(students, treatment, school.id, match.students,
student.vars, verbose, penalty.qtile, min.keep.pctg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{students}{
	a dataframe containing student covariates, with a different row for each student.

}
  \item{treatment}{
	the column name of the binary treatment status indicator in the \code{students} dataframe.
}
  \item{school.id}{
	the column name of the unique school ID in the \code{students} dataframe.
}
  \item{match.students}{
	logical value.  If \code{TRUE}, students are matched within school pairs and some students will be excluded.  If \code{FALSE}, all students will be retained in the matched sample for each school pair.
}
  \item{student.vars}{
 	column names of variables in \code{students} on which to match students and assess balance of student matches in evaluating match quality.
}
  \item{verbose}{
	a logical value indicating whether detailed output should be printed.
}
  \item{penalty.qtile}{
	a numeric value between 0 and 1 specifying a quantile of the distribution of all student-student matching distances.  The algorithm  will prefer to exclude treated students rather than form pairs with distances exceeding this quantile.
}
  \item{min.keep.pctg}{
	a minimum percentage of students in the smaller school in a pair which must be retained, even when treated students are excluded.  
}
}
\details{
	The \code{penalty.qtile} and \code{min.keep.pctg} control the rate at which students are trimmed from the match.  If the quantile is high enough no students should be excluded in any match; if the quantile is very low the \code{min.keep.pctg} can still ensure a minimal sample size in each match.
}
\value{
A list with one element for each treated school.  Each element is a list with one element for each control school, and each element of these secondary lists is a dataframe containing the matched sample for the corresponding treated-control pairing.
}

\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
