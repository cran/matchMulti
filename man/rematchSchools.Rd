\name{rematchSchools}
\alias{rematchSchools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Repeat School Match Only
}
\description{
	After \code{matchMulti} has been called, repeats the school match (with possibly different parameters)  without repeating the more computationally intensive student match.
}
\usage{
rematchSchools(match.out, students, school.fb = NULL, verbose = FALSE,
keep.target = NULL, school.penalty = NULL, tol = 1e-3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{match.out}{
	an object returned by a call to \code{matchMulti}.
}
  \item{students}{
	a dataframe containing student and school covariates, with a different row for each student.
}
  \item{school.fb}{
	an optional list of character vectors, each containing a subset of the column names of \code{students}.  Each element of the list should contain all the names in previous elements (producing a nested structure).
}

  \item{verbose}{
	a logical value indicating whether detailed output should be printed.
}
  \item{keep.target}{
	an optional numeric value specifying the number of treated schools desired in the final match.  
}
  \item{school.penalty}{
	an optional numeric value, treated as the cost (to the objective function in the underlying optimization problem) of excluding a treated school.  If it is set lower, more schools will be excluded.
}

\item{tol}{
	a numeric tolerance value for comparing distances.  It may need to be raised above the default when matching with many levels of refined balance.
}
}
\details{
		The \code{school.fb} argument encodes a refined covariate balance constraint: the matching algorithm optimally balances the interaction of the variables in the first list element, then attempts to further balance the interaction in the second element, and so on.  As such variables should be added in order of priority for balance.
	
	The \code{keep.target} and \code{school.penalty} parameters allow optimal subset matching within the school match. When the \code{keep.target} argument is specified, the school match is repeated for different values of the \code{school.penalty} parameter in a form of binary search until an optimal match is obtained with the desired number of treated schools or a stopping rule is reached.  The \code{tol} parameter controls the stopping rule; smaller values provide a stronger guarantee of obtaining the exact number of treated schools desired but may lead to greater computational costs. 
	
	It is not recommended that users specify the \code{school.penalty} parameter directly in most cases.  Instead the \code{keep.target} parameter provides an easier way to consider excluding schools.  
}

\references{
  Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
  Springer-Verlag.

  Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
  Springer-Verlag.
  
  Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally Chosen Subset in Observational Studies."  Journal of Computational and Graphical Statistics, 21.1, 57-71.

  
}
\author{
  Luke Keele, Penn State University, \email{ljk20@psu.edu}

  Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}
}
\examples{

# Load Catholic school data
data(catholic_schools)

student.cov <- c('minority','female','ses')
school.cov <- c('minority_mean','female_mean', 'ses_mean', 'size', 'acad')

#Match schools but not students within schools
match.simple <- matchMulti(catholic_schools, treatment = 'sector',
school.id = 'school', match.students = FALSE)

#Check balance after matching - this checks both student and school balance
balanceMulti(match.simple, student.cov = student.cov, school.cov = school.cov)

#now rematch excluding 2 schools
match.trimmed <- rematchSchools(match.simple, catholic_schools, keep.target = 13)
match.trimmed$dropped$schools.t
}
\seealso{
	\code{\link{matchMulti}}.
}
