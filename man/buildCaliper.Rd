\name{buildCaliper}
\alias{buildCaliper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct propensity score caliper
}
\description{
Fits a propensity score for an individual-level or group-level treatment, computes a caliper for the propensity score (based on a fractional number of standard deviations provided by the user), and creates a matrix containing information about which treated-control pairings are excluded by the caliper.
}
\usage{
buildCaliper(data, treatment, ps.vars, group.id = NULL, caliper = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the treatment variable, the variables to be used in fitting the propensity score and (if treatment is at the group level) a group ID.
}
  \item{treatment}{
Name of the treatment indicator.
}
  \item{ps.vars}{
Vector of names of variables to use in fitting the propensity score.  
}

  \item{group.id}{
Name of group ID variable, if applicable.  
}

\item{caliper}{
	Desired size of caliper, in number of standard deviations of the fitted propensity score.
}

}
\details{
The \code{treatment} variable should be binary with 1 indicating treated units and 0 indicating controls.  When \code{group.id} is \code{NULL}, treatment is assumed to be at the individual level and the propensity score is fitted using the matrix \code{data}.  When a group ID is specified, data frame \code{data} is first aggregated into groups, with variables in \code{ps.vars} replaced by their within-group means, and the propensity score is fitted on the group matrix.
}
\value{
A matrix with \code{nrow} equal to the number of treated individuals or groups and \code{ncol} equal to the number of control individuals, with \code{0} entries indicating pairings permitted by the caliper and \code{Inf} entries indicating forbidden pairings.
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{  Luke Keele, Penn State University, \email{ljk20@psu.edu}    Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
	\dontrun{
# Load Catholic school data
data(catholic_schools)

student.cov <- c('minority','female','ses','mathach')

# Check balance student balance before matching
balanceTable(catholic_schools[c(student.cov,'sector')],  treatment = 'sector')

#fit a propensity score caliper on mean values of student covariates within schools
school.caliper <- buildCaliper(data = catholic_schools, treatment = 'sector',
	ps.vars = student.cov, group.id = 'school')

#Match schools but not students within schools
match.simple <- matchMulti(catholic_schools, treatment = 'sector', 
	school.caliper = school.caliper, school.id = 'school', match.students = FALSE)

#Check balance after matching - this checks both student and school balance
balanceMulti(match.simple, student.cov = student.cov)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ balance }

