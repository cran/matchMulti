\name{balanceMulti}
\alias{balanceMulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs balance checking after multilevel matching.
}
\description{
This function checks balance after multilevel balance. It checks balance on both level-one (student) and level-two (school) covariates.
It reports results from both the t-test and nonparametric tests including Fisher's exact test and Wilcoxon signed rank test.
}
\usage{
balanceMulti(match.obj, student.cov = NULL, school.cov = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{match.obj}{
A multilevel match object
}
  \item{student.cov}{
Names of student level covariates that you want to check balance
}
  \item{school.cov}{
Names of school level covariates that you want to check balance. If not listed balance is checked
on all school level covariates included in the match.
}
}
\details{
This function returns a list which include balance checks for before and after matching for both level-one and level-two covariates.
Balance statistics include treated and control means, standardized differences, which is the difference in means divided by
the pooled standard deviation before matching, and p-values from parametric and nonparametric tests.
}
\value{
\item{students}{Balance table for student level covariates}
\item{schools}{Balance table for school level covariates}
}
%% \references{
%% ~put references to the literature/web site here ~
%% }
\author{  Luke Keele, Penn State University, \email{ljk20@psu.edu}    Sam Pimentel, University of Pennsylvania, \email{spi@wharton.upenn.edu}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	See also \code{\link{matchMulti}}, \code{\link{matchMultisens}}, 
	\code{\link{matchMultioutcome}}, \code{\link{rematchSchools}}
}
\examples{
# Load Catholic school data
data(catholic_schools)

student.cov <- c('minority','female','ses','mathach')

# Check balance student balance before matching
balanceTable(catholic_schools[c(student.cov,'sector')],  treatment = 'sector')

#Match schools but not students within schools
match.simple <- matchMulti(catholic_schools, treatment = 'sector', 
school.id = 'school', match.students = FALSE)

#Check balance after matching - this checks both student and school balance
balanceMulti(match.simple, student.cov = student.cov)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ balance }

