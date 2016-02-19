\name{catholic_schools}
\docType{data}
\alias{catholic_schools}
\title{1980 and 1982 High School and Beyond Data}

\description{
  These data are a subset of the data used in 
  Raudenbush and Bryk (1999) for multilevel modeling.
}
\usage{
	catholic_schools
	}

\format{
A \code{data.frame} with 1595 observations on the following variables.

school: unique school level identifier
ses:  student level socio-economic status scale ranges from approx. -3.578 to 2.692
mathach: senior year mathematics test score, outcome measure
female: student level indicator for sex
minority:  student level indicator for minority
minority_mean: school level measure of percentage of student body that is minority
female_mean: school level measure of percentage of student body that is female
ses_mean: school level measure of average level of student socio-economic status
sector: treatment indicator 1 if catholic 0 if public
size: school level measure of total number of enrolled students
acad: school level measure of the percentage of students on the academic track
discrm: school level measure of disciplinary climate ranges from approx. -2.4 to 2.7
size_large: school level indicator for schools with more than 1000 students
minority_mean_large: school level indicator for schools with more than ten percent minority
}
  
\source{
   Raudenbush, S. W. and Bryk, A. (2002).
   \emph{Hierarchical Linear Models: Applications and Data Analysis Methods}. 
   Thousand Oaks, CA: Sage.
}
\references{
  United States Department of Education. National Center for Education Statistics. 
  High School and Beyond, 1980: Sophomore and Senior Cohort First Follow-Up (1982).
}
\keyword{datasets}
