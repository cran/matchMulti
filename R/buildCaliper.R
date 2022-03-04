#' Construct propensity score caliper
#' 
#' Fits a propensity score for an individual-level or group-level treatment,
#' computes a caliper for the propensity score (based on a fractional number of
#' standard deviations provided by the user), and creates a matrix containing
#' information about which treated-control pairings are excluded by the
#' caliper.
#' 
#' The \code{treatment} variable should be binary with 1 indicating treated
#' units and 0 indicating controls.  When \code{group.id} is \code{NULL},
#' treatment is assumed to be at the individual level and the propensity score
#' is fitted using the matrix \code{data}.  When a group ID is specified, data
#' frame \code{data} is first aggregated into groups, with variables in
#' \code{ps.vars} replaced by their within-group means, and the propensity
#' score is fitted on the group matrix.
#' 
#' @param data A data frame containing the treatment variable, the variables to
#' be used in fitting the propensity score and (if treatment is at the group
#' level) a group ID.
#' @param treatment Name of the treatment indicator.
#' @param ps.vars Vector of names of variables to use in fitting the propensity
#' score.
#' @param group.id Name of group ID variable, if applicable.
#' @param caliper Desired size of caliper, in number of standard deviations of
#' the fitted propensity score.
#' @return A matrix with \code{nrow} equal to the number of treated individuals
#' or groups and \code{ncol} equal to the number of control individuals, with
#' \code{0} entries indicating pairings permitted by the caliper and \code{Inf}
#' entries indicating forbidden pairings.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords balance
#' @examples
#' 
#' 	\dontrun{
#' # Load Catholic school data
#' data(catholic_schools)
#' 
#' student.cov <- c('minority','female','ses','mathach')
#' 
#' # Check balance student balance before matching
#' balanceTable(catholic_schools[c(student.cov,'sector')],  treatment = 'sector')
#' 
#' #fit a propensity score caliper on mean values of student covariates within schools
#' school.caliper <- buildCaliper(data = catholic_schools, treatment = 'sector',
#' 	ps.vars = student.cov, group.id = 'school')
#' 
#' #Match schools but not students within schools
#' match.simple <- matchMulti(catholic_schools, treatment = 'sector', 
#' 	school.caliper = school.caliper, school.id = 'school', match.students = FALSE)
#' 
#' #Check balance after matching - this checks both student and school balance
#' balanceMulti(match.simple, student.cov = student.cov)
#' }
#' 
#' @export buildCaliper
buildCaliper <- function(data, treatment, ps.vars, group.id = NULL, caliper = 0.2){
	#if we are computing a group-level propensity score, aggregate variables by group
	if(!is.null(group.id)){
		school.df <- students2schools(data, c(ps.vars, treatment), group.id)
	} else {
		school.df <- data
	}
	pscores <- predict(glm(as.formula(paste(treatment, '~', paste(ps.vars, collapse = ' + '))), data = school.df, family = binomial()))
	ps.diffs <- abs(outer(pscores[school.df[[treatment]] == 1], pscores[school.df[[treatment]] == 0], FUN = '-'))
	
	calip <- sd(pscores)*caliper
	ps.diffs[ps.diffs > calip] <- Inf
	ps.diffs[is.finite(ps.diffs)] <- 0
	ps.diffs
}
