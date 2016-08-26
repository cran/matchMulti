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