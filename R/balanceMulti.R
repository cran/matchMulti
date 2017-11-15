balanceMulti <-
function(match.obj, student.cov = NULL, school.cov = NULL){
	treatment <- match.obj$treatment
	
	#if no student covariates are provided, compute balance on all
	if (is.null(student.cov)) student.cov = setdiff(colnames(match.obj$raw), c(treatment, match.obj$school.id))
	
	#student balance
	student.bal <- balanceTable(match.obj$raw[c(student.cov,treatment)], match.obj$matched[c(student.cov, treatment)], treatment)
	
	#school balance
	school.bal <- NULL
	if(!is.null(school.cov)){
		schools.raw <- students2schools(match.obj$raw, c(school.cov, treatment), match.obj$school.id)
		schools.matched <- schools.raw[which(schools.raw[[match.obj$school.id]] %in% as.vector(match.obj$school.match)),]
		
		id.idx <- which(colnames(schools.raw) == match.obj$school.id)
	
		school.bal <- balanceTable(schools.raw[-id.idx], schools.matched[-id.idx], treatment)
	}
	return(list('students' = student.bal, 'schools' = school.bal))
}
