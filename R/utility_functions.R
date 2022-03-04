
# Check: Are variable names present in data frame?
validate_inputs = function( students, treatment, school.id, 
                            student.vars = NULL, school.vars = NULL ) {
  
  if (!(treatment %in% colnames(students))) {
    stop(paste('Treatment variable', treatment,'not found'))
  }
  if (!(school.id %in% colnames(students))) {
    stop(paste('School ID variable', school.id,'not found'))
  }	
  if (!is.null(student.vars) && !all(student.vars %in% colnames(students))) {
    bad.vars <- student.vars[!(student.vars %in% colnames(students))]
    stop(paste('Student variable', bad.vars,'not found\n'))
  }
  school.vars = unique( unlist(school.vars) )
  if (!is.null(school.vars) && !all( school.vars %in% colnames(students))) {
    bad.vars <- school.vars[!(school.vars %in% colnames(students))]
    stop(paste('School variable', bad.vars,'not found\n'))
  }	
  
}






scat = function( str, ... ) {
  cat( sprintf( str, ... ) )
}



sstop <- function( str, ... ) {
  stop( sprintf( str, ... ), call. = FALSE )
}


smessage <- function( str, ... ) {
  message( sprintf( str, ... ) )
}



swarning <- function( str, ... ) {
  warning( sprintf( str, ... ), call. = FALSE )
}

