



#' @title matchMultiResult object for results of power calculations
#' @name matchMultiResult
#'
#' @description
#' The matchMultiResult object is an S3 class that holds the results from the matchMulti call.
#'
#' matchMulti result objects have the matched datasets inside of them.
#'
#' @param x a matchMultiResult object (except for is.matchMultiResult, where it is a generic
#'   object to check).
#' @rdname matchMultiResult
NULL





#' @return is.matchMultiResult: TRUE if object is a matchMultiResult object.
#'
#' @export
#'
#' @rdname matchMultiResult
is.matchMultiResult = function( x ) {
  inherits(x, "matchMultiResult")
}



format_fine_balance = function( school.fb ) {
  pieces = sapply( school.fb, paste0, collapse = ", " )
  
  paste0( pieces, collapse = " / " )
  
}



#' Pretty print matchMulti result
#'
#' @export
#' @param ... No extra options passed.
#' @rdname matchMultiResult
print.matchMultiResult = function( x, ... ) {
  cat( "Multilevel Match Result\n" )
  school.id = x$school.id
  treatment = x$treatment
  labs = sort( unique( x$raw[[treatment]] ) )
  
  scat( "School ID: '%s'.\nTreatment: '%s' with Co = %s, Tx = %s\n", 
        school.id, treatment, labs[[1]], labs[[2]] )
  
  if ( !is.null( x$student.vars ) ) {
    scat( "Student vars: %s\n", paste0( x$student.vars, collapse = ", " ) )
  }
  if ( !is.null( x$school.fb ) ) {
    scat( "School fine balance: %s\n", format_fine_balance( x$school.fb ) )
  }
  
  scat( "Student sample size: From %d -> %d students\n", nrow( x$raw ), nrow( x$matched ) )
  
  n_raw = tally_schools(x$raw, school.id, treatment )
  n_match = tally_schools(x$matched, school.id, treatment )
  scat( "Final match: %d Co and %d Tx schools (from %d and %d schools)\n", 
        n_match[[1]], n_match[[2]],
        n_raw[[1]], n_raw[[2]] )
  
  invisible( x )
}



#' Pretty print match result
#'
#' @export
#' @param object Object to summarize.
#' @param ... Extra options passed to print.matchMultiResult
#' @rdname matchMultiResult
summary.matchMultiResult = function( object, ... ) {
  cat( "Multilevel Match\n" )

  school.id = object$school.id
  treatment = object$treatment
  
  labs = sort( unique( object$raw[[treatment]] ) )
  
  scat( "\tSchool ID: '%s'.\n\tTreatment: '%s' with Co = %s, Tx = %s\n", 
        school.id, treatment, labs[[1]], labs[[2]] )
  
  if ( !is.null( object$student.vars ) ) {
    scat( "\tMatched on: %s\n", paste0( object$student.vars, collapse = ", " ) )
  }
  if ( !is.null( object$school.fb ) ) {
    scat( "\tFine balance on: %s\n", format_fine_balance(object$school.fb) )
  }
  
  scat( "\tFrom %d -> %d students\n", nrow( object$raw ), nrow( object$matched ) )
  
  scat( "\t%d control and %d treatment schools dropped\n",
        length( object$dropped$schools.c ), length( object$dropped$schools.t ) )
  
  cat( "\nFinal matched data:\n" )
  describe_data_counts( object$matched, 
                        treatment = treatment,
                        school.id = school.id )
  
  
  invisible( object )
}


#' Tally schools and students in a given dataset
#'
#'  Returns a count of schools, without printing
#' anything.
#'
#' @param data Dataset (student level)
#' @param school.id String name of ID column in data (the grouping variable)
#' @param treatment String name of the treatment variable.
#'
#' @seealso describe_data_counts
#' 
#' @author Luke Miratrix
#'
#' @return List of two things: school and student counts (invisible).
#' 
#' @export
tally_schools = function( data, school.id, treatment ) {
  labs = sort( unique( data[[treatment]] ) )
  
  mid = data[[school.id]]
  mtx = data[[treatment]]
  nco = length( unique( mid[mtx==labs[[1]] ]) )
  ntx = NA
  if ( length( labs ) > 1 ) {
    ntx = length( unique( mid[mtx==labs[[2]] ] ) )
  }
  return( c( Co = nco, Tx = ntx, N = nco + ntx ) )
}


#' Print out summary of student and school counts
#' 
#' Given a school ID and treatment variable, count up number of schools and
#' students, print out a summary of the counts of students and
#' schools.
#' 
#' @seealso tally_schools
#' 
#' @inheritParams  tally_schools
#' @return List of three numbers, # control, # Tx, # Total
#' 
#' @export
describe_data_counts = function( data, school.id, treatment ) {
  
  stopifnot( treatment %in% names(data) )
  stopifnot( school.id %in% names(data) )
  
  #do schools nest within treatment categories?
  treat.tab <- table(data[[school.id]], data[[treatment]]) 
  if(any(apply(treat.tab, 1, min) > 0)) {
    warning('Some schools contain both treated and control students')
  }
  
  labs = colnames( treat.tab )
  
  if ( dim( treat.tab )[[2]] == 1 ) {
    stop("Only have single level of treatment" )
  }
  if ( dim( treat.tab )[[2]] > 2 ) {
    stop( sprintf( "More than 2 levels of treatment (%d levels)", dim(treat.tab)[[2]] ) )
  }
  
  nsch = apply( treat.tab > 0, 2, sum )
  nsch  
  
  tx = treat.tab[ , 1 ] == 0

  
  scat( "Schools:\n\t# Control (%s): %d\n\t# Treatment (%s): %d\n",
        labs[[1]], nsch[[1]], labs[[2]], nsch[[2]] )
  
  studtot = apply( treat.tab, 2, sum )
  scat( "\nStudents:\n\t# Control (%s): %d\n\t\tAvg %.1f / school\n\t\tRange %d - %d",
        labs[[1]], studtot[[1]], studtot[[1]]/nsch[[1]],
        min(treat.tab[!tx,1]), max(treat.tab[!tx,1] ))
  
  scat( "\n\t# Treatment (%s): %d\n\t\tAvg %.1f / school\n\t\tRange %d - %d", 
        labs[[2]], studtot[[2]], studtot[[2]]/nsch[[2]], 
        min(treat.tab[tx,2]), max(treat.tab[tx,2] ))
  
  
  invisible( list( schools = nsch,
                   students = studtot ) )
}
