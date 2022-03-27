# library(matchMulti)
# library(Hmisc)
# library( testthat )

context("general testing")


data(catholic_schools)
# Trim data to speed up example
catholic_schools <- catholic_schools[catholic_schools$female_mean >.45 &
 catholic_schools$female_mean < .60,]






#match on a single covariate 
student.cov <- c('minority')


match.simple <- matchMulti(catholic_schools, treatment = 'sector', 
                           school.id = 'school', match.students = FALSE, 
                           student.vars = student.cov, verbose=TRUE, tol=.01)

#Check balance after matching - this checks both student and school balance
bM <- balanceMulti(match.simple, student.cov = student.cov)

expect_true( is.list( bM ) )
expect_true( !is.null( bM$students ) )



#match on several covariates
student.cov <- c('minority','female','ses')

# Check balance student balance before matching
bT <- balanceTable(catholic_schools[c(student.cov,'sector')],  treatment = 'sector')
bT
expect_equal( dim(bT), c(3,3) )


#Match schools but not students within schools
match.simple <- matchMulti(catholic_schools, treatment = 'sector',
                           school.id = 'school', match.students = FALSE)




# Does school.fb work?
catholic_schools = dplyr::mutate( catholic_schools,
                                  size_cut = cut( size, 2 ),
                                  discrm_cut = cut( discrm, 2 ),
                                  acad_cut = cut( acad, 3 ) )

#table( catholic_schools$discrm_cut, catholic_schools$size_cut, catholic_schools$acad_cut )

match.simple2 <- matchMulti(catholic_schools, treatment = 'sector',
                            school.id = 'school', match.students = FALSE,
                            school.fb = list( c( "size_cut" ), c( "discrm_cut", "acad_cut" ) ) )






data <- catholic_schools
treatment = 'sector'
school.id = 'school'
match.students = FALSE
student.vars = NULL
school.caliper = NULL
school.fb = NULL
verbose = FALSE
keep.target = NULL
student.penalty.qtile = 0.05
min.keep.pctg = 0.8
school.penalty = NULL
save.first.stage = TRUE
tol = 1e-3


### Match on component pieces ####

student.matches <- matchStudents(data, treatment, school.id, match.students, 
                                 student.vars,  school.caliper, verbose, student.penalty.qtile, 
                                 min.keep.pctg)

school.match <- matchSchools(student.matches$schools.matrix, data, treatment, school.id, 
                             school.fb, school.penalty, verbose, tol = tol) 


dmat <- student.matches$schools.matrix
students <- catholic_schools
treatment = 'sector'
school.id = 'school'
school.fb = NULL
verbose = FALSE 
tol	<- 1e-3
penalty <- NULL

out.match <- assembleMatch(student.matches$student.matches, school.match, school.id, treatment)

expect_true( !is.null( out.match ) )
expect_true( is.data.frame( out.match ) )





#### Check sensitivity ####

head( match.simple2$matched )

sens <- matchMultisens(match.simple2, out.name = "mathach", schl_id_name = "school", treat.name = "sector", Gamma=1.5)

expect_true( sens$pval > 0 )



