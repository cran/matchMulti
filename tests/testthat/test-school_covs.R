library(matchMulti)
# library( testthat )

context("school.fb options")

data(catholic_schools)

# Trim data to speed up example
catholic_schools <- catholic_schools[catholic_schools$female_mean >.45 &
                                       catholic_schools$female_mean < .60,]
nrow( catholic_schools )
head( catholic_schools )

student.cov <- c('minority', 'female')


# Does school.fb work?
catholic_schools = dplyr::mutate( catholic_schools,
                                  size_cut = cut( size, 2 ),
                                  discrm_cut = cut( discrm, 2 ),
                                  acad_cut = cut( acad, 3 ) )

catholic_schools = dplyr::filter( catholic_schools, 
                                  acad > quantile( acad, 0.25 ) & acad < quantile( acad, 0.75 ) )
nrow( catholic_schools )
table(  catholic_schools$sector, catholic_schools$school  )

#table( catholic_schools$discrm_cut, catholic_schools$size_cut, catholic_schools$acad_cut )


match.simple2 <- matchMulti(catholic_schools, treatment = 'sector',
                            school.id = 'school', match.students = FALSE,
                            school.fb = list( c( "size_cut" ), c( "acad_cut" ) ) )



test_that( "simple list gets packaged right", {
  
  # These should be the same
  match.simpleA <- matchMulti(catholic_schools, treatment = 'sector',
                              school.id = 'school', match.students = FALSE,
                              school.fb = list( c( "discrm_cut", "acad_cut" ) ) )
  
  
  match.simpleB <- matchMulti(catholic_schools, treatment = 'sector',
                              school.id = 'school', match.students = FALSE,
                              school.fb =  c( "discrm_cut", "acad_cut" ) )
  
  expect_equal( match.simpleA$school.match, match.simpleB$school.match )
  
})
