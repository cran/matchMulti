
# library(matchMulti)
# library( testthat )

context("balance tables")

data(minischool)



student.cov <- c('mathach', 'minority', 'female')

# Does school.fb work?
minischool = dplyr::mutate( minischool,
                            size_cut = cut( size, 2 ),
                            discrm_cut = cut( discrm, 2 ) )

nrow( minischool )
table(  minischool$sector, minischool$school  )


school.fb <- c( "discrm_cut", "size_cut" )
school.cov = c( "discrm", "size" )


##### Check balance on raw data ######

bt <- balanceTable( minischool, 
                    treatment = "sector",
                    school.id = "school",
                    var.names = c( student.cov, school.cov ),
                    include.tests = TRUE )
bt



btNt <- balanceTable( minischool, 
                      treatment = "sector",
                      var.names = c( student.cov, school.cov ),
                      include.tests = FALSE )
btNt

expect_true( nrow( btNt ) == 5 )
expect_true( all( btNt$`Treated Mean` == bt$`Treated Mean` ) )
expect_true( all( btNt$`Control Mean` == bt$`Control Mean` ) )
expect_true( all( btNt$SDiff == bt$SDiff ) )



# Are we getting balance on categorical levels?
btN <- balanceTable( minischool, 
                     treatment = "sector",
                     include.tests = FALSE )
btN
# We should see multiple discrm_cut and size_cuts 

expect_true( nrow(btN) > ncol(minischool) )



fakedbl = balanceTable( minischool, df.match=minischool, 
                        treatment = "sector", school.id = "school", 
                        var.names = c( "minority", "discrm_cut" ),
                        include.tests = TRUE )
fakedbl
expect_true( all( fakedbl$`After Agg PValue` == fakedbl$`Before Agg PValue` ))



##### Now do matching and check balance ######

if (requireNamespace("optmatch", quietly = TRUE)){
  
  match.simpleA <- matchMulti(minischool, treatment = 'sector',
                              school.id = 'school', match.students = FALSE,
                              school.fb = school.fb )
  match.simpleA
  
  
  
  test_that( "balanceMulti works", {
    
    btab_split = balanceMulti( match.simpleA,
                               school.cov = school.cov, 
                               student.cov = student.cov )
    btab_split
    
    expect_true( nrow( btab_split$schools ) == 2 )
    expect_true( nrow( btab_split$students ) == 3 )
  } )
  
  
  
  test_that( "single.table works", {
    
    btab = balanceMulti( match.simpleA, single.table = TRUE, include.tests = TRUE )
    btab
    
    btab2 = balanceMulti( match.simpleA, single.table = TRUE, include.tests = FALSE )
    btab2
    
    expect_equal( nrow(btab), nrow(btab2) )
    expect_equal( dim(btab2), c( 14, 6 ) )
    
  } )
  
}



#### Testing weights for balance tables #####



# test_that( "weights work (Basic tests)", {
#   
#   ns = table( minischool$sector )
#   ns
#   
#   bt <- expect_error( balanceTable( minischool, 
#                                     treatment = "sector",
#                                     school.id = "school",
#                                     var.names = c( student.cov, school.cov ),
#                                     treat.wts = runif( 1,5, ns[[2]] ),
#                                     include.tests = TRUE )
#   )
#   
#   bt <- balanceTable( minischool, 
#                       treatment = "sector",
#                       school.id = "school",
#                       var.names = c( student.cov, school.cov ),
#                       treat.wts = runif( ns[[2]], 1,5 ),
#                       ctrl.wts = runif( ns[[1]], 1,5 ),
#                       include.tests = TRUE )
#   
#   bt
#   expect_true( all( !is.na( bt$`Agg PValue` ) ) )
#   expect_true( all( is.na( bt$`CRVE PValue` ) ) )
#   
#   bt <- balanceTable( minischool, 
#                       treatment = "sector",
#                       school.id = "school",
#                       var.names = c( student.cov, school.cov ),
#                       treat.wts = rep( 1, ns[[2]] ),
#                       ctrl.wts = rep( 2, ns[[1]] ),
#                       include.tests = TRUE )
#   
#   bt
#   expect_true( all( !is.na( bt$`Agg PValue` ) ) )
#   expect_true( all( !is.na( bt$`CRVE PValue` ) ) )
#   
# } )




