
# library(matchMulti)
# library( testthat )

context("calculating pvalues")

data(minischool)

p.tt <- matchMulti:::ttest.balance( "ses", treatment = "sector", orig.data = minischool )

# Handle clustering
p.CRVE <- matchMulti:::CRVE.balance( "ses", treatment = "sector", school.id = "school", data = minischool )

p.agg <- matchMulti:::agg.balance( "ses", treatment = "sector", school.id = "school", data = minischool )

expect_true( p.CRVE > p.tt )
expect_true( p.agg > p.tt )


