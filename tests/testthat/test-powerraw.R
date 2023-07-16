#context("Identical results for raw and SMD power functions")

#library("TOSTER")


test_that("powerTOSTtwo.prop is functions",{

  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }


  hush(suppressMessages({
    ## Sample size for alpha = 0.05, 90% power, assuming true effect prop1 = prop 2 = 0.5,
    ## equivalence bounds of 0.4 and 0.6 (so low_eqbound_prop = -0.1 and high_eqbound_prop = 0.1)

    powerTOSTtwo.prop(alpha = 0.05, statistical_power = 0.9, prop1 = 0.5, prop2 = 0.5,
                      low_eqbound_prop = -0.1, high_eqbound_prop = 0.1)

    ## Power for alpha = 0.05, N 542 , assuming true effect prop1 = prop 2 = 0.5,
    ## equivalence bounds of 0.4 and 0.6 (so low_eqbound_prop = -0.1 and high_eqbound_prop = 0.1)

    powerTOSTtwo.prop(alpha = 0.05, N = 542, prop1 = 0.5, prop2 = 0.5,
                      low_eqbound_prop = -0.1, high_eqbound_prop = 0.1)

    ## Equivalence bounds for alpha = 0.05, N 542 , assuming true effect prop1 = prop 2 = 0.5,
    ## and 90% power

    powerTOSTtwo.prop(alpha=0.05, statistical_power=0.9, N=542, prop1 = 0.5, prop2 = 0.5)

    #Example 4.2.4 from Chow, Wang, & Shao (2007, p. 93)
    powerTOSTtwo.prop(alpha=0.05, statistical_power=0.8, prop1 = 0.75, prop2 = 0.8,
                      low_eqbound_prop = -0.2, high_eqbound_prop = 0.2)

    # Example 5 from Julious & Campbell (2012, p. 2932)
    powerTOSTtwo.prop(alpha=0.025, statistical_power=0.9, prop1 = 0.8, prop2 = 0.8,
                      low_eqbound_prop=-0.1, high_eqbound_prop=0.1)
    # From Machin, D. (Ed.). (2008). Sample size tables for clinical studies (3rd ed).

    # Example 9.4b equivalence of two proportions (p. 113) #
    expect_equal(round(powerTOSTtwo.prop(alpha=0.010, statistical_power=0.8, prop1 = 0.5, prop2 = 0.5,
                      low_eqbound_prop = -0.2, high_eqbound_prop = 0.2)/2,0),
                 81)
  }

  ))

})
