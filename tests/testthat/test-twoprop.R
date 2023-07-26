# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("test messages and warnings", {

  # two-tailed -----
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "d")
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "o")
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "r")
  # less -----
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "d",
                       alternative = "l")
  expect_equal(test1$conf.int[1],-Inf)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "o",
                       alternative = "l")
  expect_equal(test1$conf.int[1],-Inf)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "r",
                       alternative = "l")
  expect_equal(test1$conf.int[1],-Inf)

  # greater ------
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "d",
                       alternative = "g")
  expect_equal(test1$conf.int[2],Inf)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "o",
                       alternative = "g")
  expect_equal(test1$conf.int[2],Inf)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "r",
                       alternative = "g")
  expect_equal(test1$conf.int[2],Inf)

  # EQU -----
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "d",
                       alternative = "e",
                       null = .15)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "o",
                       alternative = "e",
                       null = 1.25)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "r",
                       alternative = "e",
                       null = 1.25)
  test_desc = describe_htest(test1)

  # MET -----
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "d",
                       alternative = "m",
                       null = .15)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "o",
                       alternative = "m",
                       null = 1.25)
  test1 = twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "r",
                       alternative = "m",
                       null = 1.25)

  # errors ----

  expect_message(twoprop_test(p1 = .5, p2 = .3,
                       n1 = 48, n2 = 40,
                       effect_size = "d"))
  expect_message(twoprop_test(p1 = .5, p2 = .3,
                              n1 = 48, n2 = 40,
                              effect_size = "o"))
  expect_message(twoprop_test(p1 = .5, p2 = .3,
                              n1 = 48, n2 = 40,
                              effect_size = "r"))

  expect_error(twoprop_test(p1 = .5, p2 = .3,
                            alpha = -.05,
                              n1 = 48, n2 = 40,
                              effect_size = "r"))
  expect_error(twoprop_test(p1 = .5, p2 = 1.3,
                            #alpha = -.05,
                            n1 = 48, n2 = 40,
                            effect_size = "r"))
  expect_error(twoprop_test(p1 = .5, p2 = 3/9,
                              n1 = 48, n2 = 9,
                              effect_size = "r"))

  expect_error(twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "d",
                       alternative = "e"))
  expect_error(twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "o",
                       alternative = "e"))
  expect_error(twoprop_test(p1 = .5, p2 = .3,
                       n1 = 100, n2 = 100,
                       effect_size = "r",
                       alternative = "e"))


})

test_that("Random tests against prop_test",{

  smokers  <- c( 83, 90, 129, 70 )
  patients <- c( 86, 93, 136, 82 )

  smokers  <- c(83, 90)
  patients <- c(86, 93)
  ptest_base = suppressWarnings( prop.test(smokers, patients,
                         correct = FALSE))
  ptest_prop = twoprop_test(smokers[1]/patients[1],
                            smokers[2]/patients[2],
                            patients[1],
                            patients[2])

  expect_equal(abs(ptest_base$p.value - ptest_prop$p.value),0,tolerance=.001)

  expect_equal(abs(ptest_base$conf.int[1] - ptest_base$conf.int[1]),0,tolerance=.001)
  expect_equal(abs(ptest_base$conf.int[2] - ptest_base$conf.int[2]),0,tolerance=.001)

  smokers  <- c(129, 70)
  patients <- c(136, 82)
  ptest_base = suppressWarnings( prop.test(smokers, patients,
                                           correct = FALSE))
  ptest_prop = twoprop_test(smokers[1]/patients[1],
                            smokers[2]/patients[2],
                            patients[1],
                            patients[2])

  expect_equal(abs(ptest_base$conf.int[1] - ptest_base$conf.int[1]),0,tolerance=.001)
  expect_equal(abs(ptest_base$conf.int[2] - ptest_base$conf.int[2]),0,tolerance=.001)
  set.seed(16281940)
  for(i in 1:100){
    #print(i)
    p1 = runif(1,.05,.95)
    p2 = runif(1,.05,.95)

    n1 = round(runif(1,55,1000),0)
    n2 = n1

    s1 = round(p1*n1,0)
    s2 = round(p2*n2,0)

    p1 = s1/n1
    p2 = s2/n2

    ptest_base = prop.test(c(s1,s2),
                           c(n1,n2),
                           correct = FALSE)

    ptest_prop = twoprop_test(p1,p2,n1,n2)
    #abs(ptest_base$p.value - ptest_prop$p.value)
    #if(ptest_base$p.value < 0.2){
    #  expect_equal(abs(round(ptest_base$p.value,2) - round(ptest_prop$p.value,2)),0,tolerance=.011)
    #}
    expect_equal(abs(ptest_base$conf.int[1] - ptest_base$conf.int[1]),0,tolerance=.001)
    expect_equal(abs(ptest_base$conf.int[2] - ptest_base$conf.int[2]),0,tolerance=.001)

    ptest_odds = twoprop_test(p1,p2,n1,n2,
                          effect_size = "o")

    ptest_risk = twoprop_test(p1,p2,n1,n2,
                              effect_size = "r")

    if(ptest_prop$p.value < .01){
      expect_equal(abs(ptest_prop$p.value - ptest_odds$p.value),0,tolerance=.01)
      expect_equal(abs(ptest_risk$p.value - ptest_odds$p.value),0,tolerance=.01)
    } else if(ptest_prop$p.value < .1){
      expect_equal(abs(ptest_prop$p.value - ptest_odds$p.value),0,tolerance=.02)
      expect_equal(abs(ptest_risk$p.value - ptest_odds$p.value),0,tolerance=.02)
    } else {
      expect_equal(abs(ptest_prop$p.value - ptest_odds$p.value),0,tolerance=.1)
      expect_equal(abs(ptest_risk$p.value - ptest_odds$p.value),0,tolerance=.1)
    }


  }
})
