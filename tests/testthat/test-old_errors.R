
test_that("Errors for TOSTtwo functions",{
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  hush(suppressMessages( {

    withr::local_options(lifecycle_verbosity = "quiet")
  expect_warning(TOSTtwo.raw(m1 = 5.25, m2 = 5.22,
                     sd1 = 0.95, sd2 = 0.83,
                     n1 = 95, n2 = 89,
                     low_eqbound = 0.385, high_eqbound = 0.384, plot = FALSE))
  expect_warning(TOSTtwo(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,
                                low_eqbound_d=0.44,high_eqbound_d=0.43,
                         verbose = FALSE,
                         plot = FALSE))

  expect_error(TOSTtwo(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=1,n2=89,
                         low_eqbound_d=-0.44,high_eqbound_d=0.43,
                         verbose = FALSE,
                         plot = FALSE))
  expect_error(TOSTtwo.raw(m1 = 5.25, m2 = 5.22,
                             sd1 = 0.95, sd2 = 0.83,
                             n1 = 1, n2 = 89,
                             low_eqbound = -0.385, high_eqbound = 0.384, plot = FALSE))

  expect_error(TOSTtwo.raw(m1 = 5.25, m2 = 5.22,
                           sd1 = 0, sd2 = 0.83,
                           n1 = 95, n2 = 89,
                           low_eqbound = -0.385, high_eqbound = 0.384, plot = FALSE))

  expect_error(TOSTtwo.raw(m1 = 5.25, m2 = 5.22,
                           sd1 = 0.95, sd2 = 0.83,
                           alpha = 1,
                           n1 = 95, n2 = 89,
                           low_eqbound = -0.385, high_eqbound = 0.384, plot = FALSE))

  expect_error(TOSTtwo(m1=5.25,m2=5.22,sd1=0.95,sd2=0.83,n1=95,n2=89,
                       low_eqbound_d=-0.44,high_eqbound_d=0.43,
                       alpha = 1,
                       verbose = FALSE,
                       plot = FALSE))
  }))

  })


test_that("Errors for TOSTpaired functions",{
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  hush(suppressMessages( {
  expect_warning(TOSTpaired(n=65,m1=5.83,m2=5.75,
                            sd1=1.17,sd2=1.29,r12=0.75,
                            low_eqbound_dz=0.4,high_eqbound_dz=0.4,
                            verbose = FALSE,
                            plot = FALSE))

  expect_warning(TOSTpaired.raw(n=65,m1=5.83,m2=5.75,
                                sd1=1.17,sd2=1.30,r12=0.745,
                                low_eqbound=0.34,high_eqbound=0.34,
                                verbose = FALSE,
                                plot = FALSE))

  expect_error(TOSTpaired(n=1,m1=5.83,m2=5.75,
                          sd1=1.17,sd2=1.29,r12=0.75,
                          low_eqbound_dz=-0.4,high_eqbound_dz=0.4,
                          verbose = FALSE,
                          plot = FALSE))
  expect_error(TOSTpaired.raw(n=65,m1=5.83,m2=5.75,
                              sd1=0,sd2=1.30,r12=0.745,
                              low_eqbound=-0.34,high_eqbound=0.34,
                              verbose = FALSE,
                              plot = FALSE))

  expect_error(TOSTpaired.raw(n=65,m1=5.83,m2=5.75,
                              sd1=1.17,sd2=1.30,r12=1.5,
                              low_eqbound=-0.34,high_eqbound=0.34,
                              verbose = FALSE,
                              plot = FALSE))

  expect_error(TOSTpaired.raw(n=1,m1=5.83,m2=5.75,
                              sd1=1.17,sd2=1.30,r12=0.745,
                              low_eqbound=-0.34,high_eqbound=0.34,
                              verbose = FALSE,
                              plot = FALSE))

  expect_error(TOSTpaired(n=65,m1=5.83,m2=5.75,
                          sd1=0,sd2=1.29,r12=0.75,
                          low_eqbound_dz=-0.4,high_eqbound_dz=0.4,
                          verbose = FALSE,
                          plot = FALSE))

  expect_error(TOSTpaired(n=65,m1=5.83,m2=5.75,
                          sd1=1.17,sd2=1.29,r12=-1.1,
                          low_eqbound_dz=-0.4,high_eqbound_dz=0.4,
                          verbose = FALSE,
                          plot = FALSE))
}))
})


test_that("Errors for TOSTmeta",{
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  hush(suppressMessages( {
  res = TOSTmeta(ES=0.12, se=0.09,
                 low_eqbound_d=-0.2, high_eqbound_d=0.2,
                 verbose = FALSE,
                 plot = FALSE)

  res = TOSTmeta(ES=0.12, var=0.09,
                 low_eqbound_d=-0.2, high_eqbound_d=0.2,
                 verbose = FALSE,
                 alpha=0.05, plot = FALSE)

  expect_error(TOSTmeta(ES=0.12,
                        low_eqbound_d=-0.2, high_eqbound_d=0.2,
                        verbose = FALSE,
                        alpha=0.05, plot = FALSE))

  expect_warning(TOSTmeta(ES=0.12,
                        se = .09,
                        low_eqbound_d=0.25, high_eqbound_d=0.2,
                        verbose = FALSE,
                        alpha=0.05, plot = FALSE))

  expect_error(TOSTmeta(ES=0.12,
                          se = .09,
                          low_eqbound_d=-0.25, high_eqbound_d=0.2,
                          verbose = FALSE,
                          alpha=1.1, plot = FALSE))
}))
})

test_that("Errors for TOSTone",{
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  hush(suppressMessages( {
  res <- TOSTone(m=0.54,mu=0.5,
                 sd=1.2,n=100,
                 low_eqbound_d=-0.3, high_eqbound_d=0.3,
                 alpha=0.05,
                 verbose = FALSE,
                 plot = FALSE)
  expect_warning(TOSTone(m=0.54,mu=0.5,
                       sd=1.2,n=100,
                       low_eqbound_d=0.31, high_eqbound_d=0.3,
                       alpha=0.05,
                       verbose = FALSE,
                       plot = FALSE))
  expect_error(TOSTone(m=0.54,mu=0.5,
                       sd=0,n=100,
                       low_eqbound_d=-0.3, high_eqbound_d=0.3,
                       alpha=0.05,
                       verbose = FALSE,
                       plot = FALSE))
  expect_error(TOSTone(m=0.54,mu=0.5,
                       sd=1.2,n=1,
                       low_eqbound_d=-0.3, high_eqbound_d=0.3,
                       alpha=0.05,
                       verbose = FALSE,
                       plot = FALSE))

  res <- TOSTone.raw(m=0.52,mu=0.5,
                     sd=0.5,n=300,
                     low_eqbound=-0.1, high_eqbound=0.1,
                     alpha=0.05, verbose=FALSE,
                     plot = FALSE)

  expect_error(TOSTone.raw(m=0.52,mu=0.5,
                           sd=0.5,n=300,
                           low_eqbound=-0.1, high_eqbound=0.1,
                           alpha=1.5, verbose=FALSE,
                           plot = FALSE))

  expect_error(TOSTone.raw(m=0.52,mu=0.5,
                           sd=0,n=300,
                           low_eqbound=-0.1, high_eqbound=0.1,
                           alpha=0.05, verbose=FALSE,
                           plot = FALSE))
  expect_error(TOSTone.raw(m=0.52,mu=0.5,
                           sd=0.5,n=1,
                           low_eqbound=-0.1, high_eqbound=0.1,
                           alpha=0.05, verbose=FALSE,
                           plot = FALSE))

  expect_warning(TOSTone.raw(m=0.52,mu=0.5,
                             sd=0.5,n=300,
                             low_eqbound=0.11, high_eqbound=0.1,
                             alpha=0.05, verbose=FALSE,
                             plot = FALSE))

  }))

})

test_that("Errors for TOSTr",{
  hush = function(code) {
    sink("NUL") # use /dev/null in UNIX
    tmp = code
    sink()
    return(tmp)
  }

  hush(suppressMessages( {
  res <- TOSTr(n=100, r = 0.02,
               low_eqbound_r=-0.3,
               high_eqbound_r=0.3,
               alpha=0.05,
               verbose = FALSE,
               plot = FALSE)
  expect_error(TOSTr(n=100, r = -1.1,
                     low_eqbound_r=-0.3,
                     high_eqbound_r=0.3,
                     alpha=0.05,
                     verbose = FALSE,
                     plot = FALSE))
  expect_error(TOSTr(n=100, r = 0.02,
                     low_eqbound_r=-0.3,
                     high_eqbound_r=0.3,
                     alpha=1.05,
                     verbose = FALSE,
                     plot = FALSE))
  expect_error(TOSTr(n=1, r = 0.02,
                     low_eqbound_r=-0.3,
                     high_eqbound_r=0.3,
                     alpha=0.05,
                     verbose = FALSE,
                     plot = FALSE))
  expect_warning(TOSTr(n=100, r = 0.02,
        low_eqbound_r=0.32,
        high_eqbound_r=0.3,
        alpha=0.05,
        verbose = FALSE,
        plot = FALSE))

  }))

})
