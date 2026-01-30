#context("Run Examples for boot_t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("glass", {

  # library(TOSTER)
  # library(metafor)

  smd = TOSTER::smd_calc(x = subset(sleep,
                                     group == 1)$extra,
                          y = subset(sleep,
                                     group == 2)$extra,
                          glass = "glass2",
                          output = "data.frame")
  smc = TOSTER::smd_calc(x = subset(sleep,
                                     group == 1)$extra,
                          y = subset(sleep,
                                     group == 2)$extra,
                          glass = "glass1",
                          paired = TRUE,
                          output = "data.frame")
  # x1 = subset(sleep,
  #             group == 1)$extra
  # y1 = subset(sleep,
  #             group == 2)$extra
  # m1 = mean(x1)
  # sd1 = sd(x1)
  # n1 = length(x1)
  # m2 = mean(y1)
  # sd2 = sd(y1)
  # n2 = length(y1)
  # r12 = cor(x1,y1)
  # df1 = data.frame(
  #   m1,
  #   sd1,
  #   n1,
  #   m2,
  #   sd2,
  #   n2,
  #   ri = r12
  # )

  # meta_d = escalc(measure = "SMD1H",
  #                 data = df1,
  #                 m1i = m1,
  #                 sd1i = sd1,
  #                 n1i = n1,
  #                 m2i = m2,
  #                 sd2i = sd2,
  #                 n2i = n2)
  # meta_c = escalc(measure = "SMCRH",
  #                 data = df1,
  #                 m1i = m1,
  #                 sd1i = sd1,
  #                 ni = n1,
  #                 m2i = m2,
  #                 sd2i = sd2,
  #                 n2i = n2,
  #                 ri = ri)

  expect_equal(0.2287079,round(smd$SE^2,7))

  expect_equal(0.08871185,round(smc$SE^2,8))


})

test_that("Hedges g(s/av) and g(z)",
          {
            smd = TOSTER::smd_calc(x = subset(sleep,
                                               group == 1)$extra,
                                    y = subset(sleep,
                                               group == 2)$extra,
                                    var.equal = TRUE,
                                    output = "data.frame")
            smc = TOSTER::smd_calc(x = subset(sleep,
                                               group == 1)$extra,
                                    y = subset(sleep,
                                               group == 2)$extra,

                                    paired = TRUE,
                                    output = "data.frame")
            # x1 = subset(sleep,
            #             group == 1)$extra
            # y1 = subset(sleep,
            #             group == 2)$extra
            # m1 = mean(x1)
            # sd1 = sd(x1)
            # n1 = length(x1)
            # m2 = mean(y1)
            # sd2 = sd(y1)
            # n2 = length(y1)
            # r12 = cor(x1,y1)
            # df1 = data.frame(
            #   m1,
            #   sd1,
            #   n1,
            #   m2,
            #   sd2,
            #   n2,
            #   ri = r12
            # )
            #
            # (meta_d = escalc(measure = "SMD",
            #                  data = df1,
            #                  m1i = m1,
            #                  sd1i = sd1,
            #                  n1i = n1,
            #                  m2i = m2,
            #                  sd2i = sd2,
            #                  n2i = n2,
            #                  vtype = "LS2"))
            # meta_d$vi
            # (meta_c = escalc(measure = "SMCC",
            #                  data = df1,
            #                  m1i = m1,
            #                  sd1i = sd1,
            #                  ni = n1,
            #                  m2i = m2,
            #                  sd2i = sd2,
            #                  n2i = n2,
            #                  ri = ri,
            #                  vtype = "LS")
            # )
            # meta_c$vi
            # LS2 approximation
            # # Borenstein, 2009, equation 12.17
            expect_equal(0.2195277,round(smd$SE^2,7))

            expect_equal(0.1946978,round(smc$SE^2,7))

            smd = TOSTER::smd_calc(x = subset(sleep,
                                              group == 1)$extra,
                                   y = subset(sleep,
                                              group == 2)$extra,
                                   var.equal = FALSE,
                                   output = "data.frame")

            expect_equal(0.24,round(smd$SE^2,2))
          })

