#context("Run Examples for boot_t_TOST")

# need hush function to run print through examples

hush = function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

test_that("compare_smd", {

  # Errors ----

  expect_error(compare_smd())

  # One sample
  set.seed(7894021)
  datx1 = rnorm(20)
  datx2 = rnorm(45)

  one_boot1 = boot_compare_smd(x1 = datx1,
                               x2 = datx2)
  one_1 = compare_smd(smd1 = one_boot1$df_ci$estimate[2],
                      n1 = length(datx1),
                      smd2 = one_boot1$df_ci$estimate[3],
                      n2 = length(datx2),
                      paired = TRUE)

  expect_equal(unname(one_1$statistic),
               unname(one_boot1$statistic))

  # Paired sample -----
  daty1 = rnorm(20)
  daty2 = rnorm(45)

  one_boot2 = boot_compare_smd(x1 = datx1,
                               x2 = datx2,
                               y1 = daty1,
                               y2 = daty2,
                               paired = TRUE)
  one_2 = compare_smd(smd1 = one_boot2$df_ci$estimate[2],
                      n1 = length(datx1),
                      smd2 = one_boot2$df_ci$estimate[3],
                      n2 = length(datx2),
                      paired = TRUE)

  expect_equal(unname(one_2$statistic),
               unname(one_boot2$statistic))

  # Two sample -----

  one_boot3 = boot_compare_smd(x1 = datx1,
                               x2 = datx2,
                               y1 = daty1,
                               y2 = daty2,
                               paired = FALSE)
  one_3 = compare_smd(smd1 = one_boot3$df_ci$estimate[2],
                      n1 = c(length(datx1),length(daty1)),
                      smd2 = one_boot3$df_ci$estimate[3],
                      n2 = c(length(datx2),length(daty2)),
                      paired = FALSE)

  expect_equal(unname(one_3$statistic),
               unname(one_boot3$statistic))

})

