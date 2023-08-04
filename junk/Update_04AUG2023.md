# More Big Updates TOSTER
Aaron R. Caldwell

# version 0.8.0: August 2023

I am happy to announce that today version 0.8.0 of TOSTER has now been
merged to the main branch on GitHub, and is now available for downloads.
As always, detailed documentation of the package can be found its
[website](aaroncalwell.us/TOSTERpkg). There is a *lot* in this update,
and I want to use this post to highlight the key features.

# Brunner-Munzel Test

This is a fairly recently developed test that is thought by some to be
an improvement on the Wilcoxon-Mann-Whitney test. Essentially, this
function provides a test of the probability of superiority based on the
ranks of values. This makes it a valuable test for ordinal values (e.g.,
Likert-type items). Some have even suggested it become the “default”
non-parametric test (much like how Welch’s adjusted t-test is often the
default t-test for mean differences).

Since there is no base R support for this test, I had to create a basic
form of this test that allows for two-sided and one-sided tests. Notice
below that a warning message has been printed. The primary problem I was
able to find with the Brunner-Munzel is possible inflation in the type 1
error rate with small sample sizes.

``` r
brunner_munzel(data = sleep,
               extra ~ group)
```


        two-sample Brunner-Munzel test

    data:  extra by group
    t = -2.1447, df = 16.898, p-value = 0.04682
    alternative hypothesis: true relative effect is not equal to 0.5
    95 percent confidence interval:
     0.01387048 0.49612952
    sample estimates:
    p(X<Y) + .5*P(X=Y) 
                 0.255 

However, I have a solution which is the use of a permutation test.
Notice, it gives a slightly more conservative estimate.

``` r
set.seed(04082023)
brunner_munzel(data = sleep,
               extra ~ group,
               perm = TRUE)
```


        two-sample Brunner-Munzel permutation test

    data:  extra by group
    t = -2.1447, df = 16.898, p-value = 0.0522
    alternative hypothesis: true relative effect is not equal to 0.5
    95 percent confidence interval:
     0.008780122 0.502129784
    sample estimates:
    p(X<Y) + .5*P(X=Y) 
                 0.255 

## TOST with Brunner-Munzel

The TOST procedure for Brunner-Munzel can be accomplished with the
`simple_htest` function.

``` r
set.seed(04082023)
test1 = simple_htest(data = sleep,
             extra ~ group,
             mu = .7, # set equivalence bounds to .7 - .3
             alternative = "equ", # equivalence test
             test = "brunner", # brunner-munzel test
             perm = TRUE)

test1
```


        two-sample Brunner-Munzel permutation test

    data:  extra by group
    t = -3.8954, df = 16.898, p-value = 0.9972
    alternative hypothesis: equivalence
    null values:
    relative effect relative effect 
                0.3             0.7 
    90 percent confidence interval:
     0.05524003 0.45550152
    sample estimates:
    p(X<Y) + .5*P(X=Y) 
                 0.255 

The results can be described just like other htest objects using the
`describe_htest` function.

``` r
describe_htest(test1)
```

    [1] "The two-sample Brunner-Munzel permutation test is not statistically significant (t(16.898) = -3.9, p = 0.997, p(X<Y) + .5*P(X=Y) = 0.255, 90% C.I.[0.055, 0.456]) at a 0.05 alpha-level. The null hypothesis cannot be rejected. At the desired error rate, it cannot be stated that the true relative effect is between 0.3 and 0.7."

# Testing Two Proportions

I hesitated to develop anything new for testing differences in two
proportions, but the frequency of questions from users demanded that I
address this type of test.

Therefore, I have deprecated `TOSTtwo.prop` and have created the
`twoprop_test` function. This operates like `prop.test` from base R, but
adds more functionality/tests.

First, users are able to test hypotheses based on the difference, the
risk ratio, or the odds ratio (please read the documentation for more
details).

Second, users can test a two-tailed test, a one-tailed test, or perform
a TOST test using the alternative argument.

``` r
twoprop_test(
  p1 = 30/130,
  p2 = 37/130,
  n1 = 130,
  n2 = 130,
  null = 1.25, # null on odds ratio scale
  alpha = 0.05, # .05 alpha level
  alternative = "equ", # equivalence test
  effect_size = "odds"
)
```


        approximate Odds Ratio z-test

    data:  
    z = -0.20899, p-value = 0.5828
    alternative hypothesis: equivalence
    null values:
    Odds Ratio Odds Ratio 
          1.25       0.80 
    90 percent confidence interval:
     0.4734004 1.2010923
    sample estimates:
    Odds Ratio 
     0.7540541 

Please note, these are approximate tests which makes it a fairly liberal
test when samples sizes are very small.

## Power for Test of Two Proportions

A power analysis function is also available for differences in two
proportions. However, at this time it does not allow for power analysis
based on odds ratio or risk ratio.

``` r
# Example from documentation

power_twoprop(
  alpha = 0.01,
  power = 0.8,
  p1 = 0.5,
  p2 = 0.5,
  null = 0.2,
  alternative = "e"
)
```


         Power for Test of Differences in Two Proportions (z-test) 

                  n = 162.7117
        proportions = 0.5, 0.5
              alpha = 0.01
               beta = 0.2
              power = 0.8
               null = 0.2, -0.2
        alternative = equivalence
               NOTE = Sample sizes for EACH group

# Power for Correlations

I have also retired `powerTOSTr` and replaced it with `power_z_cor` to
allow for power analyses based on the tests in `z_cor_test`.

``` r
## Sample size for alpha = 0.05, 90% power, equivalence bounds of
## r = -0.1 and r = 0.1, assuming true effect = 0
#powerTOSTr(alpha=0.05, statistical_power=0.9, low_eqbound_r=-0.1, high_eqbound_r=0.1)
power_z_cor(
  alternative = "equivalence",
  alpha = .05,
  null = .1,
  power = .9,
  rho = 0
)
```


         Approximate Power for Pearson Product-Moment Correlation (z-test) 

                  n = 1077.993
                rho = 0
              alpha = 0.05
               beta = 0.1
              power = 0.9
               null = 0.1, -0.1
        alternative = equivalence

# lifecycle documentation

Lastly, I have used the [lifecyle](https://lifecycle.r-lib.org/) R
package to document the status of each function in the package.

# Concluding Remarks

Once again, I am pleased with the new developments I have made in
TOSTER, and I am excited to see how people react/use the new functions.
I have a “Contact Me” form on my
[website](https://aaroncaldwell.us/#contact), and please feel free to
send a message at any time. I would appreciate any feedback on the
describe functions in particular!

I hope you all enjoy the new TOSTER!

Cheers everyone,

AC
