
# 	TOSTER: Two one-sided tests (TOST) equivalence testing
 <!-- badges: start -->
  [![R-CMD-check](https://github.com/Lakens/TOSTER/workflows/R-CMD-check/badge.svg)](https://github.com/Lakens/TOSTER/actions)
  [![Coverage status](https://codecov.io/gh/Lakens/TOSTER/branch/master/graph/badge.svg)](https://codecov.io/github/Lakens/TOSTER?branch=master)
  [![CRAN_Status_Badge_version_last_release](https://www.r-pkg.org/badges/version-last-release/TOSTER)](https://cran.r-project.org/package=TOSTER)
  [![documentation](https://img.shields.io/badge/website-active-blue)](https://aaroncaldwell.us/TOSTERpkg/)
  <!-- badges: end -->

*An R package and jamovi module for equivalence testing*

Please see the package's [website](https://aaroncaldwell.us/TOSTERpkg/) for updates, vignettes, and other details about the package.

**Background**

Scientists should be able to provide support for the absence of a meaningful effect. Currently, researchers often incorrectly conclude an effect is absent based a non-significant result. A widely recommended approach within a frequentist framework is to test for equivalence. In equivalence tests, such as the two one-sided tests (TOST) procedure implemented in this package, an upper and lower equivalence bound is specified based on the smallest effect size of interest. The TOST procedure can be used to statistically reject the presence of effects large enough to be considered worthwhile. Extending your statistical tool kit with equivalence tests is an easy way to improve your statistical and theoretical inferences.

# Installation

The developmental version, and most up-to-date, can be installed from [GitHub](https://github.com/Lakens/TOSTER):

```
devtools::install_github("Lakens/TOSTER")
```

The stable release can be downloaded from [CRAN](https://cran.r-project.org/package=TOSTER):

```
install.packages("TOSTER")
```

## Educational Material

For educational material on setting the smallest effect size of interest and equivalence tests, see week 2 of my MOOC "Improving Your Statistical Questions". https://www.coursera.org/teach/improving-statistical-questions. 

For an introduction to equivalence testing and the TOSTER package (this is recommended reading material to understand the basics of equivalence tests) see: 
Lakens, D. (2017). Equivalence Tests: A Practical Primer for t Tests, Correlations, and Meta-Analyses. Social Psychological and Personality Science, 8(4), 355–362. https://doi.org/10.1177/1948550617697177

For a tutorial paper (this is the recommended reading material if you want to start using equivalence testing) see:
Lakens, D., Scheel, A. M., & Isager, P. M. (2018). Equivalence Testing for Psychological Research: A Tutorial. Advances in Methods and Practices in Psychological Science, 1(2), 259–269. https://doi.org/10.1177/2515245918770963

For a comparison of Bayes factors and equivalence test (they turn out to lead to very similar inferences when used well) see: 
Lakens, D., McLatchie, N., Isager, P. M., Scheel, A. M., & Dienes, Z. (2018). Improving Inferences about Null Effects with Bayes Factors and Equivalence Tests. The Journals of Gerontology: Series B. https://doi.org/10.1093/geronb/gby065

For a comparison of equivalence tests and second generation p-values (equivalence tests are probably a better tool) see: 
Lakens, D., & Delacre, M. (2019). Equivalence Testing and the Second Generation P-Value. Meta-Psychology. https://doi.org/10.31234/osf.io/7k6ay

For a general introduction to the importance of being able to support 'null' effects, and ways to do this, including equivalence tests, bayesian estimation, and bayes factors, see:
Harms, C., & Lakens, D. (2018). Making "null effects" informative: Statistical techniques and inferential frameworks. Journal of Clinical and Translational Research, (3), 382–393. https://doi.org/10.18053/jctres.03.2017S2.007
