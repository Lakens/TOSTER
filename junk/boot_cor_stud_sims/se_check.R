devtools::load_all("C:/GitHub/TOSTER", quiet = TRUE)
set.seed(7)
gens <- list(
  normal = function(n) { x <- rnorm(n); list(x = x, y = .6*x + .8*rnorm(n)) },
  hetero = function(n) { x <- rnorm(n); list(x = x, y = .5*x + rnorm(n)*sqrt(1+x^2)) },
  t3     = function(n) { x <- rt(n,3); list(x = x, y = .6*x + .8*rt(n,3)) },
  indep  = function(n) list(x = rnorm(n), y = rnorm(n)),
  discrete = function(n) { x <- sample(1:5, n, TRUE); list(x = x, y = pmin(5, x + sample(0:2, n, TRUE))) }
)
for (g in names(gens)) for (n in c(50, 300)) {
  res <- replicate(2000, { d <- gens[[g]](n)
    sapply(c("pearson","spearman","kendall"), function(m)
      c(cor(d$x, d$y, method = m), TOSTER:::.cor_se(d$x, d$y, m))) })
  mc_sd <- apply(res[1, , ], 1, sd); mean_se <- apply(res[2, , ], 1, mean)
  cat(sprintf("%-8s n=%3d  ", g, n),
      sprintf("%s: MCsd=%.4f meanSE=%.4f ratio=%.3f | ", c("P","S","K"), mc_sd, mean_se, mean_se/mc_sd), "\n")
}
# check kendall tau-b agreement
d <- gens$discrete(40); x <- d$x; y <- d$y; n <- 40
sx <- sign(outer(x,x,"-")); sy <- sign(outer(y,y,"-"))
cat("tau-b match:", all.equal(sum(sx*sy)/sqrt(sum(sx!=0)*sum(sy!=0)), cor(x,y,method="kendall")), "\n")
