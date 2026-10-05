devtools::load_all("C:/GitHub/TOSTER", quiet = TRUE)
set.seed(30)
gens <- list(
  normal = function(n) { x <- rnorm(n); list(x = x, y = .5*x + sqrt(.75)*rnorm(n)) },
  t5 = function(n) { x <- rt(n, 5); list(x = x, y = .6*x + .8*rt(n, 5)) },
  discrete = function(n) { x <- sample(1:5, n, TRUE); list(x = x, y = pmin(5, x + sample(0:2, n, TRUE))) })
for (g in names(gens)) for (n in c(20, 30)) {
  res <- replicate(4000, { d <- gens[[g]](n)
    c(cor(d$x, d$y, method = "kendall"), TOSTER:::.cor_se(d$x, d$y, "kendall")) })
  tau <- res[1, ]; se <- res[2, ]
  zsd <- sd(atanh(tau)); zse <- se / (1 - tau^2)
  cat(sprintf("%-8s n=%d  tau-scale: SE/MCsd=%.3f  z-scale: SE/MCsd=%.3f  cv(SE)=%.2f  cor(tau,SE)=%.2f\n",
              g, n, mean(se) / sd(tau), mean(zse) / zsd, sd(se) / mean(se), cor(tau, se)))
}
