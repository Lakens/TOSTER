devtools::load_all()
data(sleep)
x <- sleep$extra[sleep$group == 2]
y <- sleep$extra[sleep$group == 1]

cat("=== Two-sample score ===\n")
res_2s <- ses_calc(x=x, y=y, paired=FALSE, ses="cstat", se_method="score")
cat("Estimate:", as.numeric(res_2s$estimate), "\n")
cat("CI:", res_2s$conf.int, "\n\n")

cat("=== Paired score ===\n")
res_ps <- ses_calc(x=x, y=y, paired=TRUE, ses="cstat", se_method="score")
cat("Estimate:", as.numeric(res_ps$estimate), "\n")
cat("CI:", res_ps$conf.int, "\n")
cat("Label:", names(res_ps$estimate), "\n\n")

cat("=== Paired agresti (for comparison) ===\n")
res_pa <- suppressMessages(ses_calc(x=x, y=y, paired=TRUE, ses="cstat", se_method="agresti"))
cat("Estimate:", as.numeric(res_pa$estimate), "\n")
cat("CI:", res_pa$conf.int, "\n\n")

cat("=== Paired score, rb scale ===\n")
res_rb <- ses_calc(x=x, y=y, paired=TRUE, ses="rb", se_method="score")
cat("Estimate:", as.numeric(res_rb$estimate), "\n")
cat("CI:", res_rb$conf.int, "\n\n")

cat("=== Paired score p-value vs wilcox.test ===\n")
wt <- wilcox.test(x, y, paired=TRUE, exact=FALSE, correct=FALSE)
res_p <- ses_calc(x=x, y=y, paired=TRUE, ses="cstat", se_method="score",
                  alternative="two.sided", null.value=0.5, correct=FALSE)
cat("wilcox.test p:", wt$p.value, "\n")
cat("ses_calc p:", res_p$p.value, "\n")
cat("Match:", all.equal(res_p$p.value, wt$p.value, tolerance=1e-6), "\n\n")

cat("=== Simple non-boundary data ===\n")
set.seed(42)
x2 <- rnorm(20, mean=0.5)
y2 <- rnorm(20)
res_s <- ses_calc(x=x2, y=y2, paired=TRUE, ses="cstat", se_method="score")
res_a <- ses_calc(x=x2, y=y2, paired=TRUE, ses="cstat", se_method="agresti")
cat("Score estimate:", as.numeric(res_s$estimate), "\n")
cat("Agresti estimate:", as.numeric(res_a$estimate), "\n")
cat("Estimates same direction:", sign(as.numeric(res_s$estimate) - 0.5) == sign(as.numeric(res_a$estimate) - 0.5), "\n")
