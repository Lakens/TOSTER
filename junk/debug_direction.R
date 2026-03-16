devtools::load_all()
data(sleep)
x <- sleep$extra[sleep$group == 2]
y <- sleep$extra[sleep$group == 1]

cat("=== Direct rbs_calc calls ===\n")
cat("rbs_calc(x=x, y=y, paired=TRUE):", rbs_calc(x=x, y=y, mu=0, paired=TRUE), "\n")
cat("rbs_calc(x=y, y=x, paired=TRUE):", rbs_calc(x=y, y=x, mu=0, paired=TRUE), "\n")

cat("\n=== ses_calc paired, agresti ===\n")
res_agr <- suppressMessages(ses_calc(x=x, y=y, paired=TRUE, ses="cstat", se_method="agresti"))
cat("Agresti estimate:", as.numeric(res_agr$estimate), "\n")
cat("Agresti CI:", res_agr$conf.int, "\n")

cat("\n=== ses_calc paired, score ===\n")
res_sc <- ses_calc(x=x, y=y, paired=TRUE, ses="cstat", se_method="score")
cat("Score estimate:", as.numeric(res_sc$estimate), "\n")
cat("Score CI:", res_sc$conf.int, "\n")

cat("\n=== ses_calc two-sample, score ===\n")
res_2s <- ses_calc(x=x, y=y, paired=FALSE, ses="cstat", se_method="score")
cat("Two-sample estimate:", as.numeric(res_2s$estimate), "\n")
cat("Two-sample CI:", res_2s$conf.int, "\n")

cat("\n=== paired_rank_info calls ===\n")
pri_xy <- paired_rank_info(x, y, mu=0)
cat("paired_rank_info(x, y): p_hat =", pri_xy$p_hat, "\n")
pri_yx <- paired_rank_info(y, x, mu=0)
cat("paired_rank_info(y, x): p_hat =", pri_yx$p_hat, "\n")

cat("\n=== What ses_calc.default actually does ===\n")
# Line 432-433: rbs_calc(x = y, y = x, mu = mu, paired = TRUE)
r_via_swap <- rbs_calc(x=y, y=x, mu=0, paired=TRUE)
cat("rb via swap (rbs_calc(y,x)):", r_via_swap, "\n")
cat("p_hat via swap:", rb_to_cstat(r_via_swap), "\n")

# Without swap
r_no_swap <- rbs_calc(x=x, y=y, mu=0, paired=TRUE)
cat("rb no swap (rbs_calc(x,y)):", r_no_swap, "\n")
cat("p_hat no swap:", rb_to_cstat(r_no_swap), "\n")
