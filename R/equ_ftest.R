

equ_ftest <- function(Fstat,
                      pes,
                      df1,
                      df2,
                      eqbound,
                      alpha = 0.05){

  f2 = eqbound/(1 - eqbound)
  lambda = (f2 * (df1 + df2 + 1))

  conf_level = 1 - alpha

  pval = pf(Fstat,
            df1 = df1,
            df2 = df2,
            ncp = lambda,
            lower.tail = TRUE)
  cint = c(0,1)
  attr(cint,"conf.level") <- conf_level
  parameter = c(df1, df2)
  names(dfs) <- "df"

  alt_val = paste0("partial eta-squared less than ", eqbound)

  rval <- list(statistic = Fstat,
               parameter = c(df1, df2),
               p.value = pval,
               conf.int = cint,
               estimate = pes,
               null.value = eqbound,
               alternative = alt_val,
               method = "Equivalence Test from F-test",
               data.name = "Summary Statistics")
  class(rval) <- "htest"
  return(rval)
}

equ_ftest(0.3066,
          .0080039,
          1,
          38,
          0.5)

library(afex)

data(sk2011.1)
a1 <- aov_ez("id", "response", sk2011.1, between = "instruction",
             anova_table = list(es = "pes"))
anova_table = a1$anova_table
colnames(anova_table) <- c("num_Df", "den_Df", "MSE", "F", "pes", "p")
anova_table$f2 <- anova_table$pes/(1 - anova_table$pes)

anova_table$lambda <- (anova_table$f2 * (anova_table$den_Df + anova_table$num_Df + 1))

(anova_table$den_Df+2) * anova_table$pes / (1 - anova_table$pes)
