cor_test = function(x,
                    y,
                    alternative = c("two.sided", "less", "greater"),
                    method = c("pearson", "kendall", "spearman"),
                    alpha = .05,
                    null = 0,
                    TOST = FALSE){
  alternative = match.arg(alternative)
  method = match.arg(method)

  if(TOST && null <=0){
    stop("positive value for null must be supplied if using TOST.")
  }
  if(TOST){
    alternative = "less"
  }
  r_xy = cor(x,y,
             method = method)
  df = data.frame(x=x,
                  y=y)
  df = na.omit(df)
  n_obs = nrow(df)

  z_xy = rho_to_z(r_xy)
  # get se ---
  if (method == "pearson") {
    # Pearson
    z.se <- 1 / sqrt(n_obs - 3)
  }
  if (method == "spearman") {
    # Kendall
    z.se <- (0.437 / (n_obs - 4)) ^ 0.5
  }
  if (method == "kendall") {
    # Spearman
    z.se <- (1.06 / (n_obs - 3)) ^ 0.5
  }
}
