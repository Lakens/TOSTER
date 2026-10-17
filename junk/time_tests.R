# Time the test suite under each tier --------
# Usage: Rscript junk/time_tests.R [cran|full|extended] [outfile.csv]

args <- commandArgs(trailingOnly = TRUE)
tier <- if (length(args) >= 1) args[1] else "cran"
out <- if (length(args) >= 2) args[2] else file.path("junk", paste0("test_times_", tier, ".csv"))

Sys.setenv(
  NOT_CRAN = if (tier == "cran") "false" else "true",
  TOSTER_EXTENDED_TESTS = if (tier == "extended") "true" else "false"
)

t0 <- proc.time()
res <- testthat::test_local(".", reporter = "summary", stop_on_failure = FALSE)
total <- (proc.time() - t0)[["elapsed"]]

df <- as.data.frame(res)
write.csv(df[, c("file", "test", "real", "skipped", "failed", "error")], out, row.names = FALSE)

## Summaries ----
by_file <- aggregate(real ~ file, data = df, FUN = sum)
by_file <- by_file[order(-by_file$real), ]
cat("\n\nTier:", tier, "| total elapsed:", round(total, 1), "s\n")
cat("Failed:", sum(df$failed), "| Errors:", sum(df$error), "| Skipped:", sum(df$skipped), "\n\n")
print(head(by_file, 20), row.names = FALSE)
cat("\nSlowest tests:\n")
slow <- df[order(-df$real), c("file", "test", "real")]
print(head(slow, 25), row.names = FALSE)
