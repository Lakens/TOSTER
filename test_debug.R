devtools::load_all(quiet = TRUE)
x_sep <- 1:5
y_sep <- 6:10

# Check what ses_compute_agresti returns
agresti_res <- TOSTER:::ses_compute_agresti(x_sep, y_sep, paired = FALSE)
cat('agresti cstat:', agresti_res$cstat, '\n')
cat('agresti se_cstat:', agresti_res$se_cstat, '\n')
cat('agresti se_rb:', agresti_res$se_rb, '\n')
cat('boundary_corrected:', agresti_res$boundary_corrected, '\n')

# Now test the full ses_calc function
suppressMessages(res <- ses_calc(x_sep, y_sep, ses = "cstat", se_method = "agresti"))
cat('\nses_calc estimate:', res$estimate, '\n')
cat('ses_calc stderr:', res$stderr, '\n')
cat('ses_calc CI:', res$conf.int, '\n')

# Expected: estimate ≈ 0.5/26 ≈ 0.0192
cat('\nExpected estimate:', 0.5 / 26, '\n')

# Also test paired
x_paired <- c(1, 2, 3, 4, 5)
y_paired <- c(6, 7, 8, 9, 10)
suppressMessages(res_paired <- ses_calc(x_paired, y_paired, paired = TRUE, ses = "cstat", se_method = "agresti"))
cat('\nPaired estimate:', res_paired$estimate, '\n')
cat('Paired stderr:', res_paired$stderr, '\n')
