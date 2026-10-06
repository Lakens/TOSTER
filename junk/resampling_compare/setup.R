# Resampling comparison: shared setup --------
#
# perm_t_test (studentized) vs boot_t_test (boot_ci = "stud" and "bca"), with an
# analytic reference (Welch t for tr = 0; Yuen's trimmed t for tr > 0), in
# two-sample and paired designs. Sourced by run.R and validate.R.
#
# Every method is evaluated through its confidence intervals, which (for all
# methods here) agree with the method's own p-values by construction. For each
# simulated data set we record, relative to the true value theta of the
# estimand (mean difference, or trimmed-mean difference when tr > 0):
#
#   cov95      95% two-sided CI covers theta (1 - cov95 = two-sided Type I error
#              of the standard test at theta)
#   err_below  90% (equal-tailed) CI lies entirely below theta
#   err_above  90% (equal-tailed) CI lies entirely above theta
#              These are one-sided alpha = 0.05 errors. If theta were an
#              equivalence bound, they are the TOST (and minimal effect) Type I
#              errors at that bound: at an upper bound, TOST error = err_below
#              and MET error = err_above; at a lower bound the roles swap.
#   width95    width of the 95% CI
#   rej0       0 lies outside the 95% CI (standard-test rejection of H0: 0;
#              power when theta != 0)
#   tost05, tost1  90% CI lies inside (-0.5, 0.5) / (-1, 1) (TOST decision with
#              those bounds; power when theta = 0)
#   fail       the method returned an error or a non-finite interval
#
# The 90% interval uses the same equal-tailed construction that each function
# uses for TOST (perm_t_test: symmetric = FALSE; boot_t_test: two-sided
# alpha = 0.10), checked in validate.R.

# Package loading --------

load_toster <- function(pkg = Sys.getenv("SIM_PKG", ".")) {
  if (nzchar(pkg) && file.exists(file.path(pkg, "DESCRIPTION"))) {
    suppressMessages(devtools::load_all(pkg, quiet = TRUE))
  } else {
    suppressMessages(library(TOSTER))
  }
  # Guard: the CI fix (issue #120) and tie tolerance must be present,
  # otherwise perm_t_test intervals are the old percentile intervals
  ns <- asNamespace("TOSTER")
  if (!exists("perm_crit", envir = ns) || !exists("perm_count", envir = ns)) {
    stop("This TOSTER version predates the perm_t_test CI fix (issue #120). ",
         "Point SIM_PKG to the package source on the fix branch or install it.")
  }
  invisible(TRUE)
}

# Distributions --------
# Continuous: standardized to mean 0, SD 1. r = generator, q = quantile
# function (only needed for asymmetric distributions, to get trimmed means).

ln_sig <- 0.5
ln_m <- exp(ln_sig^2 / 2)
ln_s <- sqrt((exp(ln_sig^2) - 1) * exp(ln_sig^2))

cont_dists <- list(
  normal = list(r = function(n) rnorm(n), q = qnorm, symmetric = TRUE),
  lognormal = list(r = function(n) (exp(ln_sig * rnorm(n)) - ln_m) / ln_s,
                   q = function(u) (exp(ln_sig * qnorm(u)) - ln_m) / ln_s,
                   symmetric = FALSE),
  exponential = list(r = function(n) rexp(n) - 1,
                     q = function(u) -log1p(-u) - 1,
                     symmetric = FALSE),
  t3 = list(r = function(n) rt(n, df = 3) / sqrt(3),
            q = function(u) qt(u, df = 3) / sqrt(3), symmetric = TRUE),
  # 10% contamination from N(0, 5^2); variance 0.9 + 0.1 * 25 = 3.4
  contaminated = list(r = function(n) ifelse(runif(n) < 0.1, rnorm(n, 0, 5),
                                             rnorm(n)) / sqrt(3.4),
                      q = NULL, symmetric = TRUE)
)

# Discrete: raw values (not standardized)
disc_dists <- list(
  likert_mid  = list(values = 1:5, probs = c(0.10, 0.20, 0.40, 0.20, 0.10)),
  likert_wide = list(values = 1:5, probs = c(0.30, 0.10, 0.20, 0.10, 0.30)),
  likert_skew = list(values = 1:5, probs = c(0.05, 0.45, 0.25, 0.15, 0.10)),
  diff_sym    = list(values = -2:2, probs = c(0.10, 0.20, 0.40, 0.20, 0.10)),
  diff_skew   = list(values = -2:2, probs = c(0.05, 0.10, 0.35, 0.40, 0.10))
)

rdist <- function(name, n) {
  if (name %in% names(cont_dists)) return(cont_dists[[name]]$r(n))
  d <- disc_dists[[name]]
  sample(d$values, n, replace = TRUE, prob = d$probs)
}

is_discrete <- function(name) name %in% names(disc_dists)

# Population trimmed means --------
# Trimmed mean with proportion tr trimmed from each tail:
#   (1 / (1 - 2 tr)) * integral_{tr}^{1 - tr} Q(u) du
# The sample estimator mean(x, trim = tr) is consistent for this quantity.

tm_cont <- function(name, tr) {
  d <- cont_dists[[name]]
  if (tr == 0 || d$symmetric) return(0)   # standardized means are 0
  stats::integrate(d$q, tr, 1 - tr, rel.tol = 1e-10)$value / (1 - 2 * tr)
}

tm_disc <- function(name, tr) {
  d <- disc_dists[[name]]
  if (tr == 0) return(sum(d$values * d$probs))
  Fhi <- cumsum(d$probs); Flo <- c(0, head(Fhi, -1))
  # length of each value's quantile segment inside [tr, 1 - tr]
  seg <- pmax(0, pmin(Fhi, 1 - tr) - pmax(Flo, tr))
  sum(d$values * seg) / (1 - 2 * tr)
}

tm_pop <- function(name, tr) {
  if (is_discrete(name)) tm_disc(name, tr) else tm_cont(name, tr)
}

# True value of the estimand for a cell
# Two-sample: x = shift + sdx * Zx, y = Zy (continuous) or raw discrete values
# Paired: differences d = shift + Zd
cell_truth <- function(cell) {
  if (cell$design == "paired") {
    return(cell$shift + tm_pop(cell$xdist, cell$tr))
  }
  cell$shift + cell$sdx * tm_pop(cell$xdist, cell$tr) - tm_pop(cell$ydist, cell$tr)
}

gen_data <- function(cell) {
  if (cell$design == "paired") {
    return(list(x = cell$shift + rdist(cell$xdist, cell$nx), y = NULL))
  }
  list(x = cell$shift + cell$sdx * rdist(cell$xdist, cell$nx),
       y = rdist(cell$ydist, cell$ny))
}

# Analytic reference intervals --------
# tr = 0: Welch t (two-sample) / one-sample t on differences (paired)
# tr > 0: Yuen's trimmed t (two-sample, Yuen-Welch df) / one-sample trimmed t

ref_ci <- function(x, y, tr, conf) {
  ns <- asNamespace("TOSTER")
  tm <- get("trimmed_mean", ns); wv <- get("winsorized_var", ns)
  en <- get("effective_n", ns); ywdf <- get("yuen_welch_df", ns)
  q <- function(df) stats::qt(1 - (1 - conf) / 2, df)
  if (is.null(y)) {
    n <- length(x)
    if (tr == 0) return(as.numeric(stats::t.test(x, conf.level = conf)$conf.int))
    h <- en(n, tr)
    se <- sqrt((n - 1) * wv(x, tr) / (h * (h - 1)))
    return(tm(x, tr) + c(-1, 1) * q(h - 1) * se)
  }
  if (tr == 0) return(as.numeric(stats::t.test(x, y, conf.level = conf)$conf.int))
  nx <- length(x); ny <- length(y)
  vx <- wv(x, tr); vy <- wv(y, tr)
  hx <- en(nx, tr); hy <- en(ny, tr)
  se <- sqrt((nx - 1) * vx / (hx * (hx - 1)) + (ny - 1) * vy / (hy * (hy - 1)))
  (tm(x, tr) - tm(y, tr)) + c(-1, 1) * q(ywdf(nx, ny, vx, vy, tr)) * se
}

# Package intervals --------
# For paired designs both functions are called with paired = TRUE and y = 0, so
# the differences are passed through exactly (no floating-point noise that
# would break ties in discrete data).

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

safe_ci <- function(expr) {
  ci <- tryCatch(as.numeric(quiet(expr)), error = function(e) c(NA_real_, NA_real_))
  if (length(ci) != 2 || any(!is.finite(ci))) c(NA_real_, NA_real_) else ci
}

method_cis <- function(x, y, tr, R, seed) {
  paired <- is.null(y)
  yy <- if (paired) rep(0, length(x)) else y

  perm_ci <- function(alpha, symmetric) {
    set.seed(seed)
    safe_ci(perm_t_test(x, yy, paired = paired, tr = tr, R = R, alpha = alpha,
                        symmetric = symmetric, keep_perm = FALSE)$conf.int)
  }
  boot_ci <- function(type, alpha) {
    set.seed(seed + 1)
    safe_ci(boot_t_test(x, yy, paired = paired, tr = tr, R = R, alpha = alpha,
                        boot_ci = type)$conf.int)
  }
  ref <- function(conf) {
    tryCatch(ref_ci(x, if (paired) NULL else y, tr, conf),
             error = function(e) c(NA_real_, NA_real_))
  }

  list(
    # 95%: package defaults (perm_t_test symmetric two-sided interval)
    ci95 = list(perm = perm_ci(0.05, TRUE),
                boot_stud = boot_ci("stud", 0.05),
                boot_bca = boot_ci("bca", 0.05),
                ref = ref(0.95)),
    # 90% equal-tailed: the interval each function uses for TOST
    ci90 = list(perm = perm_ci(0.10, FALSE),
                boot_stud = boot_ci("stud", 0.10),
                boot_bca = boot_ci("bca", 0.10),
                ref = ref(0.90))
  )
}

methods <- c("perm", "boot_stud", "boot_bca", "ref")
metric_names <- c("fail", "cov95", "err_below", "err_above", "width95", "rej0",
                  "tost05", "tost1")

# One replication --------

one_rep <- function(cell, truth, R, rep_seed) {
  set.seed(rep_seed)
  d <- gen_data(cell)
  cis <- method_cis(d$x, d$y, cell$tr, R, seed = rep_seed + 7919L)

  out <- lapply(methods, function(m) {
    a <- cis$ci95[[m]]; b <- cis$ci90[[m]]
    if (anyNA(a) || anyNA(b)) {
      return(c(fail = 1, setNames(rep(NA_real_, length(metric_names) - 1),
                                  metric_names[-1])))
    }
    c(fail = 0,
      cov95 = as.numeric(a[1] <= truth && truth <= a[2]),
      err_below = as.numeric(b[2] < truth),
      err_above = as.numeric(b[1] > truth),
      width95 = a[2] - a[1],
      rej0 = as.numeric(0 < a[1] || 0 > a[2]),
      tost05 = as.numeric(b[1] > -0.5 && b[2] < 0.5),
      tost1 = as.numeric(b[1] > -1 && b[2] < 1))
  })
  v <- unlist(out)
  names(v) <- paste(rep(methods, each = length(metric_names)), metric_names, sep = ".")
  v
}

# Cell definitions --------

# Layouts. x is always the "distinctive" group: the more variable one in
# "unequal SD" / "unequal spread" cells, and the skewed one in "different
# shape" cells. "x smaller" puts it in the smaller group (1:2), "x larger" in
# the larger group (2:1). For identical distributions only 1:2 is needed.
two_designs <- function(ns_bal = c(8, 15, 30, 50), ns_unbal = c(8, 15, 30),
                        unequal = FALSE) {
  bal <- data.frame(nx = ns_bal, ny = ns_bal, layout = "balanced")
  small_x <- data.frame(nx = ns_unbal, ny = 2 * ns_unbal,
                        layout = if (unequal) "x smaller (1:2)" else "1:2")
  if (!unequal) return(rbind(bal, small_x))
  large_x <- data.frame(nx = 2 * ns_unbal, ny = ns_unbal,
                        layout = "x larger (2:1)")
  rbind(bal, small_x, large_x)
}

make_two <- function(xdist, ydist, structure, sdx, tr, shift = 0, designs) {
  data.frame(design = "two", structure = structure, xdist = xdist, ydist = ydist,
             sdx = sdx, shift = shift, tr = tr, designs, stringsAsFactors = FALSE)
}

make_paired <- function(dist, tr, ns, shift = 0) {
  data.frame(design = "paired", structure = "paired differences", xdist = dist,
             ydist = NA_character_, sdx = 1, shift = shift, tr = tr,
             nx = ns, ny = ns, layout = "paired", stringsAsFactors = FALSE)
}

build_cells <- function(block, sd_ratio = 2) {
  conts <- names(cont_dists)
  ns_paired <- c(8, 15, 30, 50)

  core <- function(tr, include_discrete) {
    two <- do.call(rbind, c(
      lapply(conts, function(dn) make_two(dn, dn, "identical", 1, tr,
                                          designs = two_designs())),
      lapply(conts, function(dn) make_two(dn, dn, "unequal SD", sd_ratio, tr,
                                          designs = two_designs(unequal = TRUE)))))
    if (include_discrete) {
      two <- rbind(two,
        make_two("likert_mid", "likert_mid", "identical", 1, tr,
                 designs = two_designs()),
        make_two("likert_wide", "likert_mid", "unequal spread", 1, tr,
                 designs = two_designs(unequal = TRUE)))
    }
    pd <- if (include_discrete) c(conts, "diff_sym", "diff_skew") else conts
    paired <- do.call(rbind, lapply(pd, function(dn) make_paired(dn, tr, ns_paired)))
    rbind(two, paired)
  }

  cells <- switch(block,
    # Block A: core, no trimming
    A = core(tr = 0, include_discrete = TRUE),
    # Block B: core with 20% trimming (continuous distributions only)
    B = core(tr = 0.2, include_discrete = FALSE),
    # Block C: different shapes (weak null for means) and power
    C = rbind(
      do.call(rbind, lapply(c(0, 0.2), function(tr)
        make_two("exponential", "normal", "different shape", 1, tr,
                 designs = two_designs(unequal = TRUE)))),
      make_two("likert_skew", "likert_mid", "different shape", 1, 0,
               designs = two_designs(unequal = TRUE)),
      do.call(rbind, lapply(c(0, 0.2), function(tr) do.call(rbind, lapply(
        c("normal", "exponential", "t3", "contaminated"), function(dn) rbind(
          make_two(dn, dn, "power (shift 0.5)", 1, tr, shift = 0.5,
                   designs = data.frame(nx = c(15, 30), ny = c(15, 30),
                                        layout = "balanced")),
          make_paired(dn, tr, c(15, 30), shift = 0.5))))))),
    stop("block must be 'A', 'B', or 'C'"))

  cells$block <- block
  # stable cell id: block letter index * 1000 + row
  cells$cell_id <- match(block, LETTERS) * 1000L + seq_len(nrow(cells))
  cells$truth <- vapply(seq_len(nrow(cells)), function(i) cell_truth(cells[i, ]),
                        numeric(1))
  rownames(cells) <- NULL
  cells
}
