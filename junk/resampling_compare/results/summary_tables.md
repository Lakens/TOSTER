# Resampling comparison: summary tables

nsim = 2000 per cell; MC SE near .05 = 0.0049. Rates more than 2 MC SE above .05 are **bold** (liberal); more than 2 MC SE below are _italic_ (conservative).

## Block A, two-sample: two-sided Type I error at the true value

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| identical | contaminated | 1:2 | 8 v 16 | 0 | 0.058 | **0.083** | **0.137** | 0.047 |
| identical | contaminated | 1:2 | 15 v 30 | 0 | 0.046 | **0.085** | **0.126** | _0.037_ |
| identical | contaminated | 1:2 | 30 v 60 | 0 | 0.054 | **0.086** | **0.103** | 0.049 |
| identical | contaminated | balanced | 8 v 8 | 0 | 0.053 | **0.075** | **0.150** | _0.038_ |
| identical | contaminated | balanced | 15 v 15 | 0 | 0.048 | **0.089** | **0.140** | _0.036_ |
| identical | contaminated | balanced | 30 v 30 | 0 | 0.044 | **0.079** | **0.104** | _0.040_ |
| identical | contaminated | balanced | 50 v 50 | 0 | **0.061** | **0.089** | **0.098** | 0.057 |
| identical | exponential | 1:2 | 8 v 16 | 0 | 0.044 | **0.068** | **0.110** | 0.046 |
| identical | exponential | 1:2 | 15 v 30 | 0 | 0.050 | **0.068** | **0.088** | 0.052 |
| identical | exponential | 1:2 | 30 v 60 | 0 | 0.040 | 0.059 | **0.068** | 0.043 |
| identical | exponential | balanced | 8 v 8 | 0 | 0.043 | 0.048 | **0.121** | _0.028_ |
| identical | exponential | balanced | 15 v 15 | 0 | 0.054 | **0.073** | **0.102** | 0.046 |
| identical | exponential | balanced | 30 v 30 | 0 | 0.052 | **0.068** | **0.080** | 0.049 |
| identical | exponential | balanced | 50 v 50 | 0 | 0.043 | 0.057 | **0.067** | 0.043 |
| identical | likert_mid | 1:2 | 8 v 16 | 0 | 0.046 | 0.054 | **0.082** | 0.054 |
| identical | likert_mid | 1:2 | 15 v 30 | 0 | 0.041 | _0.038_ | 0.051 | 0.040 |
| identical | likert_mid | 1:2 | 30 v 60 | 0 | 0.048 | 0.048 | 0.052 | 0.048 |
| identical | likert_mid | balanced | 8 v 8 | 0 | _0.026_ | 0.051 | **0.090** | 0.045 |
| identical | likert_mid | balanced | 15 v 15 | 0 | _0.033_ | 0.046 | **0.063** | 0.046 |
| identical | likert_mid | balanced | 30 v 30 | 0 | 0.042 | 0.053 | 0.057 | 0.058 |
| identical | likert_mid | balanced | 50 v 50 | 0 | 0.046 | 0.054 | **0.060** | 0.057 |
| identical | lognormal | 1:2 | 8 v 16 | 0 | 0.048 | **0.064** | **0.112** | 0.049 |
| identical | lognormal | 1:2 | 15 v 30 | 0 | 0.047 | **0.064** | **0.083** | 0.048 |
| identical | lognormal | 1:2 | 30 v 60 | 0 | 0.046 | 0.056 | **0.068** | 0.046 |
| identical | lognormal | balanced | 8 v 8 | 0 | 0.047 | 0.054 | **0.122** | _0.035_ |
| identical | lognormal | balanced | 15 v 15 | 0 | 0.058 | **0.065** | **0.090** | 0.053 |
| identical | lognormal | balanced | 30 v 30 | 0 | 0.045 | 0.053 | **0.068** | 0.041 |
| identical | lognormal | balanced | 50 v 50 | 0 | 0.048 | **0.062** | **0.068** | 0.049 |
| identical | normal | 1:2 | 8 v 16 | 0 | 0.056 | 0.052 | **0.093** | 0.058 |
| identical | normal | 1:2 | 15 v 30 | 0 | 0.041 | 0.042 | 0.057 | 0.040 |
| identical | normal | 1:2 | 30 v 60 | 0 | 0.044 | 0.047 | 0.059 | 0.046 |
| identical | normal | balanced | 8 v 8 | 0 | 0.049 | 0.050 | **0.096** | 0.047 |
| identical | normal | balanced | 15 v 15 | 0 | 0.051 | 0.050 | **0.068** | 0.051 |
| identical | normal | balanced | 30 v 30 | 0 | 0.055 | 0.056 | **0.069** | 0.054 |
| identical | normal | balanced | 50 v 50 | 0 | 0.047 | 0.050 | 0.054 | 0.049 |
| identical | t3 | 1:2 | 8 v 16 | 0 | 0.053 | **0.074** | **0.122** | 0.047 |
| identical | t3 | 1:2 | 15 v 30 | 0 | 0.053 | **0.073** | **0.100** | 0.049 |
| identical | t3 | 1:2 | 30 v 60 | 0 | 0.045 | **0.071** | **0.089** | 0.040 |
| identical | t3 | balanced | 8 v 8 | 0 | 0.052 | **0.069** | **0.128** | 0.041 |
| identical | t3 | balanced | 15 v 15 | 0 | 0.050 | **0.071** | **0.102** | _0.040_ |
| identical | t3 | balanced | 30 v 30 | 0 | 0.043 | **0.062** | **0.081** | 0.040 |
| identical | t3 | balanced | 50 v 50 | 0 | 0.045 | **0.061** | **0.070** | 0.042 |
| unequal SD | contaminated (SD x2) | balanced | 8 v 8 | 0 | 0.045 | 0.053 | **0.128** | _0.026_ |
| unequal SD | contaminated (SD x2) | balanced | 15 v 15 | 0 | 0.057 | **0.083** | **0.126** | 0.043 |
| unequal SD | contaminated (SD x2) | balanced | 30 v 30 | 0 | 0.055 | **0.085** | **0.119** | 0.046 |
| unequal SD | contaminated (SD x2) | balanced | 50 v 50 | 0 | 0.055 | **0.086** | **0.102** | 0.053 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.046 | **0.090** | **0.143** | _0.035_ |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.043 | **0.088** | **0.119** | _0.036_ |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.053 | **0.078** | **0.089** | 0.045 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 8 v 16 | 0 | 0.058 | **0.080** | **0.143** | 0.043 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 15 v 30 | 0 | 0.050 | **0.082** | **0.124** | _0.035_ |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 30 v 60 | 0 | 0.048 | **0.092** | **0.118** | _0.037_ |
| unequal SD | exponential (SD x2) | balanced | 8 v 8 | 0 | **0.081** | **0.084** | **0.141** | **0.071** |
| unequal SD | exponential (SD x2) | balanced | 15 v 15 | 0 | **0.068** | **0.075** | **0.103** | **0.063** |
| unequal SD | exponential (SD x2) | balanced | 30 v 30 | 0 | **0.062** | **0.067** | **0.083** | 0.057 |
| unequal SD | exponential (SD x2) | balanced | 50 v 50 | 0 | 0.053 | 0.059 | **0.068** | 0.053 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.043 | **0.064** | **0.105** | 0.046 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.053 | **0.072** | **0.090** | 0.055 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.045 | 0.057 | **0.065** | 0.049 |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.088** | **0.082** | **0.121** | **0.083** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.078** | **0.066** | **0.098** | **0.079** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 30 v 60 | 0 | **0.068** | 0.058 | **0.073** | **0.070** |
| unequal SD | lognormal (SD x2) | balanced | 8 v 8 | 0 | **0.081** | **0.084** | **0.133** | **0.073** |
| unequal SD | lognormal (SD x2) | balanced | 15 v 15 | 0 | 0.059 | **0.065** | **0.102** | 0.054 |
| unequal SD | lognormal (SD x2) | balanced | 30 v 30 | 0 | **0.060** | **0.064** | **0.077** | 0.057 |
| unequal SD | lognormal (SD x2) | balanced | 50 v 50 | 0 | 0.056 | **0.061** | **0.071** | 0.054 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 16 v 8 | 0 | _0.038_ | 0.056 | **0.086** | _0.039_ |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 30 v 15 | 0 | _0.036_ | 0.055 | **0.074** | 0.041 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.049 | **0.061** | **0.067** | 0.050 |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.075** | **0.068** | **0.112** | **0.070** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.074** | **0.064** | **0.097** | **0.070** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 30 v 60 | 0 | 0.058 | 0.056 | **0.071** | 0.057 |
| unequal SD | normal (SD x2) | balanced | 8 v 8 | 0 | 0.055 | 0.049 | **0.103** | 0.051 |
| unequal SD | normal (SD x2) | balanced | 15 v 15 | 0 | 0.055 | 0.048 | **0.072** | 0.050 |
| unequal SD | normal (SD x2) | balanced | 30 v 30 | 0 | **0.061** | **0.060** | **0.076** | 0.058 |
| unequal SD | normal (SD x2) | balanced | 50 v 50 | 0 | 0.046 | 0.042 | 0.051 | 0.045 |
| unequal SD | normal (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.045 | 0.050 | **0.076** | 0.049 |
| unequal SD | normal (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.047 | 0.050 | **0.064** | 0.052 |
| unequal SD | normal (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.050 | 0.049 | **0.061** | 0.049 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.068** | 0.058 | **0.100** | **0.060** |
| unequal SD | normal (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.064** | **0.061** | **0.089** | 0.056 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 30 v 60 | 0 | 0.052 | 0.049 | **0.066** | 0.051 |
| unequal SD | t3 (SD x2) | balanced | 8 v 8 | 0 | 0.054 | **0.064** | **0.132** | 0.040 |
| unequal SD | t3 (SD x2) | balanced | 15 v 15 | 0 | 0.050 | **0.066** | **0.112** | 0.041 |
| unequal SD | t3 (SD x2) | balanced | 30 v 30 | 0 | 0.048 | **0.068** | **0.085** | _0.038_ |
| unequal SD | t3 (SD x2) | balanced | 50 v 50 | 0 | 0.057 | **0.075** | **0.085** | 0.052 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.042 | **0.066** | **0.106** | _0.038_ |
| unequal SD | t3 (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.054 | **0.076** | **0.099** | 0.050 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.053 | **0.079** | **0.090** | 0.050 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.070** | **0.077** | **0.146** | 0.046 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 15 v 30 | 0 | 0.049 | **0.068** | **0.093** | 0.040 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 30 v 60 | 0 | 0.043 | **0.064** | **0.082** | _0.034_ |
| unequal spread | likert_wide vs likert_mid | balanced | 8 v 8 | 0 | _0.030_ | _0.025_ | **0.064** | 0.049 |
| unequal spread | likert_wide vs likert_mid | balanced | 15 v 15 | 0 | 0.041 | _0.034_ | **0.061** | 0.054 |
| unequal spread | likert_wide vs likert_mid | balanced | 30 v 30 | 0 | _0.037_ | _0.038_ | 0.053 | 0.047 |
| unequal spread | likert_wide vs likert_mid | balanced | 50 v 50 | 0 | 0.048 | 0.051 | 0.054 | 0.053 |
| unequal spread | likert_wide vs likert_mid | x larger (2:1) | 16 v 8 | 0 | 0.044 | _0.038_ | **0.074** | 0.052 |
| unequal spread | likert_wide vs likert_mid | x larger (2:1) | 30 v 15 | 0 | 0.046 | _0.040_ | 0.056 | 0.050 |
| unequal spread | likert_wide vs likert_mid | x larger (2:1) | 60 v 30 | 0 | 0.044 | 0.045 | 0.052 | 0.048 |
| unequal spread | likert_wide vs likert_mid | x smaller (1:2) | 8 v 16 | 0 | 0.050 | _0.024_ | 0.059 | 0.052 |
| unequal spread | likert_wide vs likert_mid | x smaller (1:2) | 15 v 30 | 0 | 0.052 | _0.028_ | 0.047 | 0.051 |
| unequal spread | likert_wide vs likert_mid | x smaller (1:2) | 30 v 60 | 0 | 0.054 | 0.041 | 0.054 | 0.054 |

## Block A, two-sample: one-sided error (TOST / minimal effect Type I error at a bound, worst direction)

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| identical | contaminated | 1:2 | 8 v 16 | 0 | 0.055 | **0.078** | **0.104** | 0.050 |
| identical | contaminated | 1:2 | 15 v 30 | 0 | 0.055 | **0.079** | **0.098** | 0.047 |
| identical | contaminated | 1:2 | 30 v 60 | 0 | 0.051 | **0.075** | **0.082** | 0.053 |
| identical | contaminated | balanced | 8 v 8 | 0 | 0.055 | **0.077** | **0.105** | 0.044 |
| identical | contaminated | balanced | 15 v 15 | 0 | 0.053 | **0.089** | **0.112** | 0.048 |
| identical | contaminated | balanced | 30 v 30 | 0 | 0.045 | **0.070** | **0.082** | 0.046 |
| identical | contaminated | balanced | 50 v 50 | 0 | 0.054 | **0.073** | **0.081** | 0.052 |
| identical | exponential | 1:2 | 8 v 16 | 0 | 0.049 | **0.080** | **0.096** | **0.068** |
| identical | exponential | 1:2 | 15 v 30 | 0 | 0.052 | **0.078** | **0.084** | **0.074** |
| identical | exponential | 1:2 | 30 v 60 | 0 | 0.053 | 0.056 | **0.062** | 0.054 |
| identical | exponential | balanced | 8 v 8 | 0 | 0.050 | **0.062** | **0.100** | 0.042 |
| identical | exponential | balanced | 15 v 15 | 0 | 0.050 | **0.069** | **0.087** | 0.048 |
| identical | exponential | balanced | 30 v 30 | 0 | 0.050 | 0.058 | **0.067** | 0.049 |
| identical | exponential | balanced | 50 v 50 | 0 | 0.048 | 0.056 | **0.063** | 0.048 |
| identical | likert_mid | 1:2 | 8 v 16 | 0 | 0.044 | 0.053 | **0.077** | 0.053 |
| identical | likert_mid | 1:2 | 15 v 30 | 0 | 0.045 | 0.045 | 0.058 | 0.045 |
| identical | likert_mid | 1:2 | 30 v 60 | 0 | 0.052 | 0.051 | 0.052 | 0.055 |
| identical | likert_mid | balanced | 8 v 8 | 0 | _0.029_ | 0.056 | **0.098** | 0.056 |
| identical | likert_mid | balanced | 15 v 15 | 0 | _0.036_ | 0.048 | **0.070** | 0.050 |
| identical | likert_mid | balanced | 30 v 30 | 0 | 0.042 | 0.051 | **0.060** | 0.054 |
| identical | likert_mid | balanced | 50 v 50 | 0 | 0.045 | 0.052 | **0.060** | 0.052 |
| identical | lognormal | 1:2 | 8 v 16 | 0 | 0.048 | **0.080** | **0.103** | **0.072** |
| identical | lognormal | 1:2 | 15 v 30 | 0 | 0.050 | 0.058 | **0.066** | 0.056 |
| identical | lognormal | 1:2 | 30 v 60 | 0 | 0.051 | **0.066** | **0.072** | **0.066** |
| identical | lognormal | balanced | 8 v 8 | 0 | 0.052 | **0.061** | **0.090** | 0.048 |
| identical | lognormal | balanced | 15 v 15 | 0 | 0.058 | **0.065** | **0.084** | 0.056 |
| identical | lognormal | balanced | 30 v 30 | 0 | 0.048 | 0.054 | 0.059 | 0.047 |
| identical | lognormal | balanced | 50 v 50 | 0 | 0.051 | **0.060** | **0.063** | 0.053 |
| identical | normal | 1:2 | 8 v 16 | 0 | **0.060** | 0.057 | **0.084** | **0.061** |
| identical | normal | 1:2 | 15 v 30 | 0 | 0.046 | 0.044 | 0.054 | 0.047 |
| identical | normal | 1:2 | 30 v 60 | 0 | 0.052 | 0.051 | 0.054 | 0.051 |
| identical | normal | balanced | 8 v 8 | 0 | 0.059 | 0.059 | **0.085** | 0.058 |
| identical | normal | balanced | 15 v 15 | 0 | 0.054 | 0.051 | **0.066** | 0.056 |
| identical | normal | balanced | 30 v 30 | 0 | 0.054 | 0.055 | **0.062** | 0.053 |
| identical | normal | balanced | 50 v 50 | 0 | 0.054 | 0.052 | 0.057 | 0.052 |
| identical | t3 | 1:2 | 8 v 16 | 0 | 0.056 | **0.070** | **0.092** | 0.050 |
| identical | t3 | 1:2 | 15 v 30 | 0 | 0.051 | **0.070** | **0.086** | 0.050 |
| identical | t3 | 1:2 | 30 v 60 | 0 | 0.053 | **0.064** | **0.072** | 0.050 |
| identical | t3 | balanced | 8 v 8 | 0 | 0.053 | **0.065** | **0.099** | 0.050 |
| identical | t3 | balanced | 15 v 15 | 0 | 0.051 | **0.066** | **0.083** | 0.048 |
| identical | t3 | balanced | 30 v 30 | 0 | 0.048 | **0.060** | **0.069** | 0.048 |
| identical | t3 | balanced | 50 v 50 | 0 | 0.050 | **0.065** | **0.070** | 0.048 |
| unequal SD | contaminated (SD x2) | balanced | 8 v 8 | 0 | 0.050 | **0.066** | **0.098** | _0.038_ |
| unequal SD | contaminated (SD x2) | balanced | 15 v 15 | 0 | 0.054 | **0.078** | **0.102** | 0.048 |
| unequal SD | contaminated (SD x2) | balanced | 30 v 30 | 0 | 0.053 | **0.080** | **0.094** | 0.050 |
| unequal SD | contaminated (SD x2) | balanced | 50 v 50 | 0 | 0.056 | **0.078** | **0.085** | 0.056 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.041 | **0.086** | **0.106** | 0.047 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.042 | **0.080** | **0.093** | 0.048 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.045 | **0.073** | **0.082** | 0.052 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.072** | **0.070** | **0.112** | 0.044 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.066** | **0.075** | **0.102** | 0.044 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 30 v 60 | 0 | **0.069** | **0.086** | **0.096** | 0.050 |
| unequal SD | exponential (SD x2) | balanced | 8 v 8 | 0 | **0.110** | **0.107** | **0.141** | **0.106** |
| unequal SD | exponential (SD x2) | balanced | 15 v 15 | 0 | **0.088** | **0.078** | **0.094** | **0.086** |
| unequal SD | exponential (SD x2) | balanced | 30 v 30 | 0 | **0.077** | **0.068** | **0.075** | **0.077** |
| unequal SD | exponential (SD x2) | balanced | 50 v 50 | 0 | **0.066** | 0.056 | 0.056 | **0.064** |
| unequal SD | exponential (SD x2) | x larger (2:1) | 16 v 8 | 0 | **0.077** | **0.070** | **0.092** | **0.062** |
| unequal SD | exponential (SD x2) | x larger (2:1) | 30 v 15 | 0 | **0.080** | **0.070** | **0.080** | **0.066** |
| unequal SD | exponential (SD x2) | x larger (2:1) | 60 v 30 | 0 | **0.068** | 0.059 | **0.062** | 0.059 |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.098** | **0.099** | **0.125** | **0.113** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.093** | **0.080** | **0.098** | **0.102** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 30 v 60 | 0 | **0.074** | **0.066** | **0.071** | **0.086** |
| unequal SD | lognormal (SD x2) | balanced | 8 v 8 | 0 | **0.097** | **0.090** | **0.118** | **0.092** |
| unequal SD | lognormal (SD x2) | balanced | 15 v 15 | 0 | **0.084** | **0.078** | **0.094** | **0.082** |
| unequal SD | lognormal (SD x2) | balanced | 30 v 30 | 0 | **0.072** | **0.069** | **0.071** | **0.074** |
| unequal SD | lognormal (SD x2) | balanced | 50 v 50 | 0 | **0.063** | 0.059 | **0.066** | **0.061** |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 16 v 8 | 0 | **0.066** | 0.058 | **0.075** | 0.056 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 30 v 15 | 0 | **0.064** | 0.056 | **0.066** | 0.053 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.056 | 0.052 | 0.058 | 0.050 |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.092** | **0.084** | **0.112** | **0.092** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.087** | **0.080** | **0.096** | **0.096** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 30 v 60 | 0 | **0.071** | **0.067** | **0.074** | **0.078** |
| unequal SD | normal (SD x2) | balanced | 8 v 8 | 0 | 0.054 | 0.052 | **0.085** | 0.050 |
| unequal SD | normal (SD x2) | balanced | 15 v 15 | 0 | **0.064** | 0.058 | **0.070** | **0.061** |
| unequal SD | normal (SD x2) | balanced | 30 v 30 | 0 | 0.058 | **0.060** | **0.070** | 0.058 |
| unequal SD | normal (SD x2) | balanced | 50 v 50 | 0 | 0.048 | 0.047 | 0.052 | 0.045 |
| unequal SD | normal (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.050 | 0.056 | **0.073** | 0.057 |
| unequal SD | normal (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.047 | 0.050 | 0.058 | 0.050 |
| unequal SD | normal (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.052 | 0.053 | **0.062** | 0.055 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.066** | 0.051 | **0.080** | 0.053 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 15 v 30 | 0 | **0.070** | 0.057 | **0.073** | **0.060** |
| unequal SD | normal (SD x2) | x smaller (1:2) | 30 v 60 | 0 | 0.053 | 0.048 | 0.053 | 0.050 |
| unequal SD | t3 (SD x2) | balanced | 8 v 8 | 0 | 0.054 | **0.067** | **0.100** | 0.047 |
| unequal SD | t3 (SD x2) | balanced | 15 v 15 | 0 | 0.051 | **0.070** | **0.087** | 0.046 |
| unequal SD | t3 (SD x2) | balanced | 30 v 30 | 0 | 0.051 | **0.064** | **0.074** | 0.050 |
| unequal SD | t3 (SD x2) | balanced | 50 v 50 | 0 | 0.051 | **0.065** | **0.071** | 0.050 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 16 v 8 | 0 | 0.042 | **0.064** | **0.088** | 0.052 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 30 v 15 | 0 | 0.049 | **0.077** | **0.088** | **0.060** |
| unequal SD | t3 (SD x2) | x larger (2:1) | 60 v 30 | 0 | 0.052 | **0.076** | **0.085** | 0.058 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 8 v 16 | 0 | **0.071** | **0.073** | **0.102** | 0.051 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 15 v 30 | 0 | 0.056 | 0.058 | **0.077** | 0.044 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 30 v 60 | 0 | 0.053 | 0.059 | **0.071** | 0.045 |
| unequal spread | likert_wide vs likert_mid | balanced | 8 v 8 | 0 | _0.036_ | _0.035_ | **0.078** | 0.053 |
| unequal spread | likert_wide vs likert_mid | balanced | 15 v 15 | 0 | 0.041 | 0.043 | 0.057 | 0.049 |
| unequal spread | likert_wide vs likert_mid | balanced | 30 v 30 | 0 | 0.050 | 0.049 | 0.058 | 0.055 |
| unequal spread | likert_wide vs likert_mid | balanced | 50 v 50 | 0 | 0.050 | 0.053 | **0.062** | 0.056 |
| unequal spread | likert_wide vs likert_mid | x larger (2:1) | 16 v 8 | 0 | 0.058 | 0.048 | **0.071** | 0.059 |
| unequal spread | likert_wide vs likert_mid | x larger (2:1) | 30 v 15 | 0 | 0.047 | 0.045 | 0.051 | 0.050 |
| unequal spread | likert_wide vs likert_mid | x larger (2:1) | 60 v 30 | 0 | 0.050 | 0.050 | 0.052 | 0.050 |
| unequal spread | likert_wide vs likert_mid | x smaller (1:2) | 8 v 16 | 0 | 0.048 | _0.027_ | 0.058 | 0.049 |
| unequal spread | likert_wide vs likert_mid | x smaller (1:2) | 15 v 30 | 0 | 0.052 | _0.036_ | 0.050 | 0.048 |
| unequal spread | likert_wide vs likert_mid | x smaller (1:2) | 30 v 60 | 0 | 0.053 | 0.050 | 0.055 | 0.053 |

## Block A, paired: two-sided Type I error at the true value

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| paired differences | contaminated | paired | 8 | 0 | 0.044 | **0.068** | **0.165** | _0.037_ |
| paired differences | contaminated | paired | 15 | 0 | _0.038_ | **0.074** | **0.130** | _0.034_ |
| paired differences | contaminated | paired | 30 | 0 | 0.044 | **0.092** | **0.125** | _0.032_ |
| paired differences | contaminated | paired | 50 | 0 | 0.050 | **0.093** | **0.117** | 0.044 |
| paired differences | diff_skew | paired | 8 | 0 | 0.048 | **0.065** | **0.116** | _0.036_ |
| paired differences | diff_skew | paired | 15 | 0 | 0.048 | 0.048 | **0.070** | 0.052 |
| paired differences | diff_skew | paired | 30 | 0 | 0.050 | 0.047 | 0.056 | 0.052 |
| paired differences | diff_skew | paired | 50 | 0 | 0.058 | 0.050 | 0.054 | 0.055 |
| paired differences | diff_sym | paired | 8 | 0 | _0.024_ | 0.050 | **0.111** | 0.058 |
| paired differences | diff_sym | paired | 15 | 0 | 0.057 | 0.052 | **0.090** | 0.057 |
| paired differences | diff_sym | paired | 30 | 0 | **0.061** | 0.058 | **0.070** | **0.061** |
| paired differences | diff_sym | paired | 50 | 0 | 0.048 | 0.044 | 0.056 | 0.045 |
| paired differences | exponential | paired | 8 | 0 | **0.113** | **0.064** | **0.139** | **0.109** |
| paired differences | exponential | paired | 15 | 0 | **0.091** | 0.048 | **0.087** | **0.085** |
| paired differences | exponential | paired | 30 | 0 | **0.073** | 0.054 | **0.080** | **0.072** |
| paired differences | exponential | paired | 50 | 0 | **0.079** | 0.058 | **0.075** | **0.075** |
| paired differences | lognormal | paired | 8 | 0 | **0.073** | 0.055 | **0.128** | **0.073** |
| paired differences | lognormal | paired | 15 | 0 | **0.078** | **0.061** | **0.102** | **0.073** |
| paired differences | lognormal | paired | 30 | 0 | **0.070** | 0.058 | **0.083** | **0.068** |
| paired differences | lognormal | paired | 50 | 0 | **0.063** | 0.054 | **0.068** | **0.060** |
| paired differences | normal | paired | 8 | 0 | 0.049 | 0.052 | **0.112** | 0.049 |
| paired differences | normal | paired | 15 | 0 | 0.047 | 0.048 | **0.082** | 0.047 |
| paired differences | normal | paired | 30 | 0 | 0.055 | 0.052 | **0.071** | 0.055 |
| paired differences | normal | paired | 50 | 0 | 0.051 | 0.048 | 0.059 | 0.051 |
| paired differences | t3 | paired | 8 | 0 | 0.043 | **0.066** | **0.142** | 0.041 |
| paired differences | t3 | paired | 15 | 0 | 0.049 | **0.080** | **0.121** | 0.041 |
| paired differences | t3 | paired | 30 | 0 | 0.054 | **0.084** | **0.109** | 0.046 |
| paired differences | t3 | paired | 50 | 0 | _0.039_ | 0.059 | **0.078** | _0.037_ |

## Block A, paired: one-sided error (TOST / minimal effect Type I error at a bound, worst direction)

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| paired differences | contaminated | paired | 8 | 0 | 0.050 | **0.077** | **0.119** | 0.050 |
| paired differences | contaminated | paired | 15 | 0 | 0.047 | **0.075** | **0.103** | 0.041 |
| paired differences | contaminated | paired | 30 | 0 | 0.054 | **0.087** | **0.106** | 0.050 |
| paired differences | contaminated | paired | 50 | 0 | 0.058 | **0.078** | **0.093** | 0.055 |
| paired differences | diff_skew | paired | 8 | 0 | 0.042 | **0.079** | **0.140** | **0.076** |
| paired differences | diff_skew | paired | 15 | 0 | **0.062** | 0.050 | **0.082** | **0.063** |
| paired differences | diff_skew | paired | 30 | 0 | **0.066** | 0.056 | **0.062** | **0.067** |
| paired differences | diff_skew | paired | 50 | 0 | 0.058 | 0.053 | **0.065** | **0.060** |
| paired differences | diff_sym | paired | 8 | 0 | 0.053 | 0.054 | **0.116** | 0.052 |
| paired differences | diff_sym | paired | 15 | 0 | 0.050 | 0.050 | **0.088** | 0.052 |
| paired differences | diff_sym | paired | 30 | 0 | 0.056 | 0.050 | **0.066** | 0.052 |
| paired differences | diff_sym | paired | 50 | 0 | 0.047 | 0.044 | **0.060** | 0.048 |
| paired differences | exponential | paired | 8 | 0 | **0.140** | **0.086** | **0.144** | **0.138** |
| paired differences | exponential | paired | 15 | 0 | **0.123** | **0.064** | **0.101** | **0.121** |
| paired differences | exponential | paired | 30 | 0 | **0.098** | 0.059 | **0.072** | **0.095** |
| paired differences | exponential | paired | 50 | 0 | **0.093** | **0.066** | **0.076** | **0.093** |
| paired differences | lognormal | paired | 8 | 0 | **0.110** | **0.075** | **0.128** | **0.112** |
| paired differences | lognormal | paired | 15 | 0 | **0.108** | **0.074** | **0.106** | **0.106** |
| paired differences | lognormal | paired | 30 | 0 | **0.100** | **0.068** | **0.084** | **0.098** |
| paired differences | lognormal | paired | 50 | 0 | **0.082** | **0.060** | **0.071** | **0.082** |
| paired differences | normal | paired | 8 | 0 | 0.049 | 0.051 | **0.086** | 0.050 |
| paired differences | normal | paired | 15 | 0 | 0.058 | 0.050 | **0.072** | 0.056 |
| paired differences | normal | paired | 30 | 0 | 0.052 | 0.052 | **0.062** | 0.053 |
| paired differences | normal | paired | 50 | 0 | 0.050 | 0.048 | 0.053 | 0.050 |
| paired differences | t3 | paired | 8 | 0 | 0.044 | **0.066** | **0.107** | 0.044 |
| paired differences | t3 | paired | 15 | 0 | 0.055 | **0.070** | **0.098** | 0.053 |
| paired differences | t3 | paired | 30 | 0 | 0.051 | **0.071** | **0.087** | 0.051 |
| paired differences | t3 | paired | 50 | 0 | 0.043 | 0.059 | **0.067** | 0.042 |

## Block B, two-sample: two-sided Type I error at the true value

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| identical | contaminated | 1:2 | 8 v 16 | 0.2 | 0.043 | 0.049 | 0.047 | 0.043 |
| identical | contaminated | 1:2 | 15 v 30 | 0.2 | 0.052 | 0.055 | 0.048 | 0.053 |
| identical | contaminated | 1:2 | 30 v 60 | 0.2 | 0.051 | 0.048 | 0.049 | 0.052 |
| identical | contaminated | balanced | 8 v 8 | 0.2 | 0.058 | 0.059 | 0.054 | 0.045 |
| identical | contaminated | balanced | 15 v 15 | 0.2 | 0.050 | 0.050 | 0.043 | 0.049 |
| identical | contaminated | balanced | 30 v 30 | 0.2 | 0.050 | 0.052 | 0.050 | 0.049 |
| identical | contaminated | balanced | 50 v 50 | 0.2 | 0.044 | 0.045 | 0.042 | 0.045 |
| identical | exponential | 1:2 | 8 v 16 | 0.2 | 0.058 | 0.054 | **0.080** | 0.046 |
| identical | exponential | 1:2 | 15 v 30 | 0.2 | 0.049 | 0.053 | 0.053 | 0.048 |
| identical | exponential | 1:2 | 30 v 60 | 0.2 | 0.055 | 0.057 | **0.060** | 0.055 |
| identical | exponential | balanced | 8 v 8 | 0.2 | 0.049 | _0.036_ | **0.070** | _0.029_ |
| identical | exponential | balanced | 15 v 15 | 0.2 | 0.052 | 0.047 | 0.058 | 0.044 |
| identical | exponential | balanced | 30 v 30 | 0.2 | 0.057 | **0.066** | **0.060** | 0.054 |
| identical | exponential | balanced | 50 v 50 | 0.2 | 0.053 | 0.058 | 0.051 | 0.048 |
| identical | lognormal | 1:2 | 8 v 16 | 0.2 | 0.046 | 0.045 | **0.067** | 0.042 |
| identical | lognormal | 1:2 | 15 v 30 | 0.2 | 0.047 | 0.048 | 0.051 | 0.049 |
| identical | lognormal | 1:2 | 30 v 60 | 0.2 | 0.050 | 0.052 | 0.053 | 0.050 |
| identical | lognormal | balanced | 8 v 8 | 0.2 | 0.052 | 0.049 | **0.065** | _0.039_ |
| identical | lognormal | balanced | 15 v 15 | 0.2 | 0.054 | 0.051 | 0.058 | 0.051 |
| identical | lognormal | balanced | 30 v 30 | 0.2 | 0.042 | 0.050 | 0.045 | _0.040_ |
| identical | lognormal | balanced | 50 v 50 | 0.2 | 0.048 | 0.046 | 0.050 | 0.044 |
| identical | normal | 1:2 | 8 v 16 | 0.2 | 0.051 | 0.047 | **0.060** | 0.051 |
| identical | normal | 1:2 | 15 v 30 | 0.2 | _0.039_ | 0.044 | 0.046 | 0.044 |
| identical | normal | 1:2 | 30 v 60 | 0.2 | 0.054 | 0.051 | 0.057 | 0.052 |
| identical | normal | balanced | 8 v 8 | 0.2 | 0.054 | 0.048 | **0.064** | 0.049 |
| identical | normal | balanced | 15 v 15 | 0.2 | 0.045 | 0.044 | 0.052 | 0.048 |
| identical | normal | balanced | 30 v 30 | 0.2 | 0.052 | 0.046 | 0.053 | 0.053 |
| identical | normal | balanced | 50 v 50 | 0.2 | 0.045 | 0.043 | 0.048 | 0.044 |
| identical | t3 | 1:2 | 8 v 16 | 0.2 | 0.050 | 0.055 | 0.058 | 0.044 |
| identical | t3 | 1:2 | 15 v 30 | 0.2 | 0.049 | 0.052 | 0.046 | 0.048 |
| identical | t3 | 1:2 | 30 v 60 | 0.2 | 0.048 | 0.051 | 0.048 | 0.047 |
| identical | t3 | balanced | 8 v 8 | 0.2 | 0.046 | 0.049 | 0.050 | _0.039_ |
| identical | t3 | balanced | 15 v 15 | 0.2 | 0.050 | 0.056 | 0.046 | 0.046 |
| identical | t3 | balanced | 30 v 30 | 0.2 | 0.047 | 0.048 | 0.050 | 0.046 |
| identical | t3 | balanced | 50 v 50 | 0.2 | 0.050 | 0.050 | 0.052 | 0.050 |
| unequal SD | contaminated (SD x2) | balanced | 8 v 8 | 0.2 | **0.061** | **0.060** | 0.058 | 0.046 |
| unequal SD | contaminated (SD x2) | balanced | 15 v 15 | 0.2 | 0.051 | 0.050 | 0.044 | 0.047 |
| unequal SD | contaminated (SD x2) | balanced | 30 v 30 | 0.2 | 0.052 | 0.048 | 0.051 | 0.047 |
| unequal SD | contaminated (SD x2) | balanced | 50 v 50 | 0.2 | 0.049 | 0.047 | 0.047 | 0.047 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.047 | 0.055 | 0.049 | 0.044 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.046 | 0.049 | 0.046 | 0.046 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.047 | 0.046 | 0.046 | 0.048 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.062** | 0.053 | **0.060** | 0.049 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.061** | 0.054 | 0.049 | 0.051 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | **0.064** | 0.059 | 0.058 | 0.058 |
| unequal SD | exponential (SD x2) | balanced | 8 v 8 | 0.2 | **0.062** | 0.055 | **0.072** | 0.049 |
| unequal SD | exponential (SD x2) | balanced | 15 v 15 | 0.2 | 0.057 | 0.050 | 0.051 | 0.051 |
| unequal SD | exponential (SD x2) | balanced | 30 v 30 | 0.2 | 0.054 | 0.050 | 0.047 | 0.049 |
| unequal SD | exponential (SD x2) | balanced | 50 v 50 | 0.2 | 0.058 | 0.055 | 0.055 | 0.057 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.055 | 0.058 | **0.068** | 0.048 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.043 | 0.050 | 0.053 | 0.041 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.047 | 0.052 | 0.049 | 0.048 |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.078** | **0.068** | **0.080** | **0.065** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.074** | **0.062** | **0.061** | 0.059 |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | **0.067** | 0.052 | 0.058 | 0.058 |
| unequal SD | lognormal (SD x2) | balanced | 8 v 8 | 0.2 | 0.058 | 0.055 | **0.073** | 0.050 |
| unequal SD | lognormal (SD x2) | balanced | 15 v 15 | 0.2 | **0.068** | **0.062** | 0.058 | **0.061** |
| unequal SD | lognormal (SD x2) | balanced | 30 v 30 | 0.2 | 0.053 | 0.050 | 0.052 | 0.048 |
| unequal SD | lognormal (SD x2) | balanced | 50 v 50 | 0.2 | 0.055 | 0.053 | 0.056 | 0.051 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.045 | 0.047 | 0.055 | 0.042 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | _0.036_ | 0.041 | 0.046 | 0.041 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.047 | 0.056 | 0.054 | 0.051 |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.073** | **0.060** | **0.079** | 0.058 |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.072** | **0.064** | 0.053 | **0.063** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | 0.053 | 0.049 | 0.053 | 0.051 |
| unequal SD | normal (SD x2) | balanced | 8 v 8 | 0.2 | **0.060** | 0.043 | **0.060** | 0.051 |
| unequal SD | normal (SD x2) | balanced | 15 v 15 | 0.2 | 0.058 | 0.047 | 0.052 | 0.052 |
| unequal SD | normal (SD x2) | balanced | 30 v 30 | 0.2 | 0.045 | _0.038_ | 0.048 | 0.043 |
| unequal SD | normal (SD x2) | balanced | 50 v 50 | 0.2 | 0.056 | 0.052 | 0.052 | 0.052 |
| unequal SD | normal (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.042 | 0.041 | 0.057 | 0.044 |
| unequal SD | normal (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.040 | 0.044 | 0.049 | 0.043 |
| unequal SD | normal (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.044 | 0.043 | 0.049 | 0.046 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.071** | 0.056 | **0.080** | 0.055 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.064** | 0.050 | **0.062** | 0.056 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | **0.061** | 0.049 | **0.063** | 0.057 |
| unequal SD | t3 (SD x2) | balanced | 8 v 8 | 0.2 | 0.056 | 0.053 | 0.056 | 0.045 |
| unequal SD | t3 (SD x2) | balanced | 15 v 15 | 0.2 | **0.060** | 0.056 | 0.049 | 0.051 |
| unequal SD | t3 (SD x2) | balanced | 30 v 30 | 0.2 | 0.054 | 0.051 | 0.049 | 0.049 |
| unequal SD | t3 (SD x2) | balanced | 50 v 50 | 0.2 | 0.051 | 0.049 | 0.052 | 0.049 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.052 | 0.055 | 0.048 | 0.049 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.048 | 0.054 | 0.048 | 0.049 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.053 | 0.056 | 0.052 | 0.053 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.078** | **0.066** | **0.072** | 0.055 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.062** | **0.060** | 0.048 | 0.049 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | 0.054 | 0.048 | 0.046 | 0.046 |

## Block B, two-sample: one-sided error (TOST / minimal effect Type I error at a bound, worst direction)

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| identical | contaminated | 1:2 | 8 v 16 | 0.2 | 0.046 | 0.050 | 0.048 | 0.048 |
| identical | contaminated | 1:2 | 15 v 30 | 0.2 | 0.056 | 0.058 | 0.056 | 0.056 |
| identical | contaminated | 1:2 | 30 v 60 | 0.2 | 0.051 | 0.050 | 0.050 | 0.050 |
| identical | contaminated | balanced | 8 v 8 | 0.2 | 0.051 | 0.053 | 0.051 | 0.046 |
| identical | contaminated | balanced | 15 v 15 | 0.2 | 0.053 | 0.051 | 0.042 | 0.050 |
| identical | contaminated | balanced | 30 v 30 | 0.2 | 0.056 | 0.056 | 0.056 | 0.056 |
| identical | contaminated | balanced | 50 v 50 | 0.2 | 0.051 | 0.050 | 0.053 | 0.050 |
| identical | exponential | 1:2 | 8 v 16 | 0.2 | 0.050 | **0.060** | **0.069** | 0.052 |
| identical | exponential | 1:2 | 15 v 30 | 0.2 | 0.053 | 0.056 | 0.054 | 0.052 |
| identical | exponential | 1:2 | 30 v 60 | 0.2 | 0.054 | 0.058 | 0.058 | 0.059 |
| identical | exponential | balanced | 8 v 8 | 0.2 | 0.048 | 0.047 | **0.064** | _0.039_ |
| identical | exponential | balanced | 15 v 15 | 0.2 | 0.056 | 0.058 | 0.058 | 0.046 |
| identical | exponential | balanced | 30 v 30 | 0.2 | 0.050 | 0.053 | 0.053 | 0.049 |
| identical | exponential | balanced | 50 v 50 | 0.2 | 0.052 | 0.056 | 0.054 | 0.052 |
| identical | lognormal | 1:2 | 8 v 16 | 0.2 | **0.060** | 0.059 | **0.070** | 0.051 |
| identical | lognormal | 1:2 | 15 v 30 | 0.2 | 0.045 | 0.051 | 0.047 | 0.048 |
| identical | lognormal | 1:2 | 30 v 60 | 0.2 | 0.056 | 0.055 | 0.056 | 0.057 |
| identical | lognormal | balanced | 8 v 8 | 0.2 | 0.050 | 0.046 | **0.064** | 0.042 |
| identical | lognormal | balanced | 15 v 15 | 0.2 | 0.059 | 0.059 | **0.061** | 0.054 |
| identical | lognormal | balanced | 30 v 30 | 0.2 | 0.045 | 0.046 | 0.048 | 0.044 |
| identical | lognormal | balanced | 50 v 50 | 0.2 | 0.051 | 0.054 | 0.052 | 0.051 |
| identical | normal | 1:2 | 8 v 16 | 0.2 | 0.047 | 0.047 | 0.053 | 0.048 |
| identical | normal | 1:2 | 15 v 30 | 0.2 | 0.049 | 0.047 | 0.056 | 0.048 |
| identical | normal | 1:2 | 30 v 60 | 0.2 | 0.049 | 0.047 | 0.053 | 0.050 |
| identical | normal | balanced | 8 v 8 | 0.2 | 0.053 | 0.048 | **0.061** | 0.052 |
| identical | normal | balanced | 15 v 15 | 0.2 | 0.049 | 0.046 | 0.051 | 0.048 |
| identical | normal | balanced | 30 v 30 | 0.2 | 0.053 | 0.056 | 0.057 | 0.054 |
| identical | normal | balanced | 50 v 50 | 0.2 | 0.051 | 0.050 | 0.053 | 0.050 |
| identical | t3 | 1:2 | 8 v 16 | 0.2 | 0.055 | 0.057 | 0.049 | 0.050 |
| identical | t3 | 1:2 | 15 v 30 | 0.2 | 0.050 | 0.051 | 0.048 | 0.049 |
| identical | t3 | 1:2 | 30 v 60 | 0.2 | 0.054 | 0.052 | 0.053 | 0.054 |
| identical | t3 | balanced | 8 v 8 | 0.2 | 0.050 | 0.052 | 0.047 | 0.043 |
| identical | t3 | balanced | 15 v 15 | 0.2 | 0.058 | 0.059 | 0.052 | 0.059 |
| identical | t3 | balanced | 30 v 30 | 0.2 | 0.051 | 0.051 | 0.048 | 0.050 |
| identical | t3 | balanced | 50 v 50 | 0.2 | 0.055 | 0.053 | 0.054 | 0.053 |
| unequal SD | contaminated (SD x2) | balanced | 8 v 8 | 0.2 | 0.059 | **0.062** | 0.050 | 0.051 |
| unequal SD | contaminated (SD x2) | balanced | 15 v 15 | 0.2 | 0.057 | 0.053 | 0.048 | 0.052 |
| unequal SD | contaminated (SD x2) | balanced | 30 v 30 | 0.2 | 0.052 | 0.050 | 0.048 | 0.050 |
| unequal SD | contaminated (SD x2) | balanced | 50 v 50 | 0.2 | 0.058 | 0.053 | 0.059 | 0.056 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.050 | 0.051 | 0.044 | 0.046 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.050 | 0.055 | 0.049 | 0.052 |
| unequal SD | contaminated (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.045 | 0.050 | 0.048 | 0.050 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.065** | 0.058 | 0.058 | 0.053 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.061** | 0.053 | 0.053 | 0.053 |
| unequal SD | contaminated (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | **0.060** | 0.053 | 0.059 | 0.056 |
| unequal SD | exponential (SD x2) | balanced | 8 v 8 | 0.2 | **0.073** | **0.076** | **0.068** | **0.069** |
| unequal SD | exponential (SD x2) | balanced | 15 v 15 | 0.2 | **0.068** | **0.060** | 0.054 | **0.064** |
| unequal SD | exponential (SD x2) | balanced | 30 v 30 | 0.2 | 0.059 | 0.050 | 0.048 | 0.053 |
| unequal SD | exponential (SD x2) | balanced | 50 v 50 | 0.2 | **0.066** | 0.058 | 0.058 | **0.066** |
| unequal SD | exponential (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | **0.068** | **0.065** | **0.077** | **0.065** |
| unequal SD | exponential (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.055 | 0.051 | 0.051 | 0.050 |
| unequal SD | exponential (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | **0.060** | 0.053 | 0.055 | 0.056 |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.078** | **0.066** | **0.073** | **0.075** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.075** | **0.068** | **0.065** | **0.077** |
| unequal SD | exponential (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | **0.070** | 0.055 | 0.056 | **0.070** |
| unequal SD | lognormal (SD x2) | balanced | 8 v 8 | 0.2 | **0.066** | **0.060** | **0.061** | 0.059 |
| unequal SD | lognormal (SD x2) | balanced | 15 v 15 | 0.2 | **0.072** | **0.068** | **0.062** | **0.073** |
| unequal SD | lognormal (SD x2) | balanced | 30 v 30 | 0.2 | **0.066** | 0.053 | 0.052 | **0.062** |
| unequal SD | lognormal (SD x2) | balanced | 50 v 50 | 0.2 | 0.058 | 0.052 | 0.050 | 0.058 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.050 | 0.053 | 0.058 | 0.052 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.048 | 0.052 | 0.054 | 0.048 |
| unequal SD | lognormal (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.057 | 0.056 | 0.056 | 0.055 |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.074** | **0.064** | **0.063** | **0.070** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.076** | **0.064** | **0.068** | **0.075** |
| unequal SD | lognormal (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | 0.053 | 0.048 | 0.050 | 0.054 |
| unequal SD | normal (SD x2) | balanced | 8 v 8 | 0.2 | 0.055 | 0.050 | **0.064** | 0.050 |
| unequal SD | normal (SD x2) | balanced | 15 v 15 | 0.2 | 0.052 | 0.047 | 0.055 | 0.048 |
| unequal SD | normal (SD x2) | balanced | 30 v 30 | 0.2 | 0.051 | 0.050 | 0.051 | 0.050 |
| unequal SD | normal (SD x2) | balanced | 50 v 50 | 0.2 | 0.050 | 0.050 | 0.052 | 0.050 |
| unequal SD | normal (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.050 | 0.050 | 0.058 | 0.054 |
| unequal SD | normal (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.042 | 0.046 | 0.050 | 0.047 |
| unequal SD | normal (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.044 | 0.048 | 0.052 | 0.051 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.068** | 0.053 | **0.066** | 0.056 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | **0.060** | 0.053 | 0.056 | 0.056 |
| unequal SD | normal (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | **0.068** | 0.058 | **0.064** | **0.061** |
| unequal SD | t3 (SD x2) | balanced | 8 v 8 | 0.2 | **0.060** | **0.060** | **0.062** | 0.051 |
| unequal SD | t3 (SD x2) | balanced | 15 v 15 | 0.2 | 0.054 | 0.051 | 0.048 | 0.050 |
| unequal SD | t3 (SD x2) | balanced | 30 v 30 | 0.2 | 0.051 | 0.050 | 0.050 | 0.048 |
| unequal SD | t3 (SD x2) | balanced | 50 v 50 | 0.2 | 0.048 | 0.046 | 0.048 | 0.048 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 16 v 8 | 0.2 | 0.052 | 0.056 | 0.054 | 0.053 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 30 v 15 | 0.2 | 0.050 | 0.057 | 0.053 | 0.054 |
| unequal SD | t3 (SD x2) | x larger (2:1) | 60 v 30 | 0.2 | 0.049 | 0.054 | 0.056 | 0.053 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 8 v 16 | 0.2 | **0.073** | **0.065** | **0.063** | 0.059 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 15 v 30 | 0.2 | 0.058 | 0.052 | 0.050 | 0.053 |
| unequal SD | t3 (SD x2) | x smaller (1:2) | 30 v 60 | 0.2 | 0.058 | 0.050 | 0.051 | 0.052 |

## Block B, paired: two-sided Type I error at the true value

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| paired differences | contaminated | paired | 8 | 0.2 | **0.061** | 0.057 | **0.085** | 0.046 |
| paired differences | contaminated | paired | 15 | 0.2 | 0.056 | 0.047 | 0.056 | 0.048 |
| paired differences | contaminated | paired | 30 | 0.2 | 0.050 | 0.046 | 0.049 | 0.046 |
| paired differences | contaminated | paired | 50 | 0.2 | 0.055 | 0.055 | 0.057 | 0.056 |
| paired differences | exponential | paired | 8 | 0.2 | 0.057 | 0.043 | **0.077** | **0.062** |
| paired differences | exponential | paired | 15 | 0.2 | **0.068** | 0.054 | **0.068** | **0.076** |
| paired differences | exponential | paired | 30 | 0.2 | 0.054 | 0.040 | 0.051 | **0.060** |
| paired differences | exponential | paired | 50 | 0.2 | 0.056 | 0.050 | 0.059 | 0.059 |
| paired differences | lognormal | paired | 8 | 0.2 | 0.059 | 0.046 | **0.080** | **0.060** |
| paired differences | lognormal | paired | 15 | 0.2 | **0.067** | 0.056 | **0.066** | **0.069** |
| paired differences | lognormal | paired | 30 | 0.2 | 0.052 | 0.041 | 0.053 | 0.052 |
| paired differences | lognormal | paired | 50 | 0.2 | 0.055 | 0.051 | 0.056 | 0.058 |
| paired differences | normal | paired | 8 | 0.2 | **0.068** | 0.057 | **0.095** | **0.068** |
| paired differences | normal | paired | 15 | 0.2 | **0.068** | 0.045 | 0.055 | **0.061** |
| paired differences | normal | paired | 30 | 0.2 | 0.057 | 0.053 | 0.058 | 0.056 |
| paired differences | normal | paired | 50 | 0.2 | 0.053 | 0.045 | 0.054 | 0.053 |
| paired differences | t3 | paired | 8 | 0.2 | **0.061** | 0.051 | **0.088** | 0.047 |
| paired differences | t3 | paired | 15 | 0.2 | **0.064** | 0.051 | **0.066** | 0.048 |
| paired differences | t3 | paired | 30 | 0.2 | 0.053 | 0.052 | 0.055 | 0.053 |
| paired differences | t3 | paired | 50 | 0.2 | 0.050 | 0.053 | 0.056 | 0.048 |

## Block B, paired: one-sided error (TOST / minimal effect Type I error at a bound, worst direction)

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| paired differences | contaminated | paired | 8 | 0.2 | 0.054 | **0.060** | **0.066** | 0.052 |
| paired differences | contaminated | paired | 15 | 0.2 | 0.054 | 0.047 | 0.054 | 0.052 |
| paired differences | contaminated | paired | 30 | 0.2 | 0.050 | 0.048 | 0.052 | 0.050 |
| paired differences | contaminated | paired | 50 | 0.2 | 0.051 | 0.052 | 0.053 | 0.055 |
| paired differences | exponential | paired | 8 | 0.2 | **0.066** | 0.048 | **0.076** | **0.075** |
| paired differences | exponential | paired | 15 | 0.2 | **0.082** | **0.061** | **0.070** | **0.090** |
| paired differences | exponential | paired | 30 | 0.2 | **0.075** | 0.050 | 0.056 | **0.078** |
| paired differences | exponential | paired | 50 | 0.2 | **0.068** | 0.050 | 0.054 | **0.070** |
| paired differences | lognormal | paired | 8 | 0.2 | **0.068** | 0.058 | **0.067** | **0.076** |
| paired differences | lognormal | paired | 15 | 0.2 | **0.073** | 0.056 | **0.061** | **0.076** |
| paired differences | lognormal | paired | 30 | 0.2 | **0.061** | 0.048 | 0.052 | **0.064** |
| paired differences | lognormal | paired | 50 | 0.2 | **0.060** | 0.054 | 0.058 | **0.060** |
| paired differences | normal | paired | 8 | 0.2 | 0.058 | 0.056 | **0.080** | **0.064** |
| paired differences | normal | paired | 15 | 0.2 | 0.053 | 0.042 | 0.059 | 0.051 |
| paired differences | normal | paired | 30 | 0.2 | 0.058 | 0.054 | 0.058 | 0.056 |
| paired differences | normal | paired | 50 | 0.2 | 0.059 | 0.058 | 0.059 | 0.058 |
| paired differences | t3 | paired | 8 | 0.2 | 0.057 | 0.058 | **0.077** | 0.053 |
| paired differences | t3 | paired | 15 | 0.2 | 0.056 | 0.051 | 0.054 | 0.051 |
| paired differences | t3 | paired | 30 | 0.2 | 0.058 | 0.058 | 0.059 | 0.056 |
| paired differences | t3 | paired | 50 | 0.2 | 0.052 | 0.051 | 0.051 | 0.050 |

## Block C, two-sample: two-sided Type I error at the true value

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| different shape | exponential vs normal | balanced | 8 v 8 | 0 | **0.060** | 0.048 | **0.103** | 0.053 |
| different shape | exponential vs normal | balanced | 8 v 8 | 0.2 | 0.051 | 0.048 | **0.064** | 0.050 |
| different shape | exponential vs normal | balanced | 15 v 15 | 0 | 0.047 | 0.052 | **0.071** | 0.043 |
| different shape | exponential vs normal | balanced | 15 v 15 | 0.2 | 0.045 | _0.039_ | 0.053 | 0.046 |
| different shape | exponential vs normal | balanced | 30 v 30 | 0 | 0.055 | 0.051 | **0.067** | 0.053 |
| different shape | exponential vs normal | balanced | 30 v 30 | 0.2 | 0.051 | 0.052 | 0.053 | 0.054 |
| different shape | exponential vs normal | balanced | 50 v 50 | 0 | 0.049 | 0.051 | 0.059 | 0.052 |
| different shape | exponential vs normal | balanced | 50 v 50 | 0.2 | 0.055 | 0.050 | 0.054 | 0.055 |
| different shape | exponential vs normal | x larger (2:1) | 16 v 8 | 0 | 0.055 | 0.051 | **0.096** | 0.049 |
| different shape | exponential vs normal | x larger (2:1) | 16 v 8 | 0.2 | 0.059 | **0.061** | **0.068** | 0.053 |
| different shape | exponential vs normal | x larger (2:1) | 30 v 15 | 0 | 0.056 | **0.060** | **0.080** | 0.054 |
| different shape | exponential vs normal | x larger (2:1) | 30 v 15 | 0.2 | **0.065** | **0.062** | **0.070** | **0.067** |
| different shape | exponential vs normal | x larger (2:1) | 60 v 30 | 0 | 0.052 | 0.055 | **0.071** | 0.053 |
| different shape | exponential vs normal | x larger (2:1) | 60 v 30 | 0.2 | 0.058 | 0.055 | **0.064** | 0.059 |
| different shape | exponential vs normal | x smaller (1:2) | 8 v 16 | 0 | **0.061** | 0.055 | **0.102** | **0.063** |
| different shape | exponential vs normal | x smaller (1:2) | 8 v 16 | 0.2 | 0.048 | 0.050 | **0.066** | 0.050 |
| different shape | exponential vs normal | x smaller (1:2) | 15 v 30 | 0 | 0.059 | **0.064** | **0.085** | **0.061** |
| different shape | exponential vs normal | x smaller (1:2) | 15 v 30 | 0.2 | 0.041 | 0.045 | 0.056 | 0.045 |
| different shape | exponential vs normal | x smaller (1:2) | 30 v 60 | 0 | 0.048 | 0.049 | **0.062** | 0.049 |
| different shape | exponential vs normal | x smaller (1:2) | 30 v 60 | 0.2 | 0.043 | 0.046 | 0.047 | 0.046 |
| different shape | likert_skew vs likert_mid | balanced | 8 v 8 | 0 | 0.057 | 0.046 | **0.100** | 0.053 |
| different shape | likert_skew vs likert_mid | balanced | 15 v 15 | 0 | 0.057 | 0.048 | **0.079** | 0.055 |
| different shape | likert_skew vs likert_mid | balanced | 30 v 30 | 0 | 0.058 | 0.052 | **0.064** | 0.055 |
| different shape | likert_skew vs likert_mid | balanced | 50 v 50 | 0 | 0.041 | _0.038_ | 0.044 | _0.039_ |
| different shape | likert_skew vs likert_mid | x larger (2:1) | 16 v 8 | 0 | 0.046 | 0.046 | **0.085** | 0.050 |
| different shape | likert_skew vs likert_mid | x larger (2:1) | 30 v 15 | 0 | 0.051 | 0.050 | **0.073** | 0.054 |
| different shape | likert_skew vs likert_mid | x larger (2:1) | 60 v 30 | 0 | 0.048 | 0.045 | 0.058 | 0.048 |
| different shape | likert_skew vs likert_mid | x smaller (1:2) | 8 v 16 | 0 | 0.059 | 0.053 | **0.087** | **0.060** |
| different shape | likert_skew vs likert_mid | x smaller (1:2) | 15 v 30 | 0 | 0.050 | 0.046 | 0.059 | 0.052 |
| different shape | likert_skew vs likert_mid | x smaller (1:2) | 30 v 60 | 0 | 0.050 | 0.048 | 0.054 | 0.050 |

## Block C, two-sample: one-sided error (TOST / minimal effect Type I error at a bound, worst direction)

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| different shape | exponential vs normal | balanced | 8 v 8 | 0 | **0.075** | **0.065** | **0.096** | **0.071** |
| different shape | exponential vs normal | balanced | 8 v 8 | 0.2 | 0.058 | 0.059 | **0.074** | 0.048 |
| different shape | exponential vs normal | balanced | 15 v 15 | 0 | 0.058 | 0.053 | **0.066** | 0.058 |
| different shape | exponential vs normal | balanced | 15 v 15 | 0.2 | 0.054 | 0.052 | 0.056 | 0.051 |
| different shape | exponential vs normal | balanced | 30 v 30 | 0 | 0.059 | 0.055 | **0.061** | **0.060** |
| different shape | exponential vs normal | balanced | 30 v 30 | 0.2 | **0.061** | 0.053 | 0.059 | 0.059 |
| different shape | exponential vs normal | balanced | 50 v 50 | 0 | **0.061** | 0.051 | 0.058 | **0.060** |
| different shape | exponential vs normal | balanced | 50 v 50 | 0.2 | 0.056 | 0.054 | 0.058 | 0.053 |
| different shape | exponential vs normal | x larger (2:1) | 16 v 8 | 0 | **0.076** | 0.053 | **0.076** | 0.058 |
| different shape | exponential vs normal | x larger (2:1) | 16 v 8 | 0.2 | **0.063** | **0.062** | **0.066** | **0.061** |
| different shape | exponential vs normal | x larger (2:1) | 30 v 15 | 0 | **0.074** | 0.058 | **0.073** | **0.066** |
| different shape | exponential vs normal | x larger (2:1) | 30 v 15 | 0.2 | **0.062** | **0.060** | 0.059 | **0.060** |
| different shape | exponential vs normal | x larger (2:1) | 60 v 30 | 0 | **0.075** | 0.059 | **0.062** | **0.063** |
| different shape | exponential vs normal | x larger (2:1) | 60 v 30 | 0.2 | 0.059 | **0.060** | **0.064** | **0.061** |
| different shape | exponential vs normal | x smaller (1:2) | 8 v 16 | 0 | **0.082** | **0.087** | **0.110** | **0.102** |
| different shape | exponential vs normal | x smaller (1:2) | 8 v 16 | 0.2 | 0.048 | 0.049 | **0.072** | 0.051 |
| different shape | exponential vs normal | x smaller (1:2) | 15 v 30 | 0 | **0.077** | **0.068** | **0.080** | **0.084** |
| different shape | exponential vs normal | x smaller (1:2) | 15 v 30 | 0.2 | 0.048 | 0.056 | 0.053 | 0.058 |
| different shape | exponential vs normal | x smaller (1:2) | 30 v 60 | 0 | **0.060** | 0.053 | **0.060** | **0.069** |
| different shape | exponential vs normal | x smaller (1:2) | 30 v 60 | 0.2 | 0.048 | 0.052 | 0.056 | 0.056 |
| different shape | likert_skew vs likert_mid | balanced | 8 v 8 | 0 | **0.079** | 0.058 | **0.118** | **0.066** |
| different shape | likert_skew vs likert_mid | balanced | 15 v 15 | 0 | 0.056 | 0.053 | **0.069** | 0.059 |
| different shape | likert_skew vs likert_mid | balanced | 30 v 30 | 0 | **0.066** | 0.050 | **0.060** | 0.058 |
| different shape | likert_skew vs likert_mid | balanced | 50 v 50 | 0 | **0.061** | 0.047 | 0.053 | 0.054 |
| different shape | likert_skew vs likert_mid | x larger (2:1) | 16 v 8 | 0 | **0.068** | 0.051 | **0.093** | **0.061** |
| different shape | likert_skew vs likert_mid | x larger (2:1) | 30 v 15 | 0 | 0.058 | 0.051 | **0.063** | 0.054 |
| different shape | likert_skew vs likert_mid | x larger (2:1) | 60 v 30 | 0 | 0.053 | 0.045 | 0.052 | 0.046 |
| different shape | likert_skew vs likert_mid | x smaller (1:2) | 8 v 16 | 0 | **0.066** | **0.060** | **0.088** | **0.071** |
| different shape | likert_skew vs likert_mid | x smaller (1:2) | 15 v 30 | 0 | 0.057 | 0.048 | **0.064** | 0.059 |
| different shape | likert_skew vs likert_mid | x smaller (1:2) | 30 v 60 | 0 | 0.056 | 0.048 | 0.053 | 0.056 |

## Block C, paired: two-sided Type I error at the true value

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| paired differences | contaminated | paired | 15 | 0 | 0.041 | **0.068** | **0.127** | _0.033_ |
| paired differences | contaminated | paired | 15 | 0.2 | 0.058 | 0.044 | 0.052 | 0.051 |
| paired differences | contaminated | paired | 30 | 0 | 0.050 | **0.094** | **0.128** | 0.042 |
| paired differences | contaminated | paired | 30 | 0.2 | 0.056 | 0.051 | 0.059 | 0.052 |
| paired differences | exponential | paired | 15 | 0 | **0.103** | **0.064** | **0.109** | **0.096** |
| paired differences | exponential | paired | 15 | 0.2 | **0.066** | 0.049 | 0.059 | **0.068** |
| paired differences | exponential | paired | 30 | 0 | **0.078** | 0.056 | **0.078** | **0.075** |
| paired differences | exponential | paired | 30 | 0.2 | 0.055 | 0.044 | 0.055 | **0.063** |
| paired differences | normal | paired | 15 | 0 | 0.054 | 0.052 | **0.094** | 0.054 |
| paired differences | normal | paired | 15 | 0.2 | 0.056 | 0.046 | **0.064** | 0.054 |
| paired differences | normal | paired | 30 | 0 | _0.040_ | 0.044 | 0.057 | _0.039_ |
| paired differences | normal | paired | 30 | 0.2 | 0.053 | 0.041 | 0.049 | 0.052 |
| paired differences | t3 | paired | 15 | 0 | 0.042 | **0.072** | **0.113** | _0.037_ |
| paired differences | t3 | paired | 15 | 0.2 | **0.065** | 0.054 | 0.058 | 0.053 |
| paired differences | t3 | paired | 30 | 0 | 0.053 | **0.084** | **0.110** | 0.046 |
| paired differences | t3 | paired | 30 | 0.2 | 0.047 | 0.049 | 0.052 | 0.046 |

## Block C, paired: one-sided error (TOST / minimal effect Type I error at a bound, worst direction)

| structure | scenario | layout | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| paired differences | contaminated | paired | 15 | 0 | 0.043 | **0.076** | **0.101** | _0.039_ |
| paired differences | contaminated | paired | 15 | 0.2 | 0.056 | 0.049 | 0.052 | 0.052 |
| paired differences | contaminated | paired | 30 | 0 | 0.054 | **0.092** | **0.103** | 0.050 |
| paired differences | contaminated | paired | 30 | 0.2 | 0.058 | 0.055 | 0.059 | 0.057 |
| paired differences | exponential | paired | 15 | 0 | **0.130** | **0.080** | **0.107** | **0.126** |
| paired differences | exponential | paired | 15 | 0.2 | **0.077** | 0.052 | 0.058 | **0.082** |
| paired differences | exponential | paired | 30 | 0 | **0.096** | **0.064** | **0.082** | **0.098** |
| paired differences | exponential | paired | 30 | 0.2 | **0.075** | 0.053 | 0.058 | **0.078** |
| paired differences | normal | paired | 15 | 0 | **0.063** | **0.062** | **0.085** | **0.061** |
| paired differences | normal | paired | 15 | 0.2 | 0.054 | 0.048 | 0.057 | 0.054 |
| paired differences | normal | paired | 30 | 0 | 0.051 | 0.050 | 0.058 | 0.050 |
| paired differences | normal | paired | 30 | 0.2 | 0.051 | 0.050 | 0.051 | 0.053 |
| paired differences | t3 | paired | 15 | 0 | 0.046 | **0.065** | **0.088** | 0.045 |
| paired differences | t3 | paired | 15 | 0.2 | 0.058 | 0.054 | **0.064** | 0.054 |
| paired differences | t3 | paired | 30 | 0 | 0.059 | **0.074** | **0.086** | 0.055 |
| paired differences | t3 | paired | 30 | 0.2 | 0.050 | 0.048 | 0.057 | 0.048 |

## Block C: power of the standard test (true shift 0.5)

| design | scenario | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|
| two | contaminated | 15 v 15 | 0 | 0.406 | 0.431 | 0.480 | 0.362 |
| two | contaminated | 30 v 30 | 0 | 0.546 | 0.531 | 0.551 | 0.530 |
| two | contaminated | 15 v 15 | 0.2 | 0.514 | 0.502 | 0.496 | 0.504 |
| two | contaminated | 30 v 30 | 0.2 | 0.814 | 0.809 | 0.816 | 0.812 |
| two | exponential | 15 v 15 | 0 | 0.318 | 0.343 | 0.391 | 0.304 |
| two | exponential | 30 v 30 | 0 | 0.524 | 0.513 | 0.542 | 0.520 |
| two | exponential | 15 v 15 | 0.2 | 0.344 | 0.325 | 0.323 | 0.318 |
| two | exponential | 30 v 30 | 0.2 | 0.569 | 0.563 | 0.562 | 0.557 |
| two | normal | 15 v 15 | 0 | 0.260 | 0.252 | 0.311 | 0.257 |
| two | normal | 30 v 30 | 0 | 0.488 | 0.476 | 0.505 | 0.489 |
| two | normal | 15 v 15 | 0.2 | 0.213 | 0.194 | 0.247 | 0.212 |
| two | normal | 30 v 30 | 0.2 | 0.404 | 0.397 | 0.422 | 0.405 |
| two | t3 | 15 v 15 | 0 | 0.352 | 0.366 | 0.416 | 0.326 |
| two | t3 | 30 v 30 | 0 | 0.566 | 0.556 | 0.593 | 0.554 |
| two | t3 | 15 v 15 | 0.2 | 0.427 | 0.408 | 0.415 | 0.418 |
| two | t3 | 30 v 30 | 0.2 | 0.733 | 0.730 | 0.730 | 0.733 |

## TOST power at a true difference of 0, bounds +/-0.5

| block | design | scenario | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| A | paired | contaminated | 8 | 0 | 0.295 | 0.250 | 0.396 | 0.300 |
| A | paired | contaminated | 15 | 0 | 0.463 | 0.398 | 0.450 | 0.454 |
| A | paired | contaminated | 30 | 0 | 0.706 | 0.572 | 0.615 | 0.689 |
| A | paired | contaminated | 50 | 0 | 0.899 | 0.801 | 0.834 | 0.896 |
| A | paired | diff_sym | 8 | 0 | 0.018 | 0.030 | 0.035 | 0.017 |
| A | paired | diff_sym | 15 | 0 | 0.120 | 0.079 | 0.166 | 0.093 |
| A | paired | diff_sym | 30 | 0 | 0.575 | 0.579 | 0.600 | 0.569 |
| A | paired | diff_sym | 50 | 0 | 0.884 | 0.888 | 0.877 | 0.883 |
| A | paired | exponential | 8 | 0 | 0.061 | 0.048 | 0.171 | 0.058 |
| A | paired | exponential | 15 | 0 | 0.311 | 0.278 | 0.422 | 0.305 |
| A | paired | exponential | 30 | 0 | 0.743 | 0.633 | 0.705 | 0.740 |
| A | paired | exponential | 50 | 0 | 0.921 | 0.855 | 0.875 | 0.923 |
| A | paired | lognormal | 8 | 0 | 0.053 | 0.047 | 0.146 | 0.058 |
| A | paired | lognormal | 15 | 0 | 0.289 | 0.254 | 0.368 | 0.280 |
| A | paired | lognormal | 30 | 0 | 0.727 | 0.659 | 0.720 | 0.720 |
| A | paired | lognormal | 50 | 0 | 0.931 | 0.884 | 0.899 | 0.929 |
| A | paired | normal | 8 | 0 | 0.034 | 0.026 | 0.089 | 0.034 |
| A | paired | normal | 15 | 0 | 0.202 | 0.180 | 0.293 | 0.196 |
| A | paired | normal | 30 | 0 | 0.684 | 0.674 | 0.716 | 0.686 |
| A | paired | normal | 50 | 0 | 0.933 | 0.933 | 0.940 | 0.933 |
| A | paired | t3 | 8 | 0 | 0.194 | 0.146 | 0.280 | 0.196 |
| A | paired | t3 | 15 | 0 | 0.444 | 0.380 | 0.468 | 0.440 |
| A | paired | t3 | 30 | 0 | 0.760 | 0.686 | 0.719 | 0.758 |
| A | paired | t3 | 50 | 0 | 0.930 | 0.882 | 0.898 | 0.928 |
| A | two | contaminated | 8 v 8 | 0 | 0.036 | 0.034 | 0.075 | 0.038 |
| A | two | contaminated | 8 v 16 | 0 | 0.086 | 0.066 | 0.118 | 0.080 |
| A | two | contaminated | 15 v 15 | 0 | 0.136 | 0.120 | 0.155 | 0.133 |
| A | two | contaminated | 15 v 30 | 0 | 0.205 | 0.169 | 0.199 | 0.190 |
| A | two | contaminated | 30 v 30 | 0 | 0.332 | 0.278 | 0.306 | 0.328 |
| A | two | contaminated | 30 v 60 | 0 | 0.487 | 0.407 | 0.444 | 0.483 |
| A | two | contaminated | 50 v 50 | 0 | 0.588 | 0.520 | 0.547 | 0.589 |
| A | two | exponential | 8 v 8 | 0 | 0.022 | 0.019 | 0.040 | 0.020 |
| A | two | exponential | 8 v 16 | 0 | 0.025 | 0.020 | 0.039 | 0.021 |
| A | two | exponential | 15 v 15 | 0 | 0.050 | 0.044 | 0.064 | 0.049 |
| A | two | exponential | 15 v 30 | 0 | 0.103 | 0.096 | 0.130 | 0.110 |
| A | two | exponential | 30 v 30 | 0 | 0.267 | 0.245 | 0.286 | 0.267 |
| A | two | exponential | 30 v 60 | 0 | 0.457 | 0.437 | 0.471 | 0.455 |
| A | two | exponential | 50 v 50 | 0 | 0.606 | 0.580 | 0.596 | 0.611 |
| A | two | likert_mid | 8 v 8 | 0 | 0.000 | 0.000 | 0.000 | 0.000 |
| A | two | likert_mid | 8 v 16 | 0 | 0.000 | 0.000 | 0.000 | 0.000 |
| A | two | likert_mid | 15 v 15 | 0 | 0.001 | 0.001 | 0.002 | 0.002 |
| A | two | likert_mid | 15 v 30 | 0 | 0.006 | 0.008 | 0.011 | 0.007 |
| A | two | likert_mid | 30 v 30 | 0 | 0.117 | 0.092 | 0.087 | 0.091 |
| A | two | likert_mid | 30 v 60 | 0 | 0.308 | 0.312 | 0.326 | 0.308 |
| A | two | likert_mid | 50 v 50 | 0 | 0.474 | 0.463 | 0.453 | 0.462 |
| A | two | lognormal | 8 v 8 | 0 | 0.007 | 0.007 | 0.018 | 0.007 |
| A | two | lognormal | 8 v 16 | 0 | 0.012 | 0.010 | 0.026 | 0.013 |
| A | two | lognormal | 15 v 15 | 0 | 0.040 | 0.036 | 0.057 | 0.036 |
| A | two | lognormal | 15 v 30 | 0 | 0.103 | 0.098 | 0.134 | 0.102 |
| A | two | lognormal | 30 v 30 | 0 | 0.291 | 0.268 | 0.312 | 0.292 |
| A | two | lognormal | 30 v 60 | 0 | 0.475 | 0.447 | 0.493 | 0.466 |
| A | two | lognormal | 50 v 50 | 0 | 0.602 | 0.576 | 0.596 | 0.605 |
| A | two | normal | 8 v 8 | 0 | 0.002 | 0.000 | 0.002 | 0.000 |
| A | two | normal | 8 v 16 | 0 | 0.002 | 0.001 | 0.006 | 0.001 |
| A | two | normal | 15 v 15 | 0 | 0.008 | 0.007 | 0.013 | 0.007 |
| A | two | normal | 15 v 30 | 0 | 0.032 | 0.032 | 0.063 | 0.032 |
| A | two | normal | 30 v 30 | 0 | 0.208 | 0.207 | 0.252 | 0.212 |
| A | two | normal | 30 v 60 | 0 | 0.422 | 0.422 | 0.466 | 0.420 |
| A | two | normal | 50 v 50 | 0 | 0.597 | 0.590 | 0.611 | 0.595 |
| A | two | t3 | 8 v 8 | 0 | 0.017 | 0.017 | 0.038 | 0.016 |
| A | two | t3 | 8 v 16 | 0 | 0.037 | 0.029 | 0.066 | 0.034 |
| A | two | t3 | 15 v 15 | 0 | 0.094 | 0.083 | 0.121 | 0.086 |
| A | two | t3 | 15 v 30 | 0 | 0.175 | 0.156 | 0.203 | 0.176 |
| A | two | t3 | 30 v 30 | 0 | 0.367 | 0.337 | 0.378 | 0.366 |
| A | two | t3 | 30 v 60 | 0 | 0.535 | 0.496 | 0.528 | 0.537 |
| A | two | t3 | 50 v 50 | 0 | 0.675 | 0.637 | 0.655 | 0.675 |
| B | paired | contaminated | 8 | 0.2 | 0.330 | 0.272 | 0.353 | 0.340 |
| B | paired | contaminated | 15 | 0.2 | 0.671 | 0.625 | 0.733 | 0.703 |
| B | paired | contaminated | 30 | 0.2 | 0.976 | 0.974 | 0.979 | 0.979 |
| B | paired | contaminated | 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | paired | normal | 8 | 0.2 | 0.054 | 0.038 | 0.066 | 0.045 |
| B | paired | normal | 15 | 0.2 | 0.162 | 0.136 | 0.200 | 0.161 |
| B | paired | normal | 30 | 0.2 | 0.570 | 0.562 | 0.631 | 0.579 |
| B | paired | normal | 50 | 0.2 | 0.882 | 0.887 | 0.897 | 0.884 |
| B | paired | t3 | 8 | 0.2 | 0.230 | 0.197 | 0.248 | 0.227 |
| B | paired | t3 | 15 | 0.2 | 0.576 | 0.523 | 0.612 | 0.581 |
| B | paired | t3 | 30 | 0.2 | 0.931 | 0.923 | 0.940 | 0.935 |
| B | paired | t3 | 50 | 0.2 | 1.000 | 0.999 | 0.999 | 1.000 |
| B | two | contaminated | 8 v 8 | 0.2 | 0.065 | 0.061 | 0.047 | 0.057 |
| B | two | contaminated | 8 v 16 | 0.2 | 0.169 | 0.163 | 0.138 | 0.159 |
| B | two | contaminated | 15 v 15 | 0.2 | 0.318 | 0.314 | 0.304 | 0.318 |
| B | two | contaminated | 15 v 30 | 0.2 | 0.528 | 0.517 | 0.511 | 0.530 |
| B | two | contaminated | 30 v 30 | 0.2 | 0.791 | 0.784 | 0.786 | 0.795 |
| B | two | contaminated | 30 v 60 | 0.2 | 0.909 | 0.909 | 0.916 | 0.912 |
| B | two | contaminated | 50 v 50 | 0.2 | 0.972 | 0.971 | 0.971 | 0.971 |
| B | two | exponential | 8 v 8 | 0.2 | 0.032 | 0.026 | 0.018 | 0.025 |
| B | two | exponential | 8 v 16 | 0.2 | 0.051 | 0.047 | 0.028 | 0.045 |
| B | two | exponential | 15 v 15 | 0.2 | 0.096 | 0.085 | 0.070 | 0.085 |
| B | two | exponential | 15 v 30 | 0.2 | 0.185 | 0.172 | 0.162 | 0.172 |
| B | two | exponential | 30 v 30 | 0.2 | 0.401 | 0.383 | 0.389 | 0.397 |
| B | two | exponential | 30 v 60 | 0.2 | 0.582 | 0.568 | 0.582 | 0.583 |
| B | two | exponential | 50 v 50 | 0.2 | 0.732 | 0.710 | 0.720 | 0.731 |
| B | two | lognormal | 8 v 8 | 0.2 | 0.018 | 0.013 | 0.008 | 0.014 |
| B | two | lognormal | 8 v 16 | 0.2 | 0.032 | 0.027 | 0.025 | 0.027 |
| B | two | lognormal | 15 v 15 | 0.2 | 0.058 | 0.049 | 0.048 | 0.052 |
| B | two | lognormal | 15 v 30 | 0.2 | 0.146 | 0.139 | 0.138 | 0.140 |
| B | two | lognormal | 30 v 30 | 0.2 | 0.350 | 0.348 | 0.347 | 0.345 |
| B | two | lognormal | 30 v 60 | 0.2 | 0.539 | 0.534 | 0.550 | 0.542 |
| B | two | lognormal | 50 v 50 | 0.2 | 0.706 | 0.698 | 0.706 | 0.702 |
| B | two | normal | 8 v 8 | 0.2 | 0.004 | 0.003 | 0.001 | 0.004 |
| B | two | normal | 8 v 16 | 0.2 | 0.004 | 0.003 | 0.002 | 0.002 |
| B | two | normal | 15 v 15 | 0.2 | 0.015 | 0.016 | 0.010 | 0.014 |
| B | two | normal | 15 v 30 | 0.2 | 0.043 | 0.042 | 0.039 | 0.041 |
| B | two | normal | 30 v 30 | 0.2 | 0.141 | 0.138 | 0.149 | 0.134 |
| B | two | normal | 30 v 60 | 0.2 | 0.299 | 0.309 | 0.332 | 0.309 |
| B | two | normal | 50 v 50 | 0.2 | 0.487 | 0.491 | 0.498 | 0.495 |
| B | two | t3 | 8 v 8 | 0.2 | 0.040 | 0.034 | 0.028 | 0.032 |
| B | two | t3 | 8 v 16 | 0.2 | 0.098 | 0.104 | 0.077 | 0.093 |
| B | two | t3 | 15 v 15 | 0.2 | 0.213 | 0.208 | 0.169 | 0.202 |
| B | two | t3 | 15 v 30 | 0.2 | 0.383 | 0.365 | 0.358 | 0.373 |
| B | two | t3 | 30 v 30 | 0.2 | 0.677 | 0.675 | 0.675 | 0.676 |
| B | two | t3 | 30 v 60 | 0.2 | 0.838 | 0.830 | 0.839 | 0.835 |
| B | two | t3 | 50 v 50 | 0.2 | 0.912 | 0.912 | 0.912 | 0.912 |

## TOST power at a true difference of 0, bounds +/-1

| block | design | scenario | n | tr | perm | boot stud | boot BCa | Welch/Yuen |
|---|---|---|---|---|---|---|---|---|
| A | paired | contaminated | 8 | 0 | 0.759 | 0.650 | 0.733 | 0.744 |
| A | paired | contaminated | 15 | 0 | 0.894 | 0.716 | 0.794 | 0.876 |
| A | paired | contaminated | 30 | 0 | 0.990 | 0.923 | 0.961 | 0.990 |
| A | paired | contaminated | 50 | 0 | 1.000 | 0.995 | 1.000 | 1.000 |
| A | paired | diff_sym | 8 | 0 | 0.468 | 0.372 | 0.569 | 0.498 |
| A | paired | diff_sym | 15 | 0 | 0.915 | 0.924 | 0.937 | 0.921 |
| A | paired | diff_sym | 30 | 0 | 0.998 | 0.999 | 0.998 | 0.999 |
| A | paired | diff_sym | 50 | 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| A | paired | exponential | 8 | 0 | 0.764 | 0.540 | 0.749 | 0.750 |
| A | paired | exponential | 15 | 0 | 0.934 | 0.793 | 0.891 | 0.933 |
| A | paired | exponential | 30 | 0 | 0.995 | 0.960 | 0.983 | 0.994 |
| A | paired | exponential | 50 | 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| A | paired | lognormal | 8 | 0 | 0.745 | 0.581 | 0.771 | 0.734 |
| A | paired | lognormal | 15 | 0 | 0.926 | 0.814 | 0.891 | 0.922 |
| A | paired | lognormal | 30 | 0 | 0.994 | 0.976 | 0.987 | 0.993 |
| A | paired | lognormal | 50 | 0 | 1.000 | 0.998 | 0.998 | 0.999 |
| A | paired | normal | 8 | 0 | 0.593 | 0.504 | 0.751 | 0.611 |
| A | paired | normal | 15 | 0 | 0.951 | 0.942 | 0.969 | 0.954 |
| A | paired | normal | 30 | 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| A | paired | normal | 50 | 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| A | paired | t3 | 8 | 0 | 0.769 | 0.623 | 0.789 | 0.765 |
| A | paired | t3 | 15 | 0 | 0.931 | 0.858 | 0.897 | 0.924 |
| A | paired | t3 | 30 | 0 | 0.986 | 0.946 | 0.961 | 0.983 |
| A | paired | t3 | 50 | 0 | 0.997 | 0.979 | 0.986 | 0.997 |
| A | two | contaminated | 8 v 8 | 0 | 0.500 | 0.412 | 0.496 | 0.483 |
| A | two | contaminated | 8 v 16 | 0 | 0.637 | 0.516 | 0.589 | 0.596 |
| A | two | contaminated | 15 v 15 | 0 | 0.678 | 0.551 | 0.609 | 0.660 |
| A | two | contaminated | 15 v 30 | 0 | 0.798 | 0.682 | 0.740 | 0.797 |
| A | two | contaminated | 30 v 30 | 0 | 0.934 | 0.863 | 0.894 | 0.934 |
| A | two | contaminated | 30 v 60 | 0 | 0.974 | 0.931 | 0.953 | 0.977 |
| A | two | contaminated | 50 v 50 | 0 | 0.993 | 0.984 | 0.987 | 0.992 |
| A | two | exponential | 8 v 8 | 0 | 0.373 | 0.308 | 0.433 | 0.362 |
| A | two | exponential | 8 v 16 | 0 | 0.495 | 0.436 | 0.544 | 0.500 |
| A | two | exponential | 15 v 15 | 0 | 0.686 | 0.606 | 0.669 | 0.677 |
| A | two | exponential | 15 v 30 | 0 | 0.816 | 0.761 | 0.816 | 0.836 |
| A | two | exponential | 30 v 30 | 0 | 0.958 | 0.931 | 0.943 | 0.957 |
| A | two | exponential | 30 v 60 | 0 | 0.980 | 0.968 | 0.979 | 0.987 |
| A | two | exponential | 50 v 50 | 0 | 0.999 | 0.996 | 0.997 | 0.999 |
| A | two | likert_mid | 8 v 8 | 0 | 0.089 | 0.095 | 0.148 | 0.128 |
| A | two | likert_mid | 8 v 16 | 0 | 0.281 | 0.258 | 0.364 | 0.282 |
| A | two | likert_mid | 15 v 15 | 0 | 0.605 | 0.570 | 0.596 | 0.567 |
| A | two | likert_mid | 15 v 30 | 0 | 0.757 | 0.764 | 0.795 | 0.759 |
| A | two | likert_mid | 30 v 30 | 0 | 0.933 | 0.938 | 0.936 | 0.935 |
| A | two | likert_mid | 30 v 60 | 0 | 0.984 | 0.984 | 0.986 | 0.983 |
| A | two | likert_mid | 50 v 50 | 0 | 0.998 | 0.997 | 0.998 | 0.997 |
| A | two | lognormal | 8 v 8 | 0 | 0.330 | 0.285 | 0.409 | 0.323 |
| A | two | lognormal | 8 v 16 | 0 | 0.514 | 0.454 | 0.570 | 0.506 |
| A | two | lognormal | 15 v 15 | 0 | 0.713 | 0.649 | 0.712 | 0.709 |
| A | two | lognormal | 15 v 30 | 0 | 0.845 | 0.809 | 0.847 | 0.847 |
| A | two | lognormal | 30 v 30 | 0 | 0.967 | 0.946 | 0.953 | 0.965 |
| A | two | lognormal | 30 v 60 | 0 | 0.983 | 0.973 | 0.980 | 0.988 |
| A | two | lognormal | 50 v 50 | 0 | 0.998 | 0.994 | 0.995 | 0.997 |
| A | two | normal | 8 v 8 | 0 | 0.231 | 0.198 | 0.351 | 0.227 |
| A | two | normal | 8 v 16 | 0 | 0.401 | 0.384 | 0.523 | 0.398 |
| A | two | normal | 15 v 15 | 0 | 0.692 | 0.679 | 0.735 | 0.689 |
| A | two | normal | 15 v 30 | 0 | 0.861 | 0.862 | 0.887 | 0.862 |
| A | two | normal | 30 v 30 | 0 | 0.965 | 0.963 | 0.970 | 0.965 |
| A | two | normal | 30 v 60 | 0 | 0.996 | 0.996 | 0.997 | 0.996 |
| A | two | normal | 50 v 50 | 0 | 1.000 | 1.000 | 1.000 | 1.000 |
| A | two | t3 | 8 v 8 | 0 | 0.476 | 0.408 | 0.529 | 0.457 |
| A | two | t3 | 8 v 16 | 0 | 0.630 | 0.550 | 0.643 | 0.608 |
| A | two | t3 | 15 v 15 | 0 | 0.767 | 0.694 | 0.745 | 0.759 |
| A | two | t3 | 15 v 30 | 0 | 0.862 | 0.792 | 0.828 | 0.854 |
| A | two | t3 | 30 v 30 | 0 | 0.956 | 0.917 | 0.928 | 0.953 |
| A | two | t3 | 30 v 60 | 0 | 0.978 | 0.949 | 0.959 | 0.978 |
| A | two | t3 | 50 v 50 | 0 | 0.991 | 0.971 | 0.975 | 0.989 |
| B | paired | contaminated | 8 | 0.2 | 0.901 | 0.796 | 0.807 | 0.913 |
| B | paired | contaminated | 15 | 0.2 | 0.998 | 0.985 | 0.993 | 0.998 |
| B | paired | contaminated | 30 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | paired | contaminated | 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | paired | normal | 8 | 0.2 | 0.503 | 0.420 | 0.687 | 0.560 |
| B | paired | normal | 15 | 0.2 | 0.877 | 0.854 | 0.945 | 0.905 |
| B | paired | normal | 30 | 0.2 | 0.998 | 0.999 | 1.000 | 0.999 |
| B | paired | normal | 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | paired | t3 | 8 | 0.2 | 0.807 | 0.698 | 0.813 | 0.833 |
| B | paired | t3 | 15 | 0.2 | 0.990 | 0.968 | 0.988 | 0.992 |
| B | paired | t3 | 30 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | paired | t3 | 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | two | contaminated | 8 v 8 | 0.2 | 0.740 | 0.697 | 0.583 | 0.727 |
| B | two | contaminated | 8 v 16 | 0.2 | 0.876 | 0.854 | 0.772 | 0.873 |
| B | two | contaminated | 15 v 15 | 0.2 | 0.974 | 0.973 | 0.964 | 0.974 |
| B | two | contaminated | 15 v 30 | 0.2 | 0.994 | 0.991 | 0.988 | 0.992 |
| B | two | contaminated | 30 v 30 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | two | contaminated | 30 v 60 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | two | contaminated | 50 v 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | two | exponential | 8 v 8 | 0.2 | 0.423 | 0.373 | 0.375 | 0.394 |
| B | two | exponential | 8 v 16 | 0.2 | 0.573 | 0.517 | 0.534 | 0.560 |
| B | two | exponential | 15 v 15 | 0.2 | 0.773 | 0.720 | 0.754 | 0.761 |
| B | two | exponential | 15 v 30 | 0.2 | 0.878 | 0.838 | 0.885 | 0.886 |
| B | two | exponential | 30 v 30 | 0.2 | 0.981 | 0.973 | 0.978 | 0.979 |
| B | two | exponential | 30 v 60 | 0.2 | 0.993 | 0.989 | 0.993 | 0.993 |
| B | two | exponential | 50 v 50 | 0.2 | 0.999 | 0.999 | 1.000 | 1.000 |
| B | two | lognormal | 8 v 8 | 0.2 | 0.356 | 0.316 | 0.341 | 0.333 |
| B | two | lognormal | 8 v 16 | 0.2 | 0.548 | 0.506 | 0.539 | 0.541 |
| B | two | lognormal | 15 v 15 | 0.2 | 0.725 | 0.700 | 0.722 | 0.719 |
| B | two | lognormal | 15 v 30 | 0.2 | 0.878 | 0.861 | 0.891 | 0.885 |
| B | two | lognormal | 30 v 30 | 0.2 | 0.984 | 0.979 | 0.982 | 0.982 |
| B | two | lognormal | 30 v 60 | 0.2 | 0.996 | 0.994 | 0.995 | 0.996 |
| B | two | lognormal | 50 v 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | two | normal | 8 v 8 | 0.2 | 0.184 | 0.167 | 0.237 | 0.180 |
| B | two | normal | 8 v 16 | 0.2 | 0.359 | 0.341 | 0.415 | 0.364 |
| B | two | normal | 15 v 15 | 0.2 | 0.578 | 0.575 | 0.628 | 0.583 |
| B | two | normal | 15 v 30 | 0.2 | 0.755 | 0.741 | 0.794 | 0.761 |
| B | two | normal | 30 v 30 | 0.2 | 0.938 | 0.940 | 0.944 | 0.940 |
| B | two | normal | 30 v 60 | 0.2 | 0.986 | 0.986 | 0.987 | 0.985 |
| B | two | normal | 50 v 50 | 0.2 | 0.998 | 0.998 | 0.998 | 0.998 |
| B | two | t3 | 8 v 8 | 0.2 | 0.608 | 0.553 | 0.542 | 0.591 |
| B | two | t3 | 8 v 16 | 0.2 | 0.774 | 0.731 | 0.738 | 0.756 |
| B | two | t3 | 15 v 15 | 0.2 | 0.937 | 0.926 | 0.933 | 0.935 |
| B | two | t3 | 15 v 30 | 0.2 | 0.975 | 0.968 | 0.976 | 0.976 |
| B | two | t3 | 30 v 30 | 0.2 | 0.999 | 0.998 | 0.999 | 0.999 |
| B | two | t3 | 30 v 60 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |
| B | two | t3 | 50 v 50 | 0.2 | 1.000 | 1.000 | 1.000 | 1.000 |

## Method failures (proportion of data sets)

| block | cell_id | design | structure | xdist | nx | ny | tr | method | fail |
|---|---|---|---|---|---|---|---|---|---|
| A | 1123 | paired | paired differences | diff_sym | 8 | 8 | 0 | boot_stud | 0.0010 |
| A | 1123 | paired | paired differences | diff_sym | 8 | 8 | 0 | boot_bca | 0.0010 |
| A | 1123 | paired | paired differences | diff_sym | 8 | 8 | 0 | ref | 0.0010 |
| A | 1127 | paired | paired differences | diff_skew | 8 | 8 | 0 | perm | 0.0005 |
| A | 1127 | paired | paired differences | diff_skew | 8 | 8 | 0 | boot_stud | 0.0010 |
| A | 1127 | paired | paired differences | diff_skew | 8 | 8 | 0 | boot_bca | 0.0010 |
| A | 1127 | paired | paired differences | diff_skew | 8 | 8 | 0 | ref | 0.0010 |

