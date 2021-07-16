bootf2
================
2021-07-17

-   [Installation](#installation)
-   [Introduction](#introduction)
-   [Examples](#examples)
    -   [Function `sim.dp()`](#function-simdp)
    -   [Function `calcf2()`](#function-calcf2)
    -   [Function `sim.dp.byf2()`](#function-simdpbyf2)
    -   [Function `bootf2()`](#function-bootf2)
-   [References](#references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

This is the development version, so it is not available on CRAN yet. You
can install it from [GitHub](https://github.com/) with
package`devtools`.

``` r
# install.packages("devtools") 
devtools::install_github("zhengguoxu/bootf2")
```

## Introduction

<!-- badges: start -->
<!-- badges: end -->

Initially, the package `bootf2` was developed to compare the dissolution
profiles using bootstrap *f*<sub>2</sub> method, as recommended recently
by several regulatory
agencies<sup>[1](#ref-EMA-2018-09-QA.MSD.DISSO)–[4](#ref-Mandula-2019-05-WS)</sup>.
Several additional functions were added later to simulate the
dissolution profiles.

There are 4 main functions:

-   `sim.dp()` to simulate dissolution profile.
-   `calcf2()` to calculate similarity factor *f*<sub>2</sub> according
    to different regulatory rules.
-   `sim.dp.byf2()` to find a dissolution profile that, when compared to
    a given reference profile, has *f*<sub>2</sub> equal to the
    predefined target *f*<sub>2</sub>.
-   `bootf2()` for the estimation of confidence intervals of
    *f*<sub>2</sub>.

The basic usage is given below as brief demonstration. The details of
functions are explained in their corresponding vignettes, while some
common topics are included in the vignette [Introduction to
bootf2](vignettes/introduction.md).

## Examples

### Function `sim.dp()`

For the most basic use, the minimum required argument is `tp`, the
vector of time points for the dissolution. To avoid changing the content
of the README every time it is run, a seed number was specified. In
practice. if it is missing, a random seed will be generated each time
when you run the function and recorded in the output for reproducibility
purpose.

``` r
library(bootf2)
# time points
tp <- c(5, 10, 15, 20, 30, 45, 60)

# simulation. simple as that. 
d.ref <- sim.dp(tp, seed = 1234)
```

The output of `sim.dp()` is a list of at least 3 components:

1.  `sim.summary`: a data frame with summary statistics of the simulated
    profiles.

``` r
print(d.ref$sim.summary)
#   product time       dp dp.cv sim.mean sim.median    sim.cv   sim.var   sim.sd
# 1 VCI3901    0  0.00000    NA  0.00000    0.00000  0.000000  0.000000 0.000000
# 2 VCI3901    5 42.55948    NA 49.53235   48.82096 16.111382 63.685913 7.980345
# 3 VCI3901   10 66.29877    NA 72.54820   73.65294  9.699700 49.518780 7.036958
# 4 VCI3901   15 79.65966    NA 83.61915   84.72814  7.743776 41.929253 6.475280
# 5 VCI3901   20 87.19983    NA 89.00044   90.38211  6.026375 28.767133 5.363500
# 6 VCI3901   30 93.87387    NA 93.30767   94.08669  3.873792 13.064937 3.614545
# 7 VCI3901   45 96.45335    NA 95.07656   95.36508  2.823389  7.205900 2.684381
# 8 VCI3901   60 96.92426    NA 95.54486   95.41997  2.622047  6.276184 2.505231
#    sim.min  sim.max sim.qt05 sim.qt25 sim.qt75 sim.qt95
# 1  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
# 2 36.11408 62.97400 37.42962 44.92917 54.99613 60.46630
# 3 54.69761 81.81401 61.36394 70.53716 76.08199 80.48728
# 4 66.59627 90.85062 73.52548 80.95547 87.72950 89.92357
# 5 74.63857 94.68192 80.92610 86.40894 92.23759 94.17383
# 6 84.21522 98.62731 87.08164 93.13322 95.01842 96.96281
# 7 90.32152 99.76309 90.53187 94.07949 96.43587 98.68313
# 8 90.42323 99.85618 92.02013 94.29310 97.15213 99.01369
```

2.  `sim.disso`: a data frame of dissolution profiles that has the
    correct format to be used as the input for the function `calcf2()`.

``` r
print(d.ref$sim.disso)
#   time  unit.01  unit.02  unit.03  unit.04  unit.05  unit.06  unit.07  unit.08
# 1    0  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
# 2    5 49.29328 44.96255 54.17509 36.11408 38.50597 48.34864 47.36911 58.41454
# 3   10 71.14519 66.81821 74.64047 54.69761 75.74942 72.66540 68.71305 79.40177
# 4   15 81.29123 79.19484 84.95419 66.59627 90.85062 84.50209 79.94818 88.28366
# 5   20 86.07044 86.50320 90.56083 74.63857 94.68192 90.20340 86.12617 92.22997
# 6   30 89.42690 93.64105 95.60094 84.21522 95.44273 94.23164 91.60972 94.87699
# 7   45 90.32152 97.06950 97.79953 90.70398 95.45178 95.27837 93.86045 95.56065
# 8   60 90.42323 97.95304 98.32438 93.32669 95.45179 95.38816 94.33537 95.63962
#    unit.09  unit.10  unit.11  unit.12
# 1  0.00000  0.00000  0.00000  0.00000
# 2 51.94259 62.97400 57.45928 44.82905
# 3 77.07971 81.81401 75.46822 72.38533
# 4 87.54478 89.16509 84.37901 86.71988
# 5 91.67441 92.26044 89.29779 93.75812
# 6 93.83528 94.24255 93.94175 98.62731
# 7 94.15250 94.73265 96.22466 99.76309
# 8 94.16632 94.78840 96.88516 99.85618
```

3.  `sim.info`: a data frame with information of the simulation.

``` r
print(d.ref$sim.info)
#    fmax fmax.cv  mdt mdt.cv     beta beta.cv tlag tlag.cv seed n.units
# 1 97.03       3 8.69     30 0.993764      20    0       0 1234      12
#   max.disso   model time.unit
# 1       105 Weibull       min
```

Depending on the argument settings, there might be additional 2
components:

4.  `model.par.ind`: a data frame of individual model parameters that
    are used to simulate the individual dissolution profiles if
    mathematical models are used for the simulation.

``` r
print(d.ref$model.par.ind)
#    fmax.ind   mdt.ind  beta.ind tlag.ind
# 1  90.43662  6.392185 0.9720844        0
# 2  98.28721  8.650624 0.8972171        0
# 3  98.51432  6.562610 0.8282050        0
# 4  95.37133 12.095847 0.8405575        0
# 5  95.45179  7.534519 1.6110909        0
# 6  95.40077  7.024057 1.0207750        0
# 7  94.47347  7.476730 0.9008715        0
# 8  95.65084  5.330499 0.9099505        0
# 9  94.16689  6.121990 1.0894384        0
# 10 94.79651  4.518414 0.8650227        0
# 11 97.21781  5.811711 0.7438644        0
# 12 99.86358  7.955669 1.1148230        0
```

5.  `sim.plot`: a plot if the argument `plot` is set to be `TRUE`.

``` r
print(d.ref$sim.plot)
```

<img src="man/figures/README-simdp-plot-1.png" width="100%" />

The complete arguments are shown below. More details can be viewed in
the package vignette by `vignette("sim.dp", package = "bootf2")`.

``` r
dat <- sim.dp(tp, model = c("Weibull", "first-order"), model.par,
              seed, product, dp, dp.cv, ascending = FALSE, n.units = 12L,
              max.disso = 105, message = FALSE, plot = TRUE,
              time.unit = c("min", "h"), plot.max.unit = 36L)
```

### Function `calcf2()`

The minimum required arguments are dissolution profiles for test and
reference, e.g., `calcf2(test, ref)`. Data can be read from an Excel
file. For interactive use, the `test` and `ref` can be data frames with
the time as the first column and individual units as the rest columns.
The `sim.disso` data frame in the output of `sim.dp()` comes with the
correct format, as shown above.

``` r
# simulate a test data
d.test <- sim.dp(tp, seed = 100)

# calculate f2 with default settings
calcf2(d.test$sim.disso, d.ref$sim.disso)
# The f2 method was applied according to EMA's BE guideline.
# 
# Individual data was provided with option 'cv.rule = TRUE',
# therefore, CV has been calculated and checked accordingly.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         46.82      49.53      17.44      16.11 
# 10        66.28      72.55      10.85       9.70 
# 15        77.54      83.62       7.41       7.74 
# 20        84.51      89.00       5.59       6.03 
# ------------------------------------------------ 
# 30        91.92      93.31       4.02       3.87 
# 45        96.19      95.08       2.90       2.82 
# 60        97.71      95.54       2.24       2.62 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV criteria fulfilled; therefore, f2 method can be applied.
# 
# The time points above the dashed line are used in f2 calculation.
# 
# Estimated f2 = 64.24
```

<img src="man/figures/README-calcf2-good-1.png" width="100%" />

When the conditions to apply *f*<sub>2</sub> are not fulfilled, the
function will stop and show some messages.

``` r
# simulate reference with CV% not fulfil the criterion 
d.r2 <- sim.dp(tp, seed = 456)
d.r2$sim.summary
```

<div data-pagedtable="false">

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["product"],"name":[1],"type":["chr"],"align":["left"]},{"label":["time"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["dp"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["dp.cv"],"name":[4],"type":["lgl"],"align":["right"]},{"label":["sim.mean"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["sim.median"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["sim.cv"],"name":[7],"type":["dbl"],"align":["right"]},{"label":["sim.var"],"name":[8],"type":["dbl"],"align":["right"]},{"label":["sim.sd"],"name":[9],"type":["dbl"],"align":["right"]},{"label":["sim.min"],"name":[10],"type":["dbl"],"align":["right"]},{"label":["sim.max"],"name":[11],"type":["dbl"],"align":["right"]},{"label":["sim.qt05"],"name":[12],"type":["dbl"],"align":["right"]},{"label":["sim.qt25"],"name":[13],"type":["dbl"],"align":["right"]},{"label":["sim.qt75"],"name":[14],"type":["dbl"],"align":["right"]},{"label":["sim.qt95"],"name":[15],"type":["dbl"],"align":["right"]}],"data":[{"1":"LWL4216","2":"0","3":"0.00000","4":"NA","5":"0.00000","6":"0.00000","7":"0.000000","8":"0.00000","9":"0.000000","10":"0.00000","11":"0.00000","12":"0.00000","13":"0.00000","14":"0.00000","15":"0.00000"},{"1":"LWL4216","2":"5","3":"40.34113","4":"NA","5":"40.58672","6":"41.80670","7":"26.432697","8":"115.09354","9":"10.728166","10":"19.53981","11":"61.22352","12":"23.47979","13":"36.71898","14":"45.57355","15":"54.47309"},{"1":"LWL4216","2":"10","3":"62.35050","4":"NA","5":"64.46633","6":"66.34211","7":"17.836940","8":"132.22290","9":"11.498822","10":"41.70382","11":"81.85290","12":"43.97261","13":"60.47527","14":"71.65643","15":"77.09973"},{"1":"LWL4216","2":"15","3":"75.43591","4":"NA","5":"77.95486","6":"82.49419","7":"14.125419","8":"121.25205","9":"11.011451","10":"59.68043","11":"90.92986","12":"59.78431","13":"70.18733","14":"85.40015","15":"90.46134"},{"1":"LWL4216","2":"20","3":"83.40617","4":"NA","5":"85.36077","6":"89.56867","7":"11.079036","8":"89.43769","9":"9.457150","10":"70.24499","11":"96.35407","12":"71.57080","13":"76.68573","14":"93.20216","15":"95.73646"},{"1":"LWL4216","2":"30","3":"91.39903","4":"NA","5":"91.98773","6":"93.49207","7":"6.849871","8":"39.70312","9":"6.301041","10":"83.11474","11":"99.76324","12":"83.41645","13":"86.46009","14":"97.33652","15":"99.13142"},{"1":"LWL4216","2":"45","3":"95.31945","4":"NA","5":"95.20269","6":"94.56145","7":"4.193580","8":"15.93926","9":"3.992401","10":"88.70081","11":"101.45909","12":"90.06261","13":"92.38262","14":"98.78765","15":"100.36438"},{"1":"LWL4216","2":"60","3":"96.33039","4":"NA","5":"96.26318","6":"96.59041","7":"3.443959","8":"10.99098","9":"3.315264","10":"91.25553","11":"101.60724","12":"91.75410","13":"93.91981","14":"98.91278","15":"100.52078"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

</div>

``` r
calcf2(d.test$sim.disso, d.r2$sim.disso)
# The f2 method was applied according to EMA's BE guideline.
# 
# Individual data was provided with option 'cv.rule = TRUE',
# therefore, CV has been calculated and checked accordingly.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         46.82      40.59      17.44      26.43 
# 10        66.28      64.47      10.85      17.84 
# 15        77.54      77.95       7.41      14.13 
# 20        84.51      85.36       5.59      11.08 
# ------------------------------------------------ 
# 30        91.92      91.99       4.02       6.85 
# 45        96.19      95.20       2.90       4.19 
# 60        97.71      96.26       2.24       3.44 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV criteria not fulfilled; therefore, f2 method cannot be applied.
# Error in calcf2(d.test$sim.disso, d.r2$sim.disso): You should consider alternative methods such as bootstrap f2.
```

The complete arguments are show below. More details can be viewed in the
package vignette by `vignette("calcf2", package = "bootf2")`.

``` r
dat <- calcf2(test, ref, regulation = c("EMA", "FDA", "WHO"),
              path.in, file.in, path.out, file.out, digits = 2L,
              cv.rule = TRUE, min.points = 3L, both.TR.85 = FALSE,
              f2.type = c("est.f2", "exp.f2", "bc.f2", "vc.exp.f2",
                          "vc.bc.f2", "all"), plot = TRUE,
              message = TRUE, time.unit = c("min", "h"),
              plot.start.time = 0, plot.max.unit = 24L)
```

### Function `sim.dp.byf2()`

Given any dissolution profile `dp` at time points `tp`, and target
*f*<sub>2</sub> value, this function will find another dissolution
profile such that when the newly simulated profile is compared to the
`dp`, the calculated *f*<sub>2</sub> will be equal to the target
*f*<sub>2</sub>.

``` r
# mean dissolution profile fo tp
dp <- c(51, 66, 75, 81, 88, 92, 95)

# find another profile with target f2 = 60
d.t2 <- sim.dp.byf2(tp, dp, target.f2 = 60, seed = 123)
```

<img src="man/figures/README-simdpbyf2-dat-1.png" width="100%" />

    # Obtained model parameters and calculated f2 are:
    #     model seed     fmax tlag      mdt      beta f2 time.point regulation
    # 1 Weibull  123 96.95502    0 4.424869 0.3611346 60          5        EMA
    # 
    # And the difference between simulated test and reference is:
    #   time ref     test   diff.tr
    # 1    0   0  0.00000  0.000000
    # 2    5  51 62.86075 11.860752
    # 3   10  66 71.62829  5.628286
    # 4   15  75 76.46007  1.460075
    # 5   20  81 79.66658 -1.333418
    # 6   30  88 83.78233 -4.217670
    # 7   45  92 87.33950 -4.660501
    # 8   60  95 89.48884 -5.511159

If you don’t like the shape of the profile, you can run the function
many times with different `seed` number until you get the shape that you
are willing to accept.

The model parameters in the output are more useful in simulation studies
since they can be used as initial model parameter input to the function
`sim.dp()` to simulate a large population of dissolution profiles that
have known *f*<sub>2</sub> value when compared to target dissolution
profile.

The complete arguments are show below. More details can be viewed in the
package vignette by `vignette("sim.dp.byf2", package = "bootf2")`.

``` r
dat <- sim.dp.byf2(tp, dp, sim.dp.out, target.f2, seed,
                   regulation = c("EMA", "FDA", "WHO"),
                   model = c("Weibull", "first-order"), digits = 2L,
                   min.points = 3L, both.TR.85 = FALSE, max.disso = 105,
                   model.par.cv = 50, fix.fmax.cv = 0, random.factor = 3,
                   sim.target = c("ref.pop", "ref.median", "ref.mean"),
                   time.unit = c("min", "h"), message = TRUE, plot = TRUE)
```

### Function `bootf2()`

The minimum required arguments are dissolution profiles of `test` and
`ref`. The function can output many different 90% confidence intervals
for several *f*<sub>2</sub> estimators.

``` r
# get test and reference data set with correct format
test <- d.test$sim.disso
ref  <- d.ref$sim.disso

# use most default settings (output all) but small number of bootstrap
# to have shorter run time for the example. default n.boots = 10000L
t_vs_r <- bootf2(test, ref, n.boots = 100L)
# =================================================================
# |                                                               |
# |  Comparison of Dissolution Profiles by Bootstrap f2 Method.   |
# |_______________________________________________________________|
# |                                                               |
# | Smimilarity Criterion:                                        |
# | ----------------------                                        |
# |     To conclude similarity, the lower limit of 90% confidence |
# | interval should be greater than or equal to 50.               |
# |                                                               |
# =================================================================
# 
# ============================================
# |              Main Results:               |
# |  f2 and Its Confidence Intervals (CI)s   |
# ============================================
# 
# ----------------------
# * Estimated f2 Values
# ----------------------
#   - with original data                :   64.24
#   - with original data (by Jackknife) :   64.26
#   - with bootstrapped data (Mean)     :   65.58
#   - with bootstrapped data (Median)   :   64.24
# 
# -----------------------
# * Confidence Intervals
# -----------------------
#          Types of         Lower   Upper
#    Confidence Intervals   <----------->
#                  Normal   46.94   78.86
#                   Basic   43.90   75.45
#     Percentile (Type 1)   53.01   82.77
#     Percentile (Type 2)   53.23   83.72
#     Percentile (Type 3)   53.01   82.77
#     Percentile (Type 4)   53.01   82.77
#     Percentile (Type 5)   53.23   83.72
#     Percentile (Type 6)   53.03   84.58
#     Percentile (Type 7)   53.42   82.86
#     Percentile (Type 8)   53.16   84.01
#     Percentile (Type 9)   53.18   83.94
#     Percentile (boot)     53.03   84.57
#         BCa (Jackknife)   53.20   86.19
#              BCa (boot)   53.39   88.46
#   ----------------------------------------------------
#   Out of 100 bootstrapped data sets,
#   - Number of f2 calculated with 1 time point :   0
#   - Number of f2 calculated with 2 time point :   0
#   - Number of f2 cannot be calculated (NA)    :   0
#   ----------------------------------------------------
# ______________________________________________________________________
# 
# ---------------------
# * Expected f2 Values
# ---------------------
#   - with original data                :   61.59
#   - with original data (by Jackknife) :   61.46
#   - with bootstrapped data (Mean)     :   62.02
#   - with bootstrapped data (Median)   :   62.08
# 
# -----------------------
# * Confidence Intervals
# -----------------------
#          Types of         Lower   Upper
#    Confidence Intervals   <----------->
#                  Normal   50.08   72.24
#                   Basic   50.57   71.73
#     Percentile (Type 1)   51.39   71.92
#     Percentile (Type 2)   51.87   72.28
#     Percentile (Type 3)   51.39   71.92
#     Percentile (Type 4)   51.39   71.92
#     Percentile (Type 5)   51.87   72.28
#     Percentile (Type 6)   51.44   72.61
#     Percentile (Type 7)   52.30   71.95
#     Percentile (Type 8)   51.73   72.39
#     Percentile (Type 9)   51.77   72.36
#     Percentile (boot)     51.45   72.61
#         BCa (Jackknife)   48.22   71.27
#              BCa (boot)   47.95   71.16
#   ----------------------------------------------------
#   Out of 100 bootstrapped data sets,
#   - Number of f2 calculated with 1 time point :   0
#   - Number of f2 calculated with 2 time point :   0
#   - Number of f2 cannot be calculated (NA)    :   0
#   ----------------------------------------------------
# ______________________________________________________________________
# 
# ---------------------------
# * Bias-Corrected f2 Values
# ---------------------------
#   - with original data                :   67.75
#   - with original data (by Jackknife) :   68.09
#   - with bootstrapped data (Mean)     :   68.32
#   - with bootstrapped data (Median)   :   66.14
# 
# -----------------------
# * Confidence Intervals
# -----------------------
#          Types of         Lower   Upper
#    Confidence Intervals   <----------->
#                  Normal   47.44   86.92
#                   Basic   42.87   82.41
#     Percentile (Type 1)   54.34   91.57
#     Percentile (Type 2)   54.34   91.57
#     Percentile (Type 3)   54.34   90.19
#     Percentile (Type 4)   52.82   90.75
#     Percentile (Type 5)   54.34   91.44
#     Percentile (Type 6)   53.01   92.70
#     Percentile (Type 7)   54.34   90.81
#     Percentile (Type 8)   54.15   91.73
#     Percentile (Type 9)   54.29   91.61
#     Percentile (boot)     53.09   92.63
#         BCa (Jackknife)   55.66   100.33
#              BCa (boot)   56.31   103.99
#   ----------------------------------------------------
#   Out of 100 bootstrapped data sets,
#   - Number of f2 calculated with 1 time point :   0
#   - Number of f2 calculated with 2 time point :   0
#   - Number of f2 cannot be calculated (NA)    :   8
#   ----------------------------------------------------
# ______________________________________________________________________
# 
# ----------------------------------------
# * Variance-corrected Expected f2 Values
# ----------------------------------------
#   - with original data                :   61.58
#   - with original data (by Jackknife) :   61.44
#   - with bootstrapped data (Mean)     :   61.82
#   - with bootstrapped data (Median)   :   62.01
# 
# -----------------------
# * Confidence Intervals
# -----------------------
#          Types of         Lower   Upper
#    Confidence Intervals   <----------->
#                  Normal   50.44   72.24
#                   Basic   50.97   71.74
#     Percentile (Type 1)   51.37   70.97
#     Percentile (Type 2)   51.83   71.62
#     Percentile (Type 3)   51.37   70.97
#     Percentile (Type 4)   51.37   70.97
#     Percentile (Type 5)   51.83   71.62
#     Percentile (Type 6)   51.42   72.20
#     Percentile (Type 7)   52.24   71.04
#     Percentile (Type 8)   51.69   71.81
#     Percentile (Type 9)   51.73   71.76
#     Percentile (boot)     51.42   72.20
#         BCa (Jackknife)   48.30   70.90
#              BCa (boot)   47.75   70.86
#   ----------------------------------------------------
#   Out of 100 bootstrapped data sets,
#   - Number of f2 calculated with 1 time point :   0
#   - Number of f2 calculated with 2 time point :   0
#   - Number of f2 cannot be calculated (NA)    :   0
#   ----------------------------------------------------
# ______________________________________________________________________
# 
# -----------------------------------------
# * Variance- and Bias-Corrected f2 Values
# -----------------------------------------
#   - with original data                :   67.76
#   - with original data (by Jackknife) :   68.13
#   - with bootstrapped data (Mean)     :   69.11
#   - with bootstrapped data (Median)   :   66.49
# 
# -----------------------
# * Confidence Intervals
# -----------------------
#          Types of         Lower   Upper
#    Confidence Intervals   <----------->
#                  Normal   43.85   88.99
#                   Basic   40.14   82.31
#     Percentile (Type 1)   54.47   94.89
#     Percentile (Type 2)   54.47   94.89
#     Percentile (Type 3)   54.47   92.17
#     Percentile (Type 4)   52.95   93.26
#     Percentile (Type 5)   54.48   94.62
#     Percentile (Type 6)   53.14   95.42
#     Percentile (Type 7)   54.49   93.40
#     Percentile (Type 8)   54.28   94.97
#     Percentile (Type 9)   54.43   94.91
#     Percentile (boot)     53.22   95.39
#         BCa (Jackknife)   54.84   101.37
#              BCa (boot)   55.30   101.78
#   ----------------------------------------------------
#   Out of 100 bootstrapped data sets,
#   - Number of f2 calculated with 1 time point :   0
#   - Number of f2 calculated with 2 time point :   0
#   - Number of f2 cannot be calculated (NA)    :   8
#   ----------------------------------------------------
# ______________________________________________________________________
# 
# ============================================
# | Function Settings and System Information |
# ============================================
# 
# ---------------------
# * Function Settings
# ---------------------
#   - test              :   test
#   - ref               :   ref
#   - n.boots           :   100
#   - seed              :   306
#   - digits            :   2
#   - alpha             :   0.05 (90% CI)
#   - regulation        :   EMA
#   - min.points        :   1
#   - both.TR.85        :   FALSE
#   - print.report      :   TRUE
#   - report.style      :   concise
#   - f2.type           :   all
#   - ci.type           :   all
#   - quantile.type     :   all
#   - jackknife.type    :   nt+nr
#   - time.unit         :   min
#   - output.to.screen  :   TRUE
#   - sim.data.out      :   FALSE
#   - path.in           :   NA
#   - file.in           :   NA
#   - path.out          :   /home/zhengguo/github/bootf2/
#   - file.out          :   test_vs_ref_CEST_2021-07-17_011256.txt
# 
# ---------------------
# * System Information
# ---------------------
#   - Operating System Name     :   Linux 5.4.0-77-generic
#   - Operating System Version  :   #86-Ubuntu SMP Thu Jun 17 02:35:03 UTC 2021
#   - Machine Node Name         :   MyHomeLinuxPC
#   - User Name                 :   zhengguo
#   - Time Zone                 :   Europe/Madrid
#   - R Version                 :   4.1.0 (2021-05-18)
# ______________________________________________________________________
# 
# The current report was generated at 01:12:56 on 2021-07-17 CEST by
# user 'zhengguo' on machine 'MyHomeLinuxPC', and saved as text file
# 'test_vs_ref_CEST_2021-07-17_011256.txt' at the location:
# '/home/zhengguo/github/bootf2/'.
# ======================================================================
```

With default settings, the function prints all confidence intervals for
all *f*<sub>2</sub> estimators, and the result will be save in a text
file.

The complete arguments are shown below. More details can be viewed in
the package vignette by `vignette("bootf2", package = "bootf2")`.

``` r
dat <- bootf2(test, ref, path.in, file.in, path.out, file.out,
                   n.boots = 10000L, seed = 306L, digits = 2L, alpha = 0.05,
                   regulation = c("EMA", "FDA", "WHO"), min.points = 1L,
                   both.TR.85 = FALSE, print.report = TRUE,
                   report.style = c("concise", "intermediate", "detailed"),
                   f2.type = c("all", "est.f2", "exp.f2", "bc.f2",
                               "vc.exp.f2", "vc.bc.f2"),
                   ci.type = c("all", "normal", "basic", "percentile",
                               "bca.jackknife", "bca.boot"),
                   quantile.type = c("all", 1:9, "boot"),
                   jackknife.type = c("nt+nr", "nt*nr", "nt=nr"),
                   time.unit = c("min", "h"), output.to.screen = TRUE,
                   sim.data.out = FALSE)
```

## References

<div id="refs" class="references csl-bib-body">

<div id="ref-EMA-2018-09-QA.MSD.DISSO" class="csl-entry">

<span class="csl-left-margin">(1) </span><span
class="csl-right-inline">European Medicines Agency. Question and answer
on the adequacy of the Mahalanobis distance to assess the comparability
of drug dissolution profiles
<https://www.ema.europa.eu/documents/scientific-guideline/question-answer-adequacy-mahalanobis-distance-assess-comparability-drug-dissolution-profiles_en.pdf>
(accessed 2018 -12 -04).</span>

</div>

<div id="ref-Davit-2013-03-BA" class="csl-entry">

<span class="csl-left-margin">(2) </span><span
class="csl-right-inline">Davit, B. M.; Stier, E.; Jiang, X.; Anand, O.
Expectations of the US-FDA Regarding Dissolution Data in Generic Drug
Regulatory Submissions. *Biopharma Asia* **2013**, *2* (2).</span>

</div>

<div id="ref-Lum-2019-05-WS" class="csl-entry">

<span class="csl-left-margin">(3) </span><span
class="csl-right-inline">Lum, S. Health Canada’s Current Practice and
Challenges in the Evaluation of Dissolution Profile Comparisons in
Support of Minor/Moderate Product Quality Changes: Case Studies. In *In
vitro dissolution profiles similarity assessment in support of drug
product quality: What, how, and when*; 2019.</span>

</div>

<div id="ref-Mandula-2019-05-WS" class="csl-entry">

<span class="csl-left-margin">(4) </span><span
class="csl-right-inline">Mandula, H. Rational Statistical Analysis
Practice in Dissolution Profile Comparison: FDA Perspective. In *In
vitro dissolution profiles similarity assessment in support of drug
product quality: What, how, and when*; 2019.</span>

</div>

</div>
