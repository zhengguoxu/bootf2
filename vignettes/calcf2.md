---
title: "Calculating Similarity Factor $f_2$"
author: "Zhengguo XU"
date: "2021-07-26"
output: 
  rmarkdown::html_vignette:
    keep_md: true
    toc: true
    toc_depth: 3
    fig_width: 7
    fig_height: 4.5
  highlight: "tango"
bibliography: ref.bib
notes-after-punctuation: false
link-citations: yes
csl: ref.csl
vignette: >
  %\VignetteIndexEntry{Calculating Similarity Factor $f_2$}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---






## Introduction

To use traditional $f_2$ method, several conditions defined in the
regulatory guidelines have to be fulfilled. Different guidelines phrase
the conditions differently, which led to confusion some times. The details
were explained in the vignette [Introduction to bootf2](introduction.md). 

This document gives some examples of the usage of the function `calcf2()`. 

## Usage

The complete list of arguments of the function is as follows:


```r
calcf2(test, ref, regulation = c("EMA", "FDA", "WHO"),
       path.in, file.in, path.out, file.out, digits = 2L,
       cv.rule = TRUE, min.points = 3L, both.TR.85 = FALSE,
       f2.type = c("est.f2", "exp.f2", "bc.f2", "vc.exp.f2",
                   "vc.bc.f2", "all"), plot = TRUE,
       message = TRUE, time.unit = c("min", "h"),
       plot.start.time = 0, plot.max.unit = 24L)
```

### Notes on function arguments

1. Data input: `test`, `ref`, `path.in`, `file.in`
    - In the typical interactive use, `test` and `ref` are data frames with
      the time as the first column, and individual units for the rest of 
      columns. In such cases, arguments `path.in` and `file.in` should not 
      be used.
    - Data can be read directly from an Excel file with extension `.xlsx` 
      or `.xls`. In this case, data of test and reference should be 
      stored in separate worksheets. The first column should be time, 
      the rest columns are individual units. *The first row* is the 
      column head indicating the names, such as 'time', 'unit 01', 
      unit 02', .... It doesn't matter what the names are as columns 
      will be renamed internally by the function. The important point 
      is that *the first row will not be read, so do not put dissolution 
      data on the first row*. 
    - When `path.in` and `file.in` are provided, the argument `test` and
      `ref` should be the worksheet names inside quotation mark, e.g., 
      `"lot ABCD1234 pH 6.8"`. 
    - `path.in` can be an absolute path such as`"/home/myname/my.project/dat/"`,
      or a relative path such as `"../dat/"` if the working directory is 
      in the folder `"/home/myname/my.project/analysis/"` and the data file
      is in the folder `"/home/myname/my.project/dat/"`. 
    - One more note for Windows user: As Windows use "\" instead of "/" to
      separate path, you have to either escape it by an additional "\", 
      e.g., `"C:\user\myname\my.project\dat\"` *cannot* be the `path.in`, 
      you have to changed it to `"C:\\user\\myname\\my.project\\dat\\"`, or
      to `"C:/user/myname/my.project/dat/"`, the same format as used in Linux
      system.
1. For such simple calculation, argument `path.out` and `file.out` are
   really overkill. By default the result will be printed to screen. But if 
   somehow you feel that you need it, the same principle for `path.in` and 
   `file.in` applies here.
1. The $f_2$ value should be reported without decimal.
1. When argument `cv.rule` is set to `TRUE`, and if argument `message` is also 
   `TRUE`, additional message will be printed to screen. One particular 
   note is the *rounding of CV before checking compliance*. Since the
   guidelines does not specify the precision for CV, when checking if CV
   condition defined in the guidelines is fulfilled, the calculated CV 
   values are rounded with `digits` decimal first before they were compared
   to the 20%/10% criterion. And refer to the vignette [Introduction to
   bootf2](introduction.md) for details of regulation rules. 
1. In terms of how many time points to be used in the calculation (especially
   those points with dissolution more than 85%), the `"EMA"` and `"FDA"` 
   are the same. See vignette [Introduction to bootf2](introduction.md) for 
   details. The reason both of them are included in `regulation` argument
   is because they have different rules for CV.
1. If the argument `both.TR.85 = TRUE`, *and* `regulation = FDA`, then the
   function will calculate $f_2$ using all time points until both test 
   and reference dissolve more than 85%. This is the **wrong interpretation**
   of the guidelines, as explained in the vignette [Introduction to
   bootf2](introduction.md). The purpose of this argument is for historical 
   reason, in case users want to check the calculation published in the old
   articles where the wrong interpretation was used, or for validation 
   purpose. It should be set to `FALSE` for the daily analysis.
1. For traditional $f_2$ calculation, do not specify the argument
   `f2.type`. This argument is only for bootstrap method. The default
   `est.f2` is the correct one to use. See vignette [Confidence Intervals of
   $f_2$ Using Bootstrap Method](bootf2.md) for details.
1. The argument `plot.max.unit` control how individual profile will be 
   represented in the plot. When the actual number of units is greater than 
   the value of `plot.max.unit`, the individual profile will be represented 
   as boxplots at each time points. 
1. For delayed-release formulation, it is typical that the product will be 
   put in acidic media for 2 hours where almost no dissolution will
   occur, then transferred to different medium (e.g., pH 6.8 medium). 
   In such special case, setting `plot.start.time = 120` (or 
   `plot.start.time = 2` if the time unit is `"h"`) for example, 
   can make plots more readable.
1. Read the function manual by `help("calcf2")` for more details on
   each argument.

## Examples

## Simple case with low variability

First, let's simulate profiles with low variability for reference and test.

```r
# time points
tp <- c(5, 10, 15, 20, 30, 45, 60)

# model.par for reference with low variability
par.r1.lv <- list(fmax = 100, fmax.cv = 3, mdt = 15, mdt.cv = 14, 
                  tlag = 0, tlag.cv = 0, beta = 1.5, beta.cv = 8)

# simulate reference data
dr1.lv <- sim.dp(tp, model.par = par.r1.lv, seed = 100, plot = FALSE)

# model.par for test 
par.t1.lv <- list(fmax = 100, fmax.cv = 3, mdt = 12.29, mdt.cv = 12,
                  tlag = 0, tlag.cv = 0, beta = 1.727, beta.cv = 9)

# simulate test data with low variability
dt1.lv <- sim.dp(tp, model.par = par.t1.lv, seed = 100, plot = FALSE)
```

Calculate f2 with default setting (following EMA's guideline)

```r
t1.lv_vs_r1.lv <- calcf2(dt1.lv$sim.disso, dr1.lv$sim.disso)
# The f2 method was applied according to EMA's BE guideline.
# 
# Individual data was provided with option 'cv.rule = TRUE',
# therefore, CV has been calculated and checked accordingly.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         19.21      17.69      19.92      19.65 
# 10        51.21      42.73      11.99      13.54 
# 15        76.81      64.36       6.60       9.27 
# 20        91.12      79.74       3.16       6.03 
# ------------------------------------------------ 
# 30        99.27      94.78       1.53       2.28 
# 45       100.01      99.57       1.53       1.50 
# 60       100.01      99.99       1.53       1.53 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV criteria fulfilled; therefore, f2 method can be applied.
# 
# The time points above the dashed line are used in f2 calculation.
# 
# Estimated f2 = 51.06
```

![](/home/zhengguo/github/bootf2/vignettes/calcf2_files/figure-html/calcf2-cvok01-a-1.png)<!-- -->

calculate f2 following WHO guideline

```r
t.lv_vs_r.lv <- calcf2(dt1.lv$sim.disso, dr1.lv$sim.disso, regulation = "WHO")
# The f2 method was applied according to WHO's BE guideline.
# 
# Individual data was provided with option 'cv.rule = TRUE',
# therefore, CV has been calculated and checked accordingly.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         19.21      17.69      19.92      19.65 
# 10        51.21      42.73      11.99      13.54 
# 15        76.81      64.36       6.60       9.27 
# 20        91.12      79.74       3.16       6.03 
# 30        99.27      94.78       1.53       2.28 
# ------------------------------------------------ 
# 45       100.01      99.57       1.53       1.50 
# 60       100.01      99.99       1.53       1.53 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV criteria fulfilled; therefore, f2 method can be applied.
# 
# The time points above the dashed line are used in f2 calculation.
# 
# Estimated f2 = 52.86
```

![](/home/zhengguo/github/bootf2/vignettes/calcf2_files/figure-html/calcf2-cvok01-b-1.png)<!-- -->


Calculate f2 following FDA guidance, to confirm "old calculation" with
*wrong interpretation*. See [Introduction to bootf2](vignettes/introduction.md)
for details.

```r
# 
t.lv_vs_r.lv <- calcf2(dt1.lv$sim.disso, dr1.lv$sim.disso, regulation = "FDA",
                       both.TR.85 = TRUE)
# The f2 method was applied according to FDA's BE guidance.
# 
# Individual data was provided with option 'cv.rule = TRUE',
# therefore, CV has been calculated and checked accordingly.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         19.21      17.69      19.92      19.65 
# 10        51.21      42.73      11.99      13.54 
# 15        76.81      64.36       6.60       9.27 
# 20        91.12      79.74       3.16       6.03 
# 30        99.27      94.78       1.53       2.28 
# ------------------------------------------------ 
# 45       100.01      99.57       1.53       1.50 
# 60       100.01      99.99       1.53       1.53 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV criteria fulfilled; therefore, f2 method can be applied.
# 
# The time points above the dashed line are used in f2 calculation.
# 
# Estimated f2 = 52.86*
# 
# *Note: Argument 'both.TR.85' is 'TRUE', which is the wrong interpretation
# of the guidance. This should only be used for cases such as checking the
# calculation published in the old literature.
```

![](/home/zhengguo/github/bootf2/vignettes/calcf2_files/figure-html/calcf2-cvok01-c-1.png)<!-- -->

# Cases with high variability
Simulate profiles with the same population parameters but with higher
variability.


```r
# model.par for reference with high variability
par.r1.hv <- list(fmax = 100, fmax.cv = 3, mdt = 15, mdt.cv = 20, 
                  tlag = 0, tlag.cv = 0, beta = 1.5, beta.cv = 10)

# simulate reference data
dr1.hv <- sim.dp(tp, model.par = par.r1.hv, seed = 100, plot = FALSE)

# model.par for test 
par.t1.hv <- list(fmax = 100, fmax.cv = 3, mdt = 12.29, mdt.cv = 15,
                  tlag = 0, tlag.cv = 0, beta = 1.727, beta.cv = 12)

# simulate test data with low variability
dt1.hv <- sim.dp(tp, model.par = par.t1.hv, seed = 100, plot = FALSE)
```

When CV conditions are not fulfilled, the function will stop running. 

```r
dt1.hv_vs_dr1.hv <- calcf2(dt1.hv$sim.disso, dr1.hv$sim.disso)
# The f2 method was applied according to EMA's BE guideline.
# 
# Individual data was provided with option 'cv.rule = TRUE',
# therefore, CV has been calculated and checked accordingly.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         19.28      18.00      24.90      25.61 
# 10        51.37      43.24      15.41      18.66 
# 15        76.99      64.80       8.66      13.46 
# 20        91.25      79.98       4.08       9.19 
# ------------------------------------------------ 
# 30        99.29      94.75       1.57       3.44 
# 45       100.00      99.54       1.53       1.45 
# 60       100.01      99.98       1.53       1.53 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV criteria not fulfilled; therefore, f2 method cannot be applied.
# Error in calcf2(dt1.hv$sim.disso, dr1.hv$sim.disso): 
# You should consider alternative methods such as bootstrap f2.
```

There might be a time that you just want to check the $f_2$ regardless
of the variability, option `cv.rule` can be set to `FALSE` in such case.

```r
dt1.hv_vs_dr1.hv <- calcf2(dt1.hv$sim.disso, dr1.hv$sim.disso,
                           cv.rule = FALSE)
# The f2 method was applied according to EMA's BE guideline.
# 
# Individual data was provided with option 'cv.rule = FALSE',
# therefore, CV has been calculated but has not been checked.
# You should really consider setting 'cv.rule = TRUE' to comply
# regulatory requirements.
# 
# Calculated mean and CV as follows:
# Time   Mean (T)   Mean (R)     CV (T)     CV (R) 
# 5         19.28      18.00      24.90      25.61 
# 10        51.37      43.24      15.41      18.66 
# 15        76.99      64.80       8.66      13.46 
# 20        91.25      79.98       4.08       9.19 
# ------------------------------------------------ 
# 30        99.29      94.75       1.57       3.44 
# 45       100.00      99.54       1.53       1.45 
# 60       100.01      99.98       1.53       1.53 
# ==================================
# Number of units for test is      : nt = 12
# Number of units for reference is : nr = 12
# 
# CV has not been checked.
# 
# The time points above the dashed line are used in f2 calculation.
# 
# Estimated f2 = 51.53
```

![](/home/zhengguo/github/bootf2/vignettes/calcf2_files/figure-html/calcf2-cvko-b-1.png)<!-- -->
