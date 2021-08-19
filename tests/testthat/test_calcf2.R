# test function calcf2() -------------------------------------------------------

# published data ---------------------------------------------------------------
# load("data/shah1998.rda")
test_that("test data from Shah et al 1998", {
  expect_equal(round(calcf2(shah1998$test1, shah1998$ref,
                            regulation = "FDA", both.TR.85 = TRUE,
                            message = FALSE, plot = FALSE)$f2.value,
                     digits = 2), 60.03)
  expect_equal(round(calcf2(shah1998$test2, shah1998$ref,
                            regulation = "FDA", both.TR.85 = TRUE,
                            message = FALSE, plot = FALSE)$f2.value,
                     digits = 2), 51.08)
  expect_equal(round(calcf2(shah1998$test3, shah1998$ref,
                            regulation = "FDA", both.TR.85 = TRUE,
                            message = FALSE, plot = FALSE)$f2.value,
                     digits = 2), 51.19)
  expect_equal(round(calcf2(shah1998$test4, shah1998$ref,
                            regulation = "FDA", both.TR.85 = TRUE,
                            message = FALSE, plot = FALSE)$f2.value,
                     digits = 2), 50.07)
  expect_equal(round(calcf2(shah1998$test5, shah1998$ref,
                            regulation = "FDA", both.TR.85 = TRUE,
                            message = FALSE, plot = FALSE)$f2.value,
                     digits = 2), 48.05)
})


# low var ----------------------------------------------------------------------
# time points
tp <- c(5, 10, 15, 20, 30, 45, 60)

# model.par for reference with low variability
par.r1.lv <- list(fmax = 100, fmax.cv = 3, mdt = 15, mdt.cv = 14,
                  tlag = 0, tlag.cv = 0, beta = 1.5, beta.cv = 8)

# simulate reference data
dr1.lv <- sim.dp(tp, model.par = par.r1.lv, seed = 100, plot = FALSE)$sim.disso

# model.par for test
par.t1.lv <- list(fmax = 100, fmax.cv = 3, mdt = 12.29, mdt.cv = 12,
                  tlag = 0, tlag.cv = 0, beta = 1.727, beta.cv = 9)

# simulate test data with low variability
dt1.lv <- sim.dp(tp, model.par = par.t1.lv, seed = 100, plot = FALSE)$sim.disso

## test should be OK -------------------------------------------------
test_that("test calcf2 with CV OK", {
  f2ema <- calcf2(dt1.lv, dr1.lv, regulation = "EMA",
                  message = FALSE, plot = FALSE)
  f2fda <- calcf2(dt1.lv, dr1.lv, regulation = "FDA",
                  message = FALSE, plot = FALSE)
  f2fda2 <- calcf2(dt1.lv, dr1.lv, regulation = "FDA",
                  message = FALSE, plot = FALSE, both.TR.85 = TRUE)
  f2who <- calcf2(dt1.lv, dr1.lv, regulation = "WHO",
                  message = FALSE, plot = FALSE)
  expect_equal(round(f2ema$f2.value, digits = 1), 51.2)
  expect_equal(round(f2fda$f2.value, digits = 1), 51.2)
  expect_equal(round(f2fda2$f2.value, digits = 1), 53.0)
  expect_equal(round(f2who$f2.value, digits = 1), 53.0)
})

## test initial checking error -------------------------------------------------
test_that("test initial checking", {
  expect_error(calcf2(dt1.lv, regulation = "EMA", message = TRUE),
               "Both 'test' and 'ref' have to be specified.")
  expect_error(calcf2(dt1.lv, dr1.lv, regulation = "EMA", both.TR.85 = TRUE,
                      message = TRUE),
               paste0("'both.TR.85 = TRUE' is only applicable when ",
                      "'regulation = FDA'."))
  expect_error(calcf2(dt1.lv, dr1.lv, path.out = "weird/directory",
                      message = TRUE),
               "You should provided both 'path.out' and 'file.out'.")
  expect_error(calcf2(dt1.lv, dr1.lv, file.out = "weird.filename",
                      message = TRUE),
               "You should provided both 'path.out' and 'file.out'.")
  expect_error(calcf2(dt1.lv, dr1.lv, file.in = "weird.filename",
                      message = TRUE),
               paste0("Please provide the directory 'path.in' where ",
                      "the file is stored."))
  expect_error(calcf2(dt1.lv, dr1.lv, path.in = getwd(),message = TRUE),
               "Please provide the name of the data file.")
  expect_error(calcf2(dt1.lv, dr1.lv, path.in = "weird/directory",
                      file.in = "weird.filename", message = TRUE),
               paste0("The directory you specified does not exist. ",
                      "Check your spelling."))
  expect_error(calcf2(dt1.lv, dr1.lv, path.in = getwd(),
                      file.in = "weird.filename", message = TRUE),
               paste0("The file you specified does not exist. Don't ",
                      "forget to include\nthe extension 'xlsx' or 'xls' ",
                      "in the file name."))
  expect_error(calcf2(dt1.lv, dr1.lv, path.out = "weird/directory",
                      file.out = "weird.filename", message = TRUE),
               paste0("The directory you specified does not exist. ",
                      "Check your spelling."))
})

## test unequal time points ------------------------------------------
test_that("test unequal time points", {
  expect_error(calcf2(dt1.lv[-2, ], dr1.lv[-3, ], message = TRUE),
               paste0("Time points of two data sets should be same! ",
                      "Please check your data.\n\n"))
})


# high var ---------------------------------------------------------------------
# model.par for reference with high variability
par.r1.hv <- list(fmax = 100, fmax.cv = 3, mdt = 15, mdt.cv = 20,
                  tlag = 0, tlag.cv = 0, beta = 1.5, beta.cv = 10)

# simulate reference data
dr1.hv <- sim.dp(tp, model.par = par.r1.hv, seed = 100, plot = FALSE)$sim.disso

# model.par for test
par.t1.hv <- list(fmax = 100, fmax.cv = 3, mdt = 12.29, mdt.cv = 15,
                  tlag = 0, tlag.cv = 0, beta = 1.727, beta.cv = 12)

# simulate test data with low variability
dt1.hv <- sim.dp(tp, model.par = par.t1.hv, seed = 100, plot = FALSE)$sim.disso


test_that("test calcf2 with CV not OK", {
  expect_error(calcf2(dt1.hv, dr1.hv, regulation = "EMA", plot = FALSE,
                      message = TRUE),
               paste0("You should consider alternative methods such as ",
                      "bootstrap f2.\n\n"))
  expect_warning(calcf2(dt1.hv, dr1.hv, regulation = "FDA", plot = FALSE,
                        message = TRUE),
                 paste0("f2 was calculated while CV criterion is not strictly ",
                        "fulfilled; you \nmight want to consider an alternative ",
                        "method such as bootstrap f2.\n\n"))
})

# fast dissolution -------------------------------------------------------------
tp1 <- c(10, 15, 20)
# model.par for test, fast
par.t1.fast <- list(fmax = 100, fmax.cv = 3, mdt = 8, mdt.cv = 25,
                    tlag = 0, tlag.cv = 0, beta = 1.727, beta.cv = 12)

# simulate test data with low variability
dt1.fast <- sim.dp(tp1, model.par = par.t1.fast, seed = 100,
                   plot = FALSE)$sim.disso

# model.par for ref, fast
par.r1.fast <- list(fmax = 100, fmax.cv = 3, mdt = 12, mdt.cv = 25,
                    tlag = 0, tlag.cv = 0, beta = 1.5, beta.cv = 10)

# simulate test data with low variability
dr1.fast <- sim.dp(tp1, model.par = par.r1.fast, seed = 100,
                   plot = FALSE)$sim.disso
test_that("test fast dissolution",{
  expect_warning(calcf2(dt1.fast, dr1.fast, cv.rule = FALSE,
                        min.points = 2, plot = FALSE, message = TRUE),
                 paste0("Warning: f2 was calculated with less than 3 time ",
                        "points for information\npurpose only. You should ",
                        "add more earlier points in your dissolution method."))
})

