# test function sim.dp() -------------------------------------------------------

# time points
tp    <- c(5, 10, 15, 20, 30, 45, 60)
dp    <- c(5, 20, 42, 63, 88, 93, 96)
dp.cv <- c(19, 15, 10, 10, 8, 5, 3)
fo.par.err <- list(fmax = 100, fmax.cv = 0, tlag = 0, tlag.cv = 0,
                   k = 1, kcv = 10)

test_that("error checking", {
  expect_error(sim.dp(dp = dp, seed = 123, plot = FALSE),
               "Time points 'tp' must be provided.")
  expect_error(sim.dp(tp, dp = c(0, dp), seed = 123, plot = FALSE),
               "Length of 'tp' and 'dp' must be equal. Check your data.")
  expect_error(sim.dp(tp, dp, dp.cv = c(0, dp.cv), seed = 123, plot = FALSE),
               "Length of 'dp' and 'dp.cv' must be equal. Check your data.")
  expect_error(sim.dp(tp, dp, n.units = 6, seed = 123, plot = FALSE),
               paste0("To simulate dissolution profiles based on multivariate ",
                      "normal\ndistribution, 'n.units' should not be less ",
                      "than number of time\npoints 'tp'. Please either ",
                      "decrease the number of time points\nor increase ",
                      "'n.units'. Alternatively, you can use the approach\n",
                      "based on mathematical models."))
  expect_error(sim.dp(model = "first-order", model.par = fo.par.err,
                      seed = 123, plot = FALSE),
               paste0("Model parameters are incorrect. Three pairs of ",
                      "parameters are\nrequired: 'fmax/fmax.cv', 'k/k.cv', ",
                      "and 'tlag/tlag.cv'.\n"))
  expect_error(sim.dp(model.par = fo.par.err,
                      seed = 123, plot = FALSE),
               paste0("Model parameters are incorrect. Four pairs of ",
                      "parameters are\nrequired: 'fmax/fmax.cv', ",
                      "'tlag/tlag.cv', 'beta/beta.cv', and,\ndepending on the ",
                      "mathematical expression of the model, either\n",
                      "'alpha/alpha.cv' or 'mdt/mdt.cv'. Check your data.\n"))
  expect_warning(sim.dp(dp.cv = dp.cv, seed = 10, plot = FALSE, message = TRUE),
                 "'dp.cv' is given without 'dp', so it is ignored.\n")
})
