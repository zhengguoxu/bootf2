tp <- c(5, 10, 15, 20, 30, 45, 60)

mod.par.r <- list(fmax = 100, fmax.cv = 2, tlag = 0, tlag.cv = 0,
                  mdt = 25, mdt.cv = 4, beta = 2.1, beta.cv = 3)

d.r <- sim.dp(tp, model.par = mod.par.r, seed = 100, n.units = 120L,
              plot = FALSE)
dp <- d.r$sim.summary$dp

test_that("error check", {
  expect_error(sim.dp.byf2(tp, dp, plot = FALSE),
               "You have to provide 'target.f2'.")
  expect_error(sim.dp.byf2(tp, dp, target.f2 = "50", plot = FALSE),
               "'target.f2' has to be numeric.")
  expect_error(sim.dp.byf2(tp, dp, target.f2 = c(50, 60, 70), plot = FALSE),
               paste0("'target.f2' should be a single number or a vector ",
                      "of 2 numbers."))
  expect_error(sim.dp.byf2(tp, sim.dp.out = d.r$sim.summary,
                           target.f2 = c(50, 51), plot = FALSE),
               paste0("The input data 'sim.dp.out' is incorrect. Use either ",
                      "the output\nof 'sim.dp' function or provide 'tp' and ",
                      "'dp'.\n"))
  expect_error(sim.dp.byf2(tp, target.f2 = 50, plot = FALSE),
               paste0("The imput data 'sim.dp.out' is missing. In this case, ",
                      "you need to\nprovide time point 'tp' and mean ",
                      "dissolution profile 'dp'."))
  expect_error(sim.dp.byf2(tp, dp = dp[-c(2,3)], target.f2 = 50, plot = FALSE),
               "The length of 'tp' should be equal to that of 'dp'.")
})

dt55 <- sim.dp.byf2(tp, dp[-1], target.f2 = 55, seed = 123, plot = FALSE)

test_that("check results", {
  expect_equal(format(dt55$model.par$fmax, digits = 10), "100.0120726")
  expect_equal(format(dt55$model.par$mdt, digits = 10), "30.1480993")
  expect_equal(format(dt55$model.par$beta, digits = 10), "1.675170536")
  expect_equal(format(dt55$model.par$f2, digits = 10), "54.99999963")
})


