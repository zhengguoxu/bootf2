test <- shah1998$test1
ref <- shah1998$ref

res <- bootf2(test, ref, n.boots = 1000L, min.points = 3L,
              print.report = FALSE, output.to.screen = FALSE)

# results are compared to output of Mendyk's bootf2BCA, v 1.3
# save results from bootf2BCa to res0, then
# all.equal(as.vector(res0$F2), res$boot.f2$exp.f2)
# shows all each individual f2 for each replicate are equal.

test_that("compare to bootf2BCA for expected f2", {
  # original expected f2
  expect_equal(round(res$boot.summary[res$boot.summary$f2.type ==
                                      "exp.f2", ]$f2o, digits = 2),
               56.83)
  # mean of bootstrap expected f2
  expect_equal(round(res$boot.summary[res$boot.summary$f2.type ==
                                        "exp.f2", ]$boot.mean, digits = 2),
               57.40)
  # normal interval lower and upper
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Normal", 3],
                     digits = 2), 47.74)
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Normal", 4],
                     digits = 2), 64.79)
  # basic interval lower and upper
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Basic", 3],
                     digits = 2), 46.62)
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Basic", 4],
                     digits = 2), 63.64)
  # percentile interval lower and upper
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Percentile (boot)", 3],
                     digits = 2), 50.03)
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Percentile (boot)", 4],
                     digits = 2), 67.04)
  # BCa (boot) interval lower and upper
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                   res$boot.ci$ci.type == "BCa (boot)", 3],
                     digits = 2), 49.82)
  expect_equal(round(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                   res$boot.ci$ci.type == "BCa (boot)", 4],
                     digits = 3), 66.695)
})


