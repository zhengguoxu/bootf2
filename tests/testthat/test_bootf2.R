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
  expect_equal(format(res$boot.summary[res$boot.summary$f2.type ==
                                      "exp.f2", ]$f2o, digits = 10),
               "56.83163678")
  # mean of bootstrap expected f2
  expect_equal(format(res$boot.summary[res$boot.summary$f2.type ==
                                        "exp.f2", ]$boot.mean, digits = 10),
               "57.39924764")
  # normal interval lower and upper
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Normal", 3],
                     digits = 10), "47.74072198")
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Normal", 4],
                     digits = 10), "64.78732985")
  # basic interval lower and upper
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Basic", 3],
                     digits = 10), "46.61923347")
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Basic", 4],
                     digits = 10), "63.63610457")
  # percentile interval lower and upper
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Percentile (boot)", 3],
                     digits = 10), "50.02716898")
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                 res$boot.ci$ci.type == "Percentile (boot)", 4],
                     digits = 10), "67.04404008")
  # BCa (boot) interval lower and upper
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                   res$boot.ci$ci.type == "BCa (boot)", 3],
                     digits = 10), "49.82224368")
  expect_equal(format(res$boot.ci[res$boot.ci$f2.type == "exp.f2" &
                                   res$boot.ci$ci.type == "BCa (boot)", 4],
                     digits = 10), "66.69508423")
})


