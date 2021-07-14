#' Helper Functions
#'
#' @param x Numeric vector
#' @name helper
NULL

# Contents:
## - bpwhisker.l
## - bpwhisker.u
## - ci.header
## - mod.ref

#' @rdname helper
# upper and lower whisker for boxplot

bpwhisker.l <- function(x) {
  qt25 <- quantile(x, probs = 0.25, na.rm = TRUE)[[1]]
  qt75 <- quantile(x, probs = 0.75, na.rm = TRUE)[[1]]
  min(x[x >= qt25 - 1.5*(qt75 - qt25)], na.rm = TRUE)
}

#' @rdname helper
bpwhisker.u <- function(x) {
  qt25 <- quantile(x, probs = 0.25, na.rm = TRUE)[[1]]
  qt75 <- quantile(x, probs = 0.75, na.rm = TRUE)[[1]]
  max(x[x <= qt75 + 1.5*(qt75 - qt25)], na.rm = TRUE)
}

#' @param digits An integer indicating the decimal points for the output.
#' @rdname helper
ci.header <- function(digits) {
  # need some trick for print format depending on user set decimal
  # for CI type: "%23s", rest "%3s", lower, "%3s", upper
  # header ----
  if (digits <= 2) {
    cat("         Types of         Lower   Upper\n",
        "   Confidence Intervals   <----------->\n", sep = "")
  } else {
    cat("         Types of         Lower", rep(" ", 2*digits - 1),
        "Upper\n", "   Confidence Intervals   <----",
        rep("-", 2*digits - 1), "---->\n", sep = "")
  }
}

#' @param tp,ref.dp Numeric vector of time points \code{tp} and their
#'   corresponding mean dissolution profiles \code{ref.dp}.
#' @param digits An integer indicating the decimal points for the output.
#' @param model Strings of model names. Currently only 'Weibull' and
#'   'first-order' models are supported.
#' @param time.unit Character strings indicating the unit of time. It should
#'   be either \code{"min"} for minute or \code{"h"} for hour. It is mainly
#'   used for checking CV rules and making plot. @seealso [calcf2()].
#' @rdname helper
mod.ref <- function(tp, ref.dp, digits = 4, model, time.unit) {
  # add random noise in case mean values comes from simulation
  ref.dp <- ref.dp + rnorm(length(ref.dp), sd = 0.01)

  fmax.ini <- max(ref.dp) + 1

  # 50% of the 1st time points with non-zero dissolution
  tlag.ini <- 0.5*tp[ref.dp > 0][[1]]

  # for IR dosage form, 5 min not considered as lag time
  if (all(time.unit == "min", tlag.ini <= 5)) {
    tlag.ini <- 0
  }

  tp.adj <- tp - tlag.ini
  tp.adj[tp.adj < 0] <- 0

  tmp1 <- data.frame(tp = tp, ref.dp = ref.dp, tp.adj = tp.adj)

  if (model == "first-order") {
    mod.ini <- lm(log(1 - ref.dp/fmax.ini) ~ 0 + tp.adj, data = tmp1)

    k1.ini <- -coef(mod.ini)[[1]]

    ctr <- nls.control(maxiter = 500)

    mod <- minpack.lm::nlsLM(
      ref.dp ~ fmax*(1 - exp(-k1*(tp - tlag))), data = tmp1, control = ctr,
      start = c(fmax = fmax.ini, tlag = tlag.ini, k1 = k1.ini)
    )

    # coef(mod)[[2]] --> tlag,
    if (any(all(time.unit == "min", coef(mod)[[2]] <= 5),
            coef(mod)[[2]] <= 0)) {
      mod.par <- data.frame(model = "first-order",
                            fmax = coef(mod)[[1]], fmax.cv = 5,
                            tlag = 0, tlag.cv = 0,
                            k1 = coef(mod)[[3]], k1.cv = 30,
                            time.unit = time.unit, stringsAsFactors = FALSE)
    } else {
      mod.par <- data.frame(model = "first-order",
                            fmax  = coef(mod)[[1]], fmax.cv = 5,
                            tlag = coef(mod)[[2]], tlag.cv = 30,
                            k1 = coef(mod)[[3]], k1.cv = 30,
                            time.unit = time.unit, stringsAsFactors = FALSE)
    }
    # end model first-order
  } else if (model == "Weibull") {
    # mean dissolution time corresponds 0.63212 F/Fmax
    t1 <- max(tmp1[tmp1[, 2] <= 0.63212*fmax.ini, 3])
    t2 <- min(tmp1[tmp1[, 2] > 0.63212*fmax.ini, 3])
    d1 <- tmp1[tmp1[, 3] == t1, 2]
    d2 <- tmp1[tmp1[, 3] == t2, 2]
    mdt.ini <- t2 - (t2 - t1)*(d2 - 0.63212*fmax.ini)/(d2 - d1) + tlag.ini

    mdt.ini <- ifelse(mdt.ini < t1, t1, mdt.ini)

    mod.ini <- lm(log10(-log(1 - ref.dp/fmax.ini)) ~ 0 + log(tp.adj/mdt.ini),
                  data = tmp1[tmp1[, 3] > 0, ])

    beta.ini <- coef(mod.ini)[[1]]

    ctr <- nls.control(maxiter = 500)

    mod <- minpack.lm::nlsLM(
      ref.dp ~ fmax*(1 - exp(-((tp - tlag)/mdt)^beta)),
      data = tmp1[tmp1[, 3] > 0, ], control = ctr,
      start = c(fmax = fmax.ini, tlag = tlag.ini, mdt = mdt.ini,
                beta = beta.ini)
    )

    if (any(all(time.unit == "min", coef(mod)[[2]] <= 5),
            coef(mod)[[2]] <= 0)) {
      mod.par <- data.frame(model = "Weibull",
                            fmax = coef(mod)[[1]], fmax.cv = 5,
                            tlag = 0, tlag.cv = 0,
                            mdt = coef(mod)[[3]], mdt.cv = 30,
                            beta = coef(mod)[[4]], beta.cv = 30,
                            time.unit = time.unit, stringsAsFactors = FALSE)
    } else {
      mod.par <- data.frame(model = "Weibull",
                            fmax = coef(mod)[[1]], fmax.cv = 5,
                            tlag = coef(mod)[[2]], tlag.cv = 30,
                            mdt = coef(mod)[[3]],  mdt.cv = 30,
                            beta = coef(mod)[[4]], beta.cv = 30,
                            time.unit = time.unit, stringsAsFactors = FALSE)
    }
  }# end Weibull model fitting
  return(mod.par)
} # end function mod.ref



