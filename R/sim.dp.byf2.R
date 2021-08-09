#' Simulate Dissolution Profiles by \eqn{f_2}{f2} Values
#'
#' Given any mean dissolution profile `dp`, this function will simulate a mean
#' dissolution profile such that when the newly simulated profile is compared
#' to `dp`, the calculated \eqn{f_2}{f2} will be equal to the predefined target
#' \eqn{f_2}{f2} value.
#'
#' @importFrom stats lm coef nls.control optim
#' @importFrom minpack.lm nlsLM
#' @usage
#' sim.dp.byf2(tp, dp, target.f2, seed = NULL, min.points = 3L,
#'             regulation = c("EMA", "FDA", "WHO", "Canada", "ANVISA"),
#'             model = c("Weibull", "first-order"), digits = 2L,
#'             max.disso = 100, message = TRUE, both.TR.85 = FALSE,
#'             time.unit = c("min", "h"), plot = TRUE, sim.dp.out,
#'             sim.target = c("ref.pop", "ref.median", "ref.mean"),
#'             model.par.cv = 50, fix.fmax.cv = 0, random.factor = 3)
#'
#' @param tp,dp *Numeric* vector of time points `tp` and the mean dissolution
#'   profiles `dp`.
#' @param target.f2 *Numeric* value of target \eqn{f_2}{f2}. It can be a
#'   *single value*, or a *vector of two values* that represent the lower and
#'   upper limit of target \eqn{f_2}{f2} value. See Details.
#' @param seed *Integer* seed value for reproducibility. If missing, a random
#'   seed will be generated for reproducibility purpose.
#' @param min.points An *integer* indicating the minimum time points to be used
#'   to calculate \eqn{f_2}{f2}. The default value 3 should be used for
#'   conventional \eqn{f_2}{f2} calculation. See Details. @seealso [calcf2()].
#' @param regulation *Character* strings indicating regulatory guidelines.
#'   @seealso [calcf2()] for the discussion of those guidelines.
#' @param model *Character* strings of model names. Currently only `"Weibull"`
#'   and `"first-order"` models are supported. @seealso [sim.dp()] for the
#'   description of models.
#' @param digits An *integer* indicating the decimal points for the output.
#' @param max.disso *Numeric* value for the maximum possible dissolution.
#'   In theory, the maximum dissolution is 100%, but in practice, it is not
#'   uncommon to see higher values, such as 105%, or much lower values,
#'   especially for products with low solubility.
#' @param message *Logical*. If `TRUE`, basic information of the simulation
#'   will be printed on screen.
#' @param both.TR.85 *Logical*. If `TRUE`, and if `regulation = "FDA"`, all
#'   measurements up to the time points at which both test and reference
#'   products dissolve more than 85% will be used to calculate \eqn{f_2}{f2}.
#'   This is the conventional, but incorrect, interpretation of the US FDA rule.
#'   Therefore, the argument should only be set to `TRUE` for validation purpose
#'   such as comparing the results from old literatures that use the wrong
#'   interpretation to calculate \eqn{f_2}{f2}. @seealso [calcf2()] for detailed
#'   discussion.
#' @param time.unit *Character* strings indicating the unit of time. It should
#'   be either `"min"` for minute or `"h"` for hour. It is mainly used for
#'   checking CV rules and making plot. See Regulation in Details.
#'   @seealso [calcf2()].
#' @param plot *Logical*. If `TRUE`, a a dissolution versus time plot will be
#'   included in the output.
#' @param sim.dp.out Output of function `sim.dp()`. If this argument is
#'   supplied by the user, then `tp/dp`, `regulation`, `model`, `max.disso`,
#'   and `time.unit` will be ignored, if they are provided by the user, since
#'   all those arguments are included in `sim.dp.out`.
#' @param sim.target *Character* strings indicating to which target profile
#'   should the newly simulated profile be compared for the calculation of
#'   \eqn{f_2}{f2}. This argument is only applicable when `sim.dp.out` is
#'   provided. See Details.
#' @param  model.par.cv,fix.fmax.cv *Numeric* values expressed as percentages
#'   used for random generation of model parameters and fmax when optimization
#'   algorithm is not used, i.e., when `target.f2` is a vector of two numbers.
#'   See Details.
#' @param random.factor *Numeric* value used for random generation of model
#'   parameters when optimization algorithm is used, i.e., when `target.f2`
#'   is a single number. See Details.
#'
#' @return A *list* of 2 components: a *data frame of model parameters* and a
#'   *data frame of mean dissolution profile* generated using the said model
#'   parameters. The output can be passed to function `sim.dp()` to further
#'   simulate multiple individual profiles.
#'
#' @details
#'
#' The main principle of the function is as follows:
#' 1. For any given mean dissolution profile `dp`, fit a suitable mathematical
#'    model and obtain model parameters.
#'    - No precise fitting is required since those parameters will be served as
#'      *initial value* for further model fitting.
#'    - If `sim.dp.out`, the output of the function `sim.dp()`, is available,
#'      no initial fitting is necessary as model parameters can be read directly
#'      from the output, unless multivariate normal distribution approach was
#'      used in `sim.dp()`. In such case, initial model fitting will be done.
#' 2. Find a suitable model parameters and simulate a new dissolution profile,
#'    comparing the new profile to the provided reference profile `dp` by
#'    calculating \eqn{f_2}{f2}. If the the obtained \eqn{f_2}{f2} is
#'    *equal to*, or *within the lower and upper limit of*, the `target.f2`,
#'    then the function will output the obtained model parameters and the new
#'    profile.
#'
#' There are two approaches used to find the suitable model parameters:
#' - If `target.f2` is a single value, optimization algorithm will be used and
#'   the newly simulated dissolution profile will have \eqn{f_2}{f2} equal to
#'   `target.f2` when compared to `dp` (within numeric precision defined by the
#'   tolerance).
#' - If `target.f2` is a vector of two numbers representing the lower and upper
#'   limit of target \eqn{f_2}{f2} value, such as `target.f2 = c(lower, upper)`,
#'   then dissolution will be obtained by random searching and the calculated
#'   \eqn{f_2}{f2} will be within the range of lower and upper.
#'
#' For example, you can set `target.f2 = c(54.95, 55.04)` to have target
#' \eqn{f_2}{f2} of 55. Since \eqn{f_2}{f2} should be normally reported without
#' decimal, in practice, this precision is enough. You might be able to do with
#' `c(54.99, 55.01)` if you really need more precision. However, the narrower
#' the range, the longer it takes the function to run. With narrow range such as
#' `c(54.999, 55.001)`, it is likely the program will fail. In such case,
#' provide single value to use optimization algorithm.
#'
#' Arguments `model.par.cv`, `fix.fmax.cv`, and `random.factor` are certain
#' numeric values used to better control the random generation of model
#' parameters. The default values should work in most scenarios. Those values
#' should be changed only when the function failed to return any value. Read
#' vignette of the function (`vignette("sim.dp.byf2", package = "bootf2")`)
#' for more details.
#'
#' The data frame `sim.summary` in `sim.dp.out`, the output of function
#' `sim.dp()`, contains `dp`, the population profile, and `sim.mean` and
#' `sim.median`, the mean and median profiles calculated with `n.units` of
#' simulated individual profiles. All these three profiles could be used as the
#' target profile that the newly simulated profile will be compare to. Argument
#' `sim.target` defines which of the three will be used: `ref.pop`, `ref.mean`,
#' and `ref.median` correspond to `dp`, `sim.mean` and `sim.median`,
#' respectively.
#'
#' @examples
#' tp <- c(5, 10, 15, 20, 30, 45, 60)
#'
#' mod.par.r <- list(fmax = 100, fmax.cv = 2, tlag = 0, tlag.cv = 0,
#'                   mdt = 25, mdt.cv = 4, beta = 2.1, beta.cv = 3)
#'
#' d.r <- sim.dp(tp, model.par = mod.par.r, seed = 100, n.units = 120L,
#'               plot = FALSE)
#'
#' model.par1 <- sim.dp.byf2(sim.dp.out = d.r, target.f2 = 60, seed = 123)
#' model.par2 <- sim.dp.byf2(sim.dp.out = d.r, target.f2 = c(59.95, 60.04),
#'                           seed = 123)
#'
#' @export
sim.dp.byf2 <- function(tp, dp, target.f2, seed = NULL, min.points = 3L,
                        regulation = c("EMA", "FDA", "WHO", "Canada", "ANVISA"),
                        model = c("Weibull", "first-order"), digits = 2L,
                        max.disso = 100, message = TRUE, both.TR.85 = FALSE,
                        time.unit = c("min", "h"), plot = TRUE, sim.dp.out,
                        sim.target = c("ref.pop", "ref.median", "ref.mean"),
                        model.par.cv = 50, fix.fmax.cv = 0, random.factor = 3) {
  # initial check --------------------------------------------------------------
  regulation <- match.arg(regulation)
  model      <- match.arg(model)
  time.unit  <- match.arg(time.unit)
  sim.target <- match.arg(sim.target)

  # if no seed is provided, generate one
  if (is.null(seed)) {
    seed <- sample(1:(.Machine$integer.max - 1), 1)
  }
  set.seed(seed)

  if (missing(target.f2)) {
    stop("You have to provide 'target.f2'.")
  } else {# not missing
    if (isFALSE(is.numeric(target.f2))) {
      stop("'target.f2' has to be numeric.")
    } else {# numeric value
      if(length(target.f2) == 1L) {
        optimization <- TRUE
      } else if (length(target.f2) == 2L) {
        optimization <- FALSE
        f2.lower     <- min(target.f2)
        f2.upper     <- max(target.f2)
      } else {
        stop("'target.f2' should be a single number or a vector of 2 numbers.")
      }
    }
  }# end target.f2 not missing

  # get model parameters for ref------------------------------------------------
  # need tp + dp or sim.dp.out but not both. former: regression; latter: read
  if (!missing(sim.dp.out)) {# read info from sim.dp.out
    if (all(attr(sim.dp.out, "come.from") != "sim.dp",
            any(missing(tp), missing(dp)))) {
      stop("The input data 'sim.dp.out' is incorrect. Use either the output\n",
           "of 'sim.dp' function or provide 'tp' and 'dp'.\n")
    }

    if (all(isTRUE(message), any(!missing(tp), !missing(dp)))) {
      cat("The output of the function 'sim.dp' was provided as input data,\n",
          "therefore, 'tp'/'dp' will be ignored.\n\n", sep = "")
    }

    # get ref data for calcf2() function
    tp0        <- sim.dp.out$sim.summary$time
    tp         <- tp0[tp0 > 0]
    ref.dp     <- sim.dp.out$sim.summary$dp[tp0 > 0]
    ref.median <- sim.dp.out$sim.summary$sim.median[tp0 > 0]
    ref.mean   <- sim.dp.out$sim.summary$sim.mean[tp0 > 0]
    max.disso  <- sim.dp.out$sim.info$max.disso
    time.unit  <- sim.dp.out$sim.info$time.unit

    if (sim.target == "ref.pop") {
      data.r <- data.frame(tp = tp, dp = ref.dp, stringsAsFactors = FALSE)
    } else if (sim.target == "ref.mean") {
      data.r <- data.frame(tp = tp, dp = ref.mean, stringsAsFactors = FALSE)
    } else {
      data.r <- data.frame(tp = tp, dp = ref.median, stringsAsFactors = FALSE)
    }

    # for optimization later
    if (!is.na(sim.dp.out$sim.info$model)) {
      fmax.t <- sim.dp.out$sim.info$fmax
      tlag.t <- sim.dp.out$sim.info$tlag
      tadj.t <- tp - tlag.t
      tadj.t[tadj.t < 0] <- 0
    }

    if (sim.dp.out$sim.info$model == "Weibull") {# Weibull
      if ("mdt" %in% names(sim.dp.out$sim.info)) {# Weibull with mdt----
        mod.par <- data.frame(model     = sim.dp.out$sim.info$model,
                              fmax      = sim.dp.out$sim.info$fmax,
                              fmax.cv   = sim.dp.out$sim.info$fmax.cv,
                              tlag      = sim.dp.out$sim.info$tlag,
                              tlag.cv   = sim.dp.out$sim.info$tlag.cv,
                              mdt       = sim.dp.out$sim.info$mdt,
                              mdt.cv    = sim.dp.out$sim.info$mdt.cv,
                              beta      = sim.dp.out$sim.info$beta,
                              beta.cv   = sim.dp.out$sim.info$beta.cv,
                              time.unit = sim.dp.out$sim.info$time.unit,
                              stringsAsFactors = FALSE)
      } else {# Weibull with alpha ----------------------------------------
        mod.par <- data.frame(model     = sim.dp.out$sim.info$model,
                              fmax      = sim.dp.out$sim.info$fmax,
                              fmax.cv   = sim.dp.out$sim.info$fmax.cv,
                              tlag      = sim.dp.out$sim.info$tlag,
                              tlag.cv   = sim.dp.out$sim.info$tlag.cv,
                              alpha     = sim.dp.out$sim.info$alpha,
                              alpha.cv  = sim.dp.out$sim.info$alpha.cv,
                              beta      = sim.dp.out$sim.info$beta,
                              beta.cv   = sim.dp.out$sim.info$beta.cv,
                              time.unit = sim.dp.out$sim.info$time.unit,
                              stringsAsFactors = FALSE)
      }# end Weibull with alpha
    } else if (sim.dp.out$sim.info$model == "first-order") {# first-order ----
      mod.par <- data.frame(model     = sim.dp.out$sim.info$model,
                            fmax      = sim.dp.out$sim.info$fmax,
                            fmax.cv   = sim.dp.out$sim.info$fmax.cv,
                            tlag      = sim.dp.out$sim.info$tlag,
                            tlag.cv   = sim.dp.out$sim.info$tlag.cv,
                            k         = sim.dp.out$sim.info$k,
                            k.cv      = sim.dp.out$sim.info$k.cv,
                            time.unit = sim.dp.out$sim.info$time.unit,
                            stringsAsFactors = FALSE)
    } else {# model == NA, simulated by multivariate normal distribution
      # use regression to get mod.par in this case. 'model' from function input.
      mod.par <- mod.ref(tp = tp, ref.dp = data.r$dp, digits = digits,
                         model = model, max.disso = max.disso,
                         time.unit = time.unit)
    }# end model == NA in sim.do.out
  } else {# for missing sim.dp.out, regression to get initial mod.par ---------
    if (any(missing(tp), missing(dp))) {
      stop("The imput data 'sim.dp.out' is missing. In this case, you need to\n",
           "provide time point 'tp' and mean dissolution profile 'dp'.")
    }

    if (length(tp) != length(dp)) {
      stop("The length of 'tp' should be equal to that of 'dp'.")
    }

    # remove time 0 if any
    ref.dp <- dp[tp > 0]
    tp     <- tp[tp > 0]
    data.r <- data.frame(tp = tp, dp = ref.dp, stringsAsFactors = FALSE)

    # function for modelling ref.pop to get parameters # add here mod.ref
    mod.par <- mod.ref(tp = data.r$tp, ref.dp = data.r$dp, digits = digits,
                       model = model, max.disso = max.disso,
                       time.unit = time.unit)


    fmax.t <- mod.par$fmax
    tlag.t <- mod.par$tlag
    tadj.t <- tp - tlag.t
    tadj.t[tadj.t < 0] <- 0
  }# end missing sim.dp.out

  # having mod.par for all scenarios, now find model.par for test---------------
  if (isTRUE(optimization)) {#use stats::optim
    if (mod.par$model == "Weibull") {# Weibull optim
      if ("mdt" %in% names(mod.par)) {# Weibull optim for mdt ----
        # model param for optim function
        mpar <- c(
          mdt = mod.par$mdt*exp(rnorm(1, 0, mod.par$mdt.cv*random.factor/100)),
          beta = mod.par$beta*exp(rnorm(1, 0, mod.par$beta.cv*random.factor/100))
        )

        opt.weibull.mdt <- function(mpar, target.f2, ...) {
          tmp.t  <- fmax.t*(1 - exp(-(tadj.t/mpar[["mdt"]])^mpar[["beta"]]))
          data.t <- data.frame(tp = tp, dp = tmp.t, stringsAsFactors = FALSE)
          tmp.f2 <- calcf2(data.t, data.r, regulation = regulation,
                           cv.rule = FALSE, min.points = min.points,
                           both.TR.85 = both.TR.85, message = FALSE,
                           plot = FALSE, time.unit = time.unit)
          return(abs(tmp.f2$f2.value - target.f2))
        }

        res.f2 <- optim(mpar, opt.weibull.mdt, target.f2 = target.f2,
                        na.rm = TRUE)

        # mean of test that has f2 = target.f2
        dp.2 <-
          fmax.t*(1 - exp(-(tadj.t/res.f2$par[["mdt"]])^res.f2$par[["beta"]]))
        dp.2[!is.finite(dp.2)] <- 0
        dp.2[dp.2 < 0] <- 0

        data.t <- data.frame(tp = tp, dp = dp.2, stringsAsFactors = FALSE)
        f2.tmp <- calcf2(data.t, data.r, regulation = regulation,
                         cv.rule = FALSE, min.points = min.points,
                         both.TR.85 = both.TR.85, message = FALSE,
                         plot = FALSE, time.unit = time.unit)

        model.par <- data.frame(model = "Weibull", seed = seed,
                                fmax = fmax.t, tlag = tlag.t,
                                mdt = res.f2$par[["mdt"]],
                                beta = res.f2$par[["beta"]],
                                f2 = f2.tmp$f2.value,
                                f2.tp = f2.tmp$f2.tp,
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      } else {# Weibull optim for alpha ----
        mpar <- c(
          alpha = mod.par$alpha*exp(rnorm(1, 0, mod.par$alpha.cv*random.factor/100)),
          beta = mod.par$beta*exp(rnorm(1, 0, mod.par$beta.cv*random.factor/100))
        )

        opt.weibull.alpha <- function(mpar, target.f2, ...) {
          tmp.t  <- fmax.t*(1 - exp(-(tadj.t^mpar[["beta"]])/mpar[["alpha"]]))
          data.t <- data.frame(tp = tp, dp = tmp.t, stringsAsFactors = FALSE)
          tmp.f2 <- calcf2(data.t, data.r, regulation = regulation,
                           cv.rule = FALSE, min.points = min.points,
                           both.TR.85 = both.TR.85, message = FALSE,
                           plot = FALSE, time.unit = time.unit)
          return(abs(tmp.f2$f2.value - target.f2))
        }

        res.f2 <- optim(mpar, opt.weibull.alpha, target.f2 = target.f2,
                        na.rm = TRUE)

        dp.2 <- fmax.t*
          (1 - exp(-(tadj.t^res.f2$par[["beta"]])/res.f2$par[["alpha"]]))
        dp.2[!is.finite(dp.2)] <- 0
        dp.2[dp.2 < 0] <- 0

        data.t <- data.frame(tp = tp, dp = dp.2, stringsAsFactors = FALSE)
        f2.tmp <- calcf2(data.t, data.r, regulation = regulation,
                         cv.rule = FALSE, min.points = min.points,
                         both.TR.85 = both.TR.85, message = FALSE,
                         plot = FALSE, time.unit = time.unit)

        model.par <- data.frame(model = "Weibull", seed = seed,
                                fmax = fmax.t, tlag = tlag.t,
                                alpha = res.f2$par[["alpha"]],
                                beta = res.f2$par[["beta"]],
                                f2 = f2.tmp$f2.value,
                                f2.tp = f2.tmp$f2.tp,
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      }
    } else {# first-order optim
      mpar <- c(
        k = mod.par$k*exp(rnorm(1, 0, mod.par$k.cv*random.factor/100)),
        k.cv = mod.par$k.cv
      )

      opt.first.order <- function(mpar, target.f2, ...) {
        tmp.t  <- fmax.t*(1 - exp(-mpar[["k"]]*tadj.t))
        data.t <- data.frame(tp = tp, dp = tmp.t, stringsAsFactors = FALSE)
        tmp.f2 <- calcf2(data.t, data.r, regulation = regulation,
                         cv.rule = FALSE, min.points = min.points,
                         both.TR.85 = both.TR.85, message = FALSE,
                         plot = FALSE, time.unit = time.unit)
        return(abs(tmp.f2$f2.value - target.f2))
      }

      res.f2 <- optim(mpar, opt.first.order, target.f2 = target.f2,
                      na.rm = TRUE)

      dp.2 <- fmax.t*(1 - exp(-res.f2$par[["k"]]*tadj.t))
      dp.2[!is.finite(dp.2)] <- 0
      dp.2[dp.2 < 0] <- 0

      data.t <- data.frame(tp = tp, dp = dp.2, stringsAsFactors = FALSE)
      f2.tmp <- calcf2(data.t, data.r, regulation = regulation,
                       cv.rule = FALSE, min.points = min.points,
                       both.TR.85 = both.TR.85, message = FALSE,
                       plot = FALSE, time.unit = time.unit)

      model.par <- data.frame(model = "first-order", seed = seed,
                              fmax = fmax.t, tlag = tlag.t,
                              k = res.f2$par[["k"]],
                              f2 = f2.tmp$f2.value,
                              f2.tp = f2.tmp$f2.tp,
                              regulation = regulation,
                              stringsAsFactors = FALSE)
    }# end optim first-order model
  } else {# random search if not optimization
    repeat{
      fmax.2 <- mod.par$fmax*exp(rnorm(1, 0, fix.fmax.cv/100))
      fmax.2 <- ifelse(fmax.2 > max.disso, max.disso, fmax.2)
      tlag.2 <- mod.par$tlag*exp(rnorm(1, 0, model.par.cv/100))
      if (all(time.unit == "min", tlag.2 <= 5)) tlag.2 <- 0

      tp.adj2 <- tp - tlag.2
      tp.adj2[tp.adj2 < 0] <- 0

      if (mod.par$model == "Weibull") {
        beta.2 <- mod.par$beta*exp(rnorm(1, 0, model.par.cv/100))
        if ("mdt" %in% names(mod.par)) {
          mdt.2 <- mod.par$mdt*exp(rnorm(1, 0, model.par.cv/100))
          dp.2  <- fmax.2*(1 - exp(-(tp.adj2/mdt.2)^beta.2))
        } else {
          alpha.2 <- mod.par$alpha*exp(rnorm(1, 0, model.par.cv/100))
          dp.2    <- fmax.2*(1 - exp(-(tp.adj2^beta.2)/alpha.2))
        }
      } else {# first-order
        k.2 <- mod.par$k*exp(rnorm(1, 0, model.par.cv/100))
        dp.2 <- fmax.2*(1 - exp(-k.2*tp.adj2))
      }

      dp.2[!is.finite(dp.2)] <- 0
      dp.2[dp.2 < 0] <- 0
      data.t <- data.frame(tp = tp, dp = dp.2, stringsAsFactors = FALSE)

      # get f2 and compare to criteria ----------------------------------
      f2.tmp <- calcf2(data.t, data.r, regulation = regulation,
                       cv.rule = FALSE, digits = digits, message = FALSE,
                       min.points = min.points, both.TR.85 = both.TR.85,
                       plot = FALSE, time.unit = time.unit)
      if (all(is.finite(f2.tmp$f2.value), f2.tmp$f2.value >= f2.lower,
              f2.tmp$f2.value <= f2.upper)) break
    }# end repeat

    if (mod.par$model == "Weibull") {
      if ("mdt" %in% names(mod.par)) {
        model.par <- data.frame(model = "Weibull", seed = seed, fmax = fmax.2,
                                tlag = tlag.2, mdt = mdt.2, beta = beta.2,
                                f2 = f2.tmp$f2.value, f2.tp = f2.tmp$f2.tp,
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      } else {
        model.par <- data.frame(model = "Weibull", seed = seed, fmax = fmax.2,
                                tlag = tlag.2, alpha = alpha.2, beta = beta.2,
                                f2 = f2.tmp$f2.value, f2.tp = f2.tmp$f2.tp,
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      }
    } else {# first-order
      model.par <- data.frame(model = "first-order", seed = seed, fmax = fmax.2,
                              tlag = tlag.2, k = k.2, f2 = f2.tmp$f2.value,
                              f2.tp = f2.tmp$f2.tp, regulation = regulation,
                              stringsAsFactors = FALSE)
    }
  }# end random search

  # output population dissolution data
  diff.tr <- dp.2 - data.r$dp
  sim.disso <- data.frame(time = c(0, data.r$tp), ref = c(0, data.r$dp),
                          test = c(0, dp.2), diff.tr = c(0, diff.tr),
                          stringsAsFactors = FALSE)

  # plot -----------------------------------------------------------------------
  if (isTRUE(plot)) {
    # need this to get rid of "no visible binding for global variable" notes:
    time <- ref <- test <- NULL

    # RColorBrewer::brewer.pal(12, "Paired")
    # 2 "#1F78B4" and 6 "#E31A1C" for T and R
    cols <- c("R" = "#1F78B4", "T" = "#E31A1C")

    y.tick.max <- ceiling(max(sim.disso$ref, sim.disso$test))

    p1 <- ggplot(sim.disso, aes(x = time)) +
      geom_hline(yintercept = 85, size = 1.2, color = "gray",
                 alpha = 0.7, linetype = "dashed") +
      geom_point(aes(y = ref, color = "R"), size = 1.5) +
      geom_line(aes(y = ref, color = "R"), size =1.1) +
      geom_point(aes(y = test, color = "T"), size = 1.5) +
      geom_line(aes(y = test, color = "T"),size =1.1) +
      theme_bw() +
      theme(legend.title = element_blank(), legend.justification = c(1, 0),
            legend.position = c(0.98, 0.02)) +
      scale_x_continuous(limits = c(min(sim.disso$time), max(sim.disso$time)),
                         breaks = sim.disso$time) +
      scale_y_continuous(breaks = seq(0, y.tick.max + 5, 5),
                         limits = c(0, y.tick.max + 5)) +
      scale_colour_manual(values = cols) +
      ylab("Dissolution (%)") +
      xlab(paste0("Time (", time.unit, ")"))+
      ggtitle(paste0("Simulated Population Dissolution Profiles with f2 = ",
                     round(f2.tmp$f2.value, digits)),
              subtitle = paste0("Blue is the target profile (reference) and ",
                                "red is the simulated profile (test) with ",
                                "the predefined f2."))
    print(p1)

    sim.out <- list(model.par = model.par, sim.disso = sim.disso,
                    sim.plot = p1)
  } else sim.out <- list(model.par = model.par, sim.disso = sim.disso)

  if (isTRUE(message)) {
    cat("Obtained model parameters and calculated f2 are:\n")
    print(sim.out$model.par)

    cat("\nAnd the difference between simulated test and reference is:\n")
    print(sim.out$sim.disso)
  }

  return(sim.out)
}
