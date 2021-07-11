#' Simulate Dissolution Profiles by f2 Values
#'
#' Given any mean dissolution profile \code{dp}, this function will simulate a
#' mean dissolution profile such that when the newly simulated profile is
#' compared to \code{dp}, the calculated f2 will be equal to the predefined
#' target f2 value (within certain numeric precision).
#'
#' @importFrom stats lm coef nls.control optim
#' @importFrom minpack.lm nlsLM
#' @usage
#' sim.dp.byf2(tp, dp, sim.dp.out, target.f2, seed,
#'             regulation = c("EMA", "FDA", "WHO"),
#'             model = c("Weibull", "first-order"), digits = 2L,
#'             min.points = 3L, both.TR.85 = FALSE, max.disso = 105,
#'             model.par.cv = 50, fix.fmax.cv = 0, random.factor = 3,
#'             sim.target = c("ref.pop", "ref.median", "ref.mean"),
#'             time.unit = c("min", "h"), message = TRUE, plot = TRUE)
#'
#' @param tp,dp Numeric vector of time points \code{tp} and their corresponding
#'   mean dissolution profiles \code{dp}.
#' @param sim.dp.out Output of function \code{sim.dp}. If this parameter is
#'   supplied by users, then \code{tp/dp} will be ignored if they are provided.
#' @param target.f2 Numeric value of target f2. It can be a single value, or a
#'   vector of two values representing lower and upper limit of target f2 value.
#'   See Details.
#' @param seed Numeric seed value for reproducibility. If missing, a random
#'   seed will be generated for reproducibility purpose.
#' @param regulation Character strings indicating regulatory guidelines.
#'   @seealso [calcf2()].
#' @param model Strings of model names. Currently only 'Weibull' and
#'   'first-order' models are supported.
#' @param digits An integer indicating the decimal points for the output.
#' @param min.points An integer indicating the minimum time points to be used
#'   to calculate f2. The default is 3 for conventional f2 calculation.
#'   This parameter is mainly used for bootstrap f2 method. @seealso [calcf2()].
#' @param both.TR.85 A logic value indicating if the old (incorrect)
#'   interpretation (that both the test and reference should release more than
#'   85%) should be used for f2 calculation when \code{regulation = 'FDA'}.
#'   @seealso [calcf2()].
#' @param max.disso Numeric value for the maximum possible dissolution.
#'   In theory, the maximum dissolution should be 100%, but in practice,
#'   it is not uncommon to see higher values, such as 105%, or even 110%.
#' @param  model.par.cv,fix.fmax.cv Numeric value expressed as percentage.
#'   Used for random generation of model parameters and fmax when optimization
#'   algorithm is not used, i.e., when \code{target.f2} is a vector of two
#'   numbers. See Details.
#' @param random.factor Numeric value used for random generation of model
#'   parameters when optimization algorithm is used, i.e., when \code{target.f2}
#'   is a single number. See Details.
#' @param sim.target Character strings indicating to which target should the
#'   newly simulated be used to calculate f2. Only applicable when
#'   \code{sim.dp.out} is not missing.
#' @param time.unit Character strings indicating the unit of time. It should
#'   be either \code{"min"} for minute or \code{"h"} for hour. It is mainly
#'   used for checking CV rules and making plot. @seealso [calcf2()].
#' @param message Logic value. If \code{TRUE}, basic information of the
#'   simulation will be printed on screen.
#' @param plot A logic value indicating if the dissolution versus time plot
#'   should be printed.
#'
#' @return A list of 2 components: a data frame of model parameters and a data
#'   frame of mean dissolution generated using the said model parameters.
#'   The output can be passed to function \code{sim.dp()} to simulate multiple
#'   individual profiles.
#'
#' @details
#'
#'   The main principle of the function is as follows:
#'
#'   1. For any given mean dissolution profile, fit a suitable mathematical
#'      model and obtain model parameters.
#'       - No precise fitting is required since those parameters will be
#'         served as initial value for further model fitting.
#'       - If the output of the function \code{sim.dp()} is available, no
#'         initial fitting is necessary as model parameters can be read directly
#'         from the output, unless multivariate normal distribution approach
#'         is used. In such case, initial model fiting will be done.
#'   2. Find a suitable model parameters and simulate a new data set, comparing
#'      the new data to the provided reference data by calculating f2. If the
#'      f2 is equal to the \code{target.f2}, or within the lower and upper
#'      limit of the \code{target.f2} (see below), then output the obtained
#'      model parameters. There are two approaches used to find the suitable
#'      model parameters:
#'       - If \code{target.f2} is a single value, optimization algorithm
#'         will be used and the simulated dissolution profile will have f2
#'         equal to \code{target.f2} when compared to \code{dp} (within
#'         numeric precision defined by the tolerance).
#'       - If \code{target.f2} is vector of two numbers representing the
#'         lower and upper limit of target f2 value, such as \code{target.f2 =
#'         c(lower, upper)}, then dissolution will be obtained by random
#'         searching and the calculated f2 will be will within range of lower
#'         and upper. For example, set \code{target.f2 = c(54.95, 55.04)} to
#'         have target f2 of 55. Since f2 should be normally reported without
#'         decimal, in practice, the precision is enough. You might be able to
#'         do with \code{c(54.99, 55.01)} if you really need more precision.
#'         However, the narrower the range, the longer time it takes to run.
#'         With narrow range such as \code{c(54.999, 55.001)}, it is likely the
#'         program will hang up. In such case, provide single value to use
#'         optimization algorithm.
#'
#' Arguments \code{model.par.cv}, \code{fix.fmax.cv}, and \code{random.factor}
#' are certain numeric values used to better control the random generation of
#' model parameters. The default values should work in most scenarios. Those
#' values should be changed only when the function failed to return any value.
#' A vignette will be produced later to give more details.
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
#' model.par2 <- sim.dp.byf2(sim.dp.out = d.r, target.f2 = c(59.95, 60.04))
#'
#' @export
sim.dp.byf2 <- function(tp, dp, sim.dp.out, target.f2, seed,
                        regulation = c("EMA", "FDA", "WHO"),
                        model = c("Weibull", "first-order"), digits = 2L,
                        min.points = 3L, both.TR.85 = FALSE, max.disso = 105,
                        model.par.cv = 50, fix.fmax.cv = 0, random.factor = 3,
                        sim.target = c("ref.pop", "ref.median", "ref.mean"),
                        time.unit = c("min", "h"), message = TRUE, plot = TRUE){
  # initial check --------------------------------------------------------------
  regulation <- match.arg(regulation)
  model      <- match.arg(model)
  sim.target <- match.arg(sim.target)
  time.unit  <- match.arg(time.unit)

  # save state of the RNG
  # if (exists(".Random.seed", .GlobalEnv))
  #   oldseed <- .GlobalEnv$.Random.seed
  # else
  #   oldseed <- NULL

  # if no seed is provided, generate one
  if (missing(seed)) {
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
    if (all(!("sim.dp.out" %in% class(sim.dp.out)),
            any(missing(tp), missing(dp)))) {
      stop("The input data 'sim.dp.out' is incorrect. Use either the output\n",
           "of 'sim.dp' function or provide 'tp' and 'dp'.\n")
    }

    if (all(isTRUE(message), any(!missing(tp), !missing(dp)))) {
      cat("The output of the function 'sim.dp' was provided as input data,\n",
          "therefore, 'tp'/'dp' will be ignored.\n\n", sep = "")
    }

    # get ref data for calcf2() function
    tp         <- sim.dp.out$sim.summary$time[sim.dp.out$sim.summary$time > 0]
    ref.dp     <- sim.dp.out$sim.summary$dp[sim.dp.out$sim.summary$time > 0]
    ref.median <- sim.dp.out$sim.summary$sim.median[sim.dp.out$sim.summary$time > 0]
    ref.mean   <- sim.dp.out$sim.summary$sim.mean[sim.dp.out$sim.summary$time > 0]
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
    fmax.t <- sim.dp.out$sim.info$fmax
    tlag.t <- sim.dp.out$sim.info$tlag
    tadj.t <- tp - tlag.t
    tadj.t[tadj.t < 0] <- 0

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
                            k1        = sim.dp.out$sim.info$k1,
                            k1.cv     = sim.dp.out$sim.info$k1.cv,
                            time.unit = sim.dp.out$sim.info$time.unit,
                            stringsAsFactors = FALSE)
    } else {# model == NA, simulated by multivariate normal distribution
      # use regression to get mod.par in this case. 'model' from function input.
      mod.par <- mod.ref(tp = tp, ref.dp = data.r$dp, digits = digits,
                         model = model, time.unit = time.unit)
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
    tp     <- tp[tp > 0]
    ref.dp <- dp[tp > 0]
    data.r <- data.frame(tp = tp, dp = ref.dp, stringsAsFactors = FALSE)

    # function for modelling ref.pop to get parameters # add here mod.ref
    mod.par <- mod.ref(tp = tp, ref.dp = data.r$dp, digits = digits,
                       model = model, time.unit = time.unit)


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
          return(abs(tmp.f2[[1]] - target.f2))
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
                                f2 = f2.tmp[[1]], time.point = f2.tmp[[2]],
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
          return(abs(tmp.f2[[1]] - target.f2))
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
                                f2 = f2.tmp[[1]], time.point = f2.tmp[[2]],
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      }
    } else {# first-order optim
      mpar <- c(
        k1 = mod.par$k1*exp(rnorm(1, 0, mod.par$k1.cv*random.factor/100)),
        k1.cv = mod.par$k1.cv
      )

      opt.first.order <- function(mpar, target.f2, ...) {
        tmp.t  <- fmax.t*(1 - exp(-mpar[["k1"]]*tadj.t))
        data.t <- data.frame(tp = tp, dp = tmp.t, stringsAsFactors = FALSE)
        tmp.f2 <- calcf2(data.t, data.r, regulation = regulation,
                         cv.rule = FALSE, min.points = min.points,
                         both.TR.85 = both.TR.85, message = FALSE,
                         plot = FALSE, time.unit = time.unit)
        return(abs(tmp.f2[[1]] - target.f2))
      }

      res.f2 <- optim(mpar, opt.first.order, target.f2 = target.f2,
                      na.rm = TRUE)

      dp.2 <- fmax.t*(1 - exp(-res.f2$par[["k1"]]*tadj.t))
      dp.2[!is.finite(dp.2)] <- 0
      dp.2[dp.2 < 0] <- 0

      data.t <- data.frame(tp = tp, dp = dp.2, stringsAsFactors = FALSE)
      f2.tmp <- calcf2(data.t, data.r, regulation = regulation,
                       cv.rule = FALSE, min.points = min.points,
                       both.TR.85 = both.TR.85, message = FALSE,
                       plot = FALSE, time.unit = time.unit)

      model.par <- data.frame(model = "first-order", seed = seed,
                              fmax = fmax.t, tlag = tlag.t,
                              k1 = res.f2$par[["k1"]],
                              f2 = f2.tmp[[1]], time.point = f2.tmp[[2]],
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
        k1.2 <- mod.par$k1*exp(rnorm(1, 0, model.par.cv/100))
        dp.2 <- fmax.2*(1 - exp(-k1.2*tp.adj2))
      }

      dp.2[!is.finite(dp.2)] <- 0
      dp.2[dp.2 < 0] <- 0
      data.t <- data.frame(tp = tp, dp = dp.2, stringsAsFactors = FALSE)

      # get f2 and compare to criteria ----------------------------------
      f2.tmp <- calcf2(data.t, data.r, regulation = regulation,
                       cv.rule = FALSE, digits = digits, message = FALSE,
                       min.points = min.points, both.TR.85 = both.TR.85,
                       plot = FALSE, time.unit = time.unit)
      if (all(is.finite(f2.tmp[[1]]), f2.tmp[[1]] >= f2.lower,
              f2.tmp[[1]] <= f2.upper)) break
    }# end repeat

    if (mod.par$model == "Weibull") {
      if ("mdt" %in% names(mod.par)) {
        model.par <- data.frame(model = "Weibull", seed = seed, fmax = fmax.2,
                                tlag = tlag.2, mdt = mdt.2, beta = beta.2,
                                f2 = f2.tmp[[1]], time.point = f2.tmp[[2]],
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      } else {
        model.par <- data.frame(model = "Weibull", seed = seed, fmax = fmax.2,
                                tlag = tlag.2, alpha = alpha.2, beta = beta.2,
                                f2 = f2.tmp[[1]], time.point = f2.tmp[[2]],
                                regulation = regulation,
                                stringsAsFactors = FALSE)
      }
    } else {# first-order
      model.par <- data.frame(model = "first-order", seed = seed, fmax = fmax.2,
                              tlag = tlag.2, k1 = k1.2, f2 = f2.tmp[[1]],
                              time.point = f2.tmp[[2]], regulation = regulation,
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
                     round(f2.tmp[[1]], digits)),
              subtitle = paste0("Blue is the target profile (reference) and ",
                                "red is the simulated profile (test) with ",
                                "the predefined f2."))
    print(p1)
  }

  sim.out <- list(model.par = model.par, sim.disso = sim.disso)

  if (isTRUE(message)) {
    cat("Obtained model parameters and calculated f2 are:\n")
    print(sim.out$model.par)

    cat("\nAnd the difference between simulated test and reference is:\n")
    print(sim.out$sim.disso)
  }

  # if (!is.null(oldseed))
  #  .GlobalEnv$.Random.seed <- oldseed
  # else
  #   rm(".Random.seed", envir = .GlobalEnv)

  return(sim.out)
}
