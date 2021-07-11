#' Simulate Dissolution Profiles
#'
#' Function to simulate dissolution profiles based on mathematical models or
#' multivariate normal distribution.
#'
#' @importFrom MASS mvrnorm
#' @usage
#' sim.dp(tp, model = c("Weibull", "first-order"), model.par,
#'        seed, product, dp, dp.cv, ascending = FALSE, n.units = 12L,
#'        max.disso = 105, message = FALSE, plot = TRUE,
#'        time.unit = c("min", "h"), plot.max.unit = 36L)
#'
#' @param tp Numeric vector of time points for the dissolution profiles.
#' @param model Character strings of model names. Currently only 'Weibull' and
#'   'first-order' models are supported.
#' @param model.par A list with model parameters. If missing, a list of random
#'   \code{model.par} will be generated. See Details.
#' @param seed Numeric seed value for reproducibility. If missing, a random
#'   seed will be generated for reproducibility purpose.
#' @param product Character strings indicating the product name of the
#'   simulated profiles. If missing, a random name with 3 letters and 4 digits
#'   will be generated.
#' @param dp,dp.cv Numeric vectors of mean dissolution profile (\code{dp}) and
#'   its corresponding CV% (\code{dp.cv}). See Details.
#' @param ascending Logical. If \code{TRUE}, simulated profiles will always
#'   increase with time. Only applicable when See Details.
#' @param n.units An integer indicating the number of individual profiles
#'   to be simulated.
#' @param max.disso Numeric value for the maximum possible dissolution.
#'   In theory, the maximum dissolution should be 100%, but in practice,
#'   it is not uncommon to see higher values, such as 105%.
#' @param message Logical. If \code{TRUE}, basic information of the
#'   simulation will be printed on screen.
#' @param plot Logical If \code{TRUE}, a plot will be included in the output.
#' @param time.unit Character strings indicating the unit of time.
#'   It should be either \code{"min"} for minute or \code{"h"} for hour.
#'   It is mainly used for generating \code{model.par} and \code{dp.cv}
#'   when they are missing, and making plot. @seealso [calcf2()].
#' @param plot.max.unit Integer. If the number of individual units is no more
#'   than this value, mean and all individual profiles will be plotted;
#'   otherwise, individual profiles will represented by boxplots.
#'   Therefore, to avoid overplotting, this value should not be too big.
#'
#' @return A list of 3 to 5 components:
#'   - \code{sim.summary}: A data frame with summary statistics of all
#'     individual profiles.
#'   - \code{sim.disso}: A data frame with all individual profiles.
#'   - \code{sim.info}: A data frame with information of simulations such as
#'     the seed number. If modelling approach is used, the data frame contains
#'     model parameters; if multivariate normal distribution approach is used,
#'     the data frame contains \code{dp/dp.cv}.
#'   - \code{model.par.ind}: A data frame of all individual model parameters
#'     that were used for the simulation of individual dissolution profiles.
#'     Available only if the modelling approach is used, i.e., when \code{dp}
#'     is not missing.
#'   - \code{sim.dp}: A plot. Available only if \code{plot = TRUE}.
#'
#' @details
#'
#' The approach used to simulate individual dissolution profiles depends on
#' if the mean dissolution profile \code{dp} is provided by user or not.
#' If it is not provided, \code{dp} and individual profiles will be simulated
#' using mathematical models, otherwise it will be simulated using multivariate
#' normal distribution. For the latter, the CV at time points \code{dp.cv}
#' should be supplier by user; if not, it will be generated automatically.
#'
#' ## Use models (Recommended Approach)
#'
#'   The first-order model is expressed as \deqn{f_t = f_\mathrm{max}%
#'   \left(1 - e^{-k_1\(t - t_\mathrm{lag}\right)}\right).}{f(t) = fmax%
#'   (1 - exp(-k1(t - tlag))),}
#'   and the Weibull model was expressed either as \deqn{f_t = f_\mathrm{max}%
#'   \left(1 - e^{-\left(\frac{t - t_\mathrm{lag}}{\mathrm{MDT}}%
#'   \right)^\beta}\right)}{f(t) = fmax (1 - exp(-((t - tlag)/MDT)^\beta))} or
#'   \deqn{f_t = f_\mathrm{max}\left(1 - e^{-\frac{(t - t_\mathrm{lag})^\beta}%
#'   {\alpha}}\right)}{f(t) = fmax (1 - exp(-(((t - tlag)^\beta)/\alpha))),}
#'   where \eqn{f_\mathrm{max}}{fmax} is the maximum dissolution, MDT is
#'   the mean dissolution time, \eqn{t_\mathrm{lag}}{tlag} is the lag time,
#'   \eqn{\alpha}{\alpha} and \eqn{\beta}{\beta} are the scale and shape factor
#'   in Weibull function, and \eqn{k_1}{k1} is the rate constant in the
#'   first-order model. Obviously, The two Weibull models are mathematically
#'   equivalent by letting \eqn{\alpha = \mathrm{MDT}^\beta}{\alpha = MDT^\beta}.
#'   The second expression of Weibull model was included to be compatible to
#'   \code{DDSolver} program.
#'
#'   Individual model parameter were simulated according to the exponential
#'   error model \deqn{P_i = P_\mu e^{N\left(0, \sigma^2\right)}}{P(i) = P(\mu)%
#'   exp(N(0, \sigma^2)),}
#'   where \eqn{N\left(0, \sigma^2\right)}{N(0, \sigma^2)} represents the
#'   normal distribution with mean zero and variance \eqn{\sigma^2} (σ = CV/100);
#'   \eqn{P_i}{P(i)} and \eqn{P_\mu}{P(\mu)} denote the individual and
#'   population model parameters.
#'
#'   Therefore, \code{model.par} should be supplied as a named list of 6
#'   parameters for the first-order model (\code{fmax/fmax.cv}, \code{k1/k1.cv},
#'   and \code{tlag/tlag.cv}), and 8 parameters for Weibull model
#'   (\code{fmax/fmax.cv}, \code{tlag/tlag.cv}, \code{beta/beta.cv}, and
#'   either \code{alpha/alpha.cv} or \code{mdt/mdt.cv}, depending on the
#'   mathematical formula used). For example:
#'
#'   - \code{model.par = list(fmax = 100, fmax.cv = 5, k1 = 0.6, k1.cv = 25,
#'     tlag = 0, tlag.cv = 0)} for the first-order model.
#'   - \code{model.par = list(fmax = 100, fmax.cv = 5, tlag = 5, tlag.cv = 10,
#'     mdt = 15, mdt.cv = 20, beta = 1.5, beta.cv = 5)}, or
#'   - \code{model.par = list(fmax = 100, fmax.cv = 5, tlag = 5, tlag.cv = 10,
#'     alpha = 60, alpha.cv = 20, beta = 1.5, beta.cv = 5)} for the Weibull
#'     model.
#'
#'   See vignette for more details.
#'
#'   ## Use multivariate normal distribution
#'
#'   If \code{dp} is supplied by user, simulation will be done based on
#'   multivariate normal distribution, and \code{n.units} of individual
#'   profiles will be simulated with mean dissolution profile equals to
#'   \code{dp}. IF CV at time points (\code{dp.cv}) is supplier, the CV% of
#'   simulated profiles will equal to \code{dp.cv}. If \code{dp.cv} is missing,
#'   it will be simulated such that the CV is 20% for time points up to 10 min,
#'   5% for time points where dissolution is greater than 90% (more realistic),
#'   and 10% for the rest time points. The length of \code{dp} and \code{dp.cv}
#'   should be the same.
#'
#'   When this method is used, depending on \code{dp/dp.cv}, there is no
#'   guarantee that all individual profiles increase with time; near the end
#'   of the time points, some profiles may decrease. This is realistic for
#'   those drug substances that are unstable in the dissolution medium.
#'   To avoid this, set \code{ascending = TRUE}.
#'
#'
#' @examples
#' # time points
#' tp <- c(5, 10, 15, 20, 30, 45, 60)
#' # using all default values to simulate profiles
#' d1 <- sim.dp(tp, plot = FALSE)
#'
#' # summary statistics
#' d1$sim.summary
#'
#' # individual profiles
#' d1$sim.disso
#'
#' # simulation information
#' d1$sim.info
#'
#' #individual model parameters
#' d1$mod.par.ind
#'
#' # plot with individual and mean profiles
#' # d1$sim.plot
#'
#' # using Weibull model to simulate 24 profiles without lag time
#' mod.par <- list(fmax = 100, fmax.cv = 2, tlag = 0, tlag.cv = 0,
#'                 mdt = 20, mdt.cv = 5, beta = 2.2, beta.cv = 5)
#' d2 <- sim.dp(tp, model.par = mod.par, seed = 100, n.units = 100L,
#'              plot = FALSE)
#'
#' # mean plot and boxplot
#' # d2$sim.plot
#'
#' # using multivariate normal distribution approach
#' # mean profile same length as tp
#' dp <- c(5, 20, 42, 63, 88, 91, 99)
#'
#' # CV% at each time point
#' dp.cv <- c(20, 20, 10, 10, 10, 10, 5)
#'
#' d3 <- sim.dp(tp, dp = dp, dp.cv = dp.cv, seed = 100, plot = FALSE)
#'
#' @export
sim.dp <- function(tp, model = c("Weibull", "first-order"), model.par,
                   seed, product, dp, dp.cv, ascending = FALSE, n.units = 12L,
                   max.disso = 105, message = FALSE, plot = TRUE,
                   time.unit = c("min", "h"), plot.max.unit = 36L){
  # initial checking -----------------------------------------------------------
  model     <- match.arg(model)
  time.unit <- match.arg(time.unit)

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

  # tp is mandatory
  if (missing(tp)) {
    stop("Time points ('tp') has to be provided.\n")
  }

  # if model.par is given by user, check for error
  if (!missing(model.par)){
    if (model == "first-order") {
      if (length(setdiff(c("fmax", "fmax.cv", "k1", "k1.cv", "tlag", "tlag.cv"),
                         names(model.par))) != 0) {
        stop("Model parameters are incorrect. Three pairs of parameters are\n",
             "required: 'fmax/fmax.cv', 'k1/k1.cv', and 'tlag/tlag.cv'.\n")
      }
    } else {
      if (all(length(setdiff(c("fmax", "fmax.cv", "tlag", "tlag.cv",
                               "alpha", "alpha.cv", "beta", "beta.cv"),
                             names(model.par))) != 0,
              length(setdiff(c("fmax", "fmax.cv", "tlag", "tlag.cv",
                               "mdt", "mdt.cv", "beta", "beta.cv"),
                             names(model.par))) != 0)) {
        stop("Model parameters are incorrect. Four pairs of parameters are\n",
             "required: 'fmax/fmax.cv', 'tlag/tlag.cv', 'beta/beta.cv', and,\n",
             "depending on the mathematical expression of the model, either\n",
             "'alpha/alpha.cv' or 'mdt/mdt.cv'. Check your data.\n")
      }
    }
  }# end checking model.par

  # check length, tp = dp = dp.cv
  if (!missing(dp)) {
    if (length(tp) != length(dp)) {
      stop("Length of 'tp' and 'dp' should be equal. Check your data.")
    }
    if (all(!missing(dp.cv), length(dp) != length(dp.cv))) {
      stop("Length of 'dp' and 'dp.cv' should be equal. Check your data.")
    }
  }

  # if dp.cv is provided without dp
  if (all(missing(dp), !missing(dp.cv))) {
    stop("Variation of mean dissolution profile at each time point ('dp.cv')",
         "\nwas given but the mean profile was missing, please provide it.")
  }

  # population approach when dp is missing -------------------------------------
  if (missing(dp)) {
    use.model <- TRUE

    # remove time 0 if any
    tp <- tp[tp > 0]

    ## generate default model parameters if missing ------------------
    if (missing(model.par)) {
      if (model == "first-order") {
        if (time.unit == "min") {# typical IR dosage form if unit is 'min'
          model.par <- list(fmax = round(100*exp(rnorm(1, 0, 0.025)), 2),
                            fmax.cv = 3,
                            k1 = round(0.15*exp(rnorm(1, 0, 0.4)), 6),
                            k1.cv = 40,
                            tlag = 0, tlag.cv = 0)
        } else {# typical ER dosage form if unit is 'h'
          model.par <- list(fmax = round(98*exp(rnorm(1, 0, 0.03)), 2),
                            fmax.cv = 3,
                            k1 = round(0.3*exp(rnorm(1, 0, 0.4)), 6),
                            k1.cv = 40,
                            tlag = round(0.3*exp(rnorm(1, 0, 0.25)), 2),
                            tlag.cv = 20)
        }# end missing model.par for first-order model
      } else {# for Weibull model
        if (time.unit == "min") {# typical IR dosage form if unit is 'min'
          model.par <- list(fmax = round(100*exp(rnorm(1, 0, 0.025)), 2),
                            fmax.cv = 3,
                            tlag = 0, tlag.cv = 0,
                            mdt = round(8*exp(rnorm(1, 0, 0.3)), 2),
                            mdt.cv = 30,
                            beta = round(0.8*exp(rnorm(1, 0, 0.2)), 6),
                            beta.cv = 20)
        } else {# typical ER dosage form if unit is 'h'
          model.par <- list(fmax = round(98*exp(rnorm(1, 0, 0.05)), 2),
                            fmax.cv = 5,
                            tlag = round(0.5*exp(rnorm(1, 0, 0.3)), 2),
                            tlag.cv = 30,
                            mdt = round(5*exp(rnorm(1, 0, 0.2)), 2),
                            mdt.cv = 20,
                            beta = round(2*exp(rnorm(1, 0, 0.3)), 6),
                            beta.cv = 30)
        }# end ER dosage form
      }# end missing model.par for Weibull model
    }# end missing(model.par)

    ## dissolution profile simulation --------------------------------
    # prepare common individual parameters for both models
    fmax.ind <- model.par$fmax*exp(rnorm(n.units, 0, model.par$fmax.cv/100))
    fmax.ind[fmax.ind > max.disso] <- max.disso
    tlag.ind <- model.par$tlag*exp(rnorm(n.units, 0, model.par$tlag.cv/100))

    # for mean dissolution profile
    tp.adj <- tp - model.par$tlag
    tp.adj[tp.adj < 0] <- 0

    if (model == "first-order") {### first-order -----------
      # target mean dissolution profile
      dp <- model.par$fmax*(1 - exp(-model.par$k1*tp.adj))

      # individual model par
      k1.ind <- model.par$k1*exp(rnorm(n.units, 0, model.par$k1.cv/100))

      # individual dissolution data. each column is one unit
      dp.ind <- matrix(NA, nrow = length(tp), ncol = n.units)

      for (i in seq_len(n.units)) {
        tp.adj.ind <- tp - tlag.ind[i]
        tp.adj.ind[tp.adj.ind < 0] <- 0
        dp.ind[, i] <- fmax.ind[i]*(1 - exp(-k1.ind[i]*tp.adj.ind))
      }

      # save individual parameters for output
      model.par.ind <- data.frame(
        fmax.ind = fmax.ind,
        k1.ind   = k1.ind,
        tlag.ind = tlag.ind,
        tringsAsFactors = FALSE
      )

      # output simulation info
      sim.info <- data.frame(
        fmax      = model.par$fmax,
        fmax.cv   = model.par$fmax.cv,
        k1        = model.par$k1,
        k1.cv     = model.par$k1.cv,
        tlag      = model.par$tlag,
        tlag.cv   = model.par$tlag.cv,
        seed      = seed,
        n.units   = n.units,
        max.disso = max.disso,
        model     = "first-order",
        time.unit = time.unit,
        stringsAsFactors = FALSE
      )
    } else {### Weibull model ------------------------------
      # individual parameters
      beta.ind <- model.par$beta*exp(rnorm(n.units, 0, model.par$beta.cv/100))

      dp.ind <- matrix(NA, nrow = length(tp), ncol = n.units)

      if ("mdt" %in% names(model.par)) {
        # target mean dissolution profile
        dp <- model.par$fmax*(1 - exp(-(tp.adj/model.par$mdt)^model.par$beta))

        mdt.ind <- model.par$mdt*exp(rnorm(n.units, 0, model.par$mdt.cv/100))

        for (i in seq_len(n.units)) {
          tp.adj.ind <- tp - tlag.ind[i]
          tp.adj.ind[tp.adj.ind < 0] <- 0
          dp.ind[, i] <-
            fmax.ind[i]*(1 - exp(-(tp.adj.ind/mdt.ind[i])^beta.ind[i]))
        }

        # save individual parameters for output
        model.par.ind <- data.frame(
          fmax.ind = fmax.ind,
          mdt.ind  = mdt.ind,
          beta.ind = beta.ind,
          tlag.ind = tlag.ind,
          stringsAsFactors = FALSE
        )

        # output simulation info
        sim.info <- data.frame(
          fmax      = model.par$fmax,
          fmax.cv   = model.par$fmax.cv,
          mdt       = model.par$mdt,
          mdt.cv    = model.par$mdt.cv,
          beta      = model.par$beta,
          beta.cv   = model.par$beta.cv,
          tlag      = model.par$tlag,
          tlag.cv   = model.par$tlag.cv,
          seed      = seed,
          n.units   = n.units,
          max.disso = max.disso,
          model     = "Weibull",
          time.unit = time.unit,
          stringsAsFactors = FALSE
        )
      } else {# compatible with model expression in DDSolver
        # target mean dissolution profile
        dp <- model.par$fmax*(1-exp(-((tp.adj^model.par$beta)/model.par$alpha)))

        alpha.ind <- model.par$alpha*exp(rnorm(n.units, 0, model.par$alpha.cv/100))

        for (i in seq_len(n.units)) {
          tp.adj.ind <- tp - tlag.ind[i]
          tp.adj.ind[tp.adj.ind < 0] <- 0
          dp.ind[, i] <-
            fmax.ind[i]*(1 - exp(-((tp.adj.ind^beta.ind[i])/alpha.ind[i])))
        }

        # save individual parameters for output
        model.par.ind <- data.frame(
          fmax.ind  = fmax.ind,
          alpha.ind = alpha.ind,
          beta.ind  = beta.ind,
          tlag.ind  = tlag.ind,
          stringsAsFactors = FALSE
        )

        # output info
        sim.info <- data.frame(
          fmax      = model.par$fmax,
          fmax.cv   = model.par$fmax.cv,
          alpha     = model.par$alpha,
          alpha.cv  = model.par$alpha.cv,
          beta      = model.par$beta,
          beta.cv   = model.par$beta.cv,
          tlag      = model.par$tlag,
          tlag.cv   = model.par$tlag.cv,
          seed      = seed,
          n.units   = n.units,
          max.disso = max.disso,
          model     = "Weibull",
          time.unit = time.unit,
          stringsAsFactors = FALSE
        )
      }# end Weibull expression with alpha
    }# end Weibull model

    ## output message ------------------------------------------------
    if (isTRUE(message)) {
      cat("Dissolution data was generated using ", model, " model, ",
          "with the \nfollowing model parameters:\n", sep = "")
      cat("- fmax      = ", model.par$fmax, "\n", sep = "")
      cat("- fmax.cv   = ", model.par$fmax.cv, "%\n", sep = "")
      cat("- tlag      = ", model.par$tlag, "\n", sep = "")
      cat("- tlag.cv   = ", model.par$tlag.cv, "%\n", sep = "")
      if (model == "first-order") {
        cat("- k1        = ", model.par$k1, "\n", sep = "")
        cat("- k1.cv     = ", model.par$k1.cv, "%\n\n", sep = "")
      } else if (model == "Weibull") {
        if ("mdt" %in% names(model.par)) {
          cat("- mdt       = ", model.par$mdt, "\n", sep = "")
          cat("- mdt.cv    = ", model.par$mdt.cv, "%\n", sep = "")
        } else {
          cat("- alpha     = ", model.par$alpha, "\n", sep = "")
          cat("- alpha.cv  = ", model.par$alpha.cv, "%\n", sep = "")
        }
        cat("- beta      = ", model.par$beta,  "\n", sep = "")
        cat("- beta.cv   = ", model.par$beta.cv,  "%\n\n", sep = "")
      }
      cat("\nSeed number used: ", seed, "\n")
    } # end message
  # end missing dp
  } else {
    # with dp, using multivariate normal distribution approach ---------
    use.model <- FALSE

    # remove time 0 from tp if any, and corresponding point from dp,
    # and if dp.cv is not missing, from dp.cv as well
    if (isTRUE((length(which(tp == 0)) != 0))) {
      tp <- tp[-(which(tp == 0))]
      dp <- dp[-(which(tp == 0))]
      if (!missing(dp.cv)) dp.cv[-(which(tp == 0))]
    }

    ## generate dp.cv if missing -------------------------------------
    # CV 20% for time points up to 10 min, 10% for the rest points,
    # unless release is > 90%, in which case CV 5%.
    if (missing(dp.cv)) {
      # initialize cv 10% for all points.
      dp.cv <- rep(10, length(tp))

      # check if there is any time point <= 10, if so, cv = 20%
      if (time.unit == "min") {
        if (length(which(tp <= 10)) != 0) {
          dp.cv[which(tp <= 10)] <- 20
        }
      } else {# convert h to min
        if (length(which(tp*60 <= 10)) != 0) {
          dp.cv[which(tp*60 <= 10)] <- 20
        }
      }# end checking time points <= 10 min

      # if dp > 90%, then cv about 5%, which is more realistic
      dp.cv <- ifelse(dp >= 90, 5, dp.cv)

      # lastly, check if there's dp = 0, if so, cv = 0.
      dp.cv[dp <= 0] <- 0

      if (isTRUE(message)) {
        cat("Multivariate normal distribution was chosen for the method of\n",
            "simulation but variation of mean dissolution profile at each\n",
            "time point ('dp.cv') was missing; therefore, it was generated\n",
            "such that the CV is 20% for time points up to 10 min, 5% for\n",
            "time points where dissolution is greater than 90%, and 10% for\n",
            "the rest time points.\n\n", sep = "")
      }
    } else {# if dp.cv not missing
      dp.cv[!is.finite(dp.cv)] <- 0
    }

    # Get sd for random number generation
    sd <- dp.cv/100*dp
    sd[!is.finite(sd)] <- 0
    covar <- sd%o%sd

    ## generate individual data dp1 based on the mean profile dp -----
    if (length(dp) <= n.units) {
      repeat {# each row is a unit. use MASS::mvrnorm
        dp1 <- mvrnorm(n = n.units, mu = dp, Sigma = covar, empirical = TRUE)

        if (isFALSE(ascending)) {
          if (max(dp1) < max.disso) break
        } else {# could be very slow! Need better solution!
          tmp0 <- matrix(rep(1:length(tp), n.units), nrow = n.units,
                         ncol = length(tp), byrow = TRUE)
          tmp1 <- matrix(NA, nrow = n.units, ncol = length(tp))
          for (i in 1:NROW(dp1)) {
            tmp1[i, ] <- order(dp1[i, ])
          }
          if (identical(tmp1, tmp0) && max(dp1) < max.disso) break
        }
      }# end repeat
    } else {# not common but still need better solution in such cases.
      stop("To simulate dissolution profiles based on multivariate normal\n",
           "distribution, the number of dissolution units ('n.units') should\n",
           "not be less than number of time points. Please either decrease\n",
           "the number of time points or increase 'n.units'. Alternatively,\n",
           "you can use population modelling approach.\n\n")
    }

    # change to each column for each unit
    dp.ind <- t(dp1)

    # output simulation info
    sim.info <- data.frame(seed = seed, n.units = n.units,
                           max.disso = max.disso, model = NA,
                           time.unit = time.unit,
                           stringsAsFactors = FALSE)
  }# end multivariate normal distribution approach

  # basic descriptive statistics -----------------------------------------------
  # add names of formulation/batch if missing. 3 letters + 4 numbers
  if (missing(product)) {
    product <-
      paste0(paste0(sample(LETTERS, 3, replace = TRUE), collapse = ""),
             sample(seq(1000, 9999, 1), 1), collapse = "")
  }

  # add time 0 for better plot
  tp    <- c(0, tp)
  dp    <- c(0, dp)
  dp.cv <- ifelse(isTRUE(use.model), rep(NA, length(tp)), c(0, dp.cv))

  dp.ind   <- rbind(rep(0, NCOL(dp.ind)), dp.ind)

  sim.mean <- rowMeans(dp.ind, na.rm = TRUE)
  sim.mean[!is.finite(sim.mean)] <- 0

  sim.var <- apply(dp.ind, 1, var)
  sim.var[!is.finite(sim.var)] <- 0

  sim.sd  <- sqrt(sim.var)
  sim.sd[!is.finite(sim.sd)] <- 0

  sim.cv  <- sim.sd/sim.mean*100
  sim.cv[!is.finite(sim.cv)] <- 0

  dp.summary <- data.frame(
    product    = rep(product, length(tp)),
    time       = tp,
    dp         = dp,
    dp.cv      = dp.cv,
    sim.mean   = sim.mean,
    sim.median = apply(dp.ind, 1, median, na.rm = TRUE),
    sim.cv     = sim.cv,
    sim.var    = sim.var,
    sim.sd     = sim.sd,
    sim.min    = apply(dp.ind, 1, min, na.rm = TRUE),
    sim.max    = apply(dp.ind, 1, max, na.rm = TRUE),
    sim.qt05   = apply(dp.ind, 1, quantile, probs = 0.05, na.rm = TRUE),
    sim.qt25   = apply(dp.ind, 1, quantile, probs = 0.25, na.rm = TRUE),
    sim.qt75   = apply(dp.ind, 1, quantile, probs = 0.75, na.rm = TRUE),
    sim.qt95   = apply(dp.ind, 1, quantile, probs = 0.95, na.rm = TRUE),
    sim.bpmin  = apply(dp.ind, 1, bpwhisker.l),#16,17 column for plot only
    sim.bpmax  = apply(dp.ind, 1, bpwhisker.u),
    stringsAsFactors = FALSE
  )

  dp.ind.out <- data.frame(time = tp, dp.ind, stringsAsFactors = FALSE)

  # generate names for each unit
  unit.name <- c(rep(NA, n.units))
  index <- seq_len(n.units)
  z0 <- nchar(format(n.units, scientific = FALSE))

  for (i in seq_len(z0)) {
    j <- index[index < 10^i & index >= 10^(i - 1)]
    unit.name[j] <- paste0(
      paste0("unit.", paste0(rep(0, z0 - i), collapse = "")), j)
  }

  names(dp.ind.out)[-1] <- unit.name

  # making plot ----------------------------------------------------------------
  if (isTRUE(plot)) {
    # need this to get rid of "no visible binding for global variable" notes:
    time <- release <- sim.bpmin <- sim.qt25 <- sim.qt75 <- sim.median <-
      sim.bpmax <- NULL

    # RColorBrewer::brewer.pal(11, "Paired") 1 "#A6CEE3" and 2 "#1F78B4"
    # for individual and mean profiles, 4 "#33A02C" for target population mean

    if (n.units <= plot.max.unit) {# print individual profiles
      p.dat <- reshape(dp.ind.out, varying = 2:(n.units + 1),
                       direction = "long", timevar = "unit",
                       v.names = "release", idvar = "time")

      y.tick.max <- ceiling(max(p.dat$release))

      p1 <- ggplot(p.dat, aes(x = time)) +
        geom_hline(yintercept = 85, size = 1.2, color = "gray",
                   alpha = 0.7, linetype = "dashed") +
        geom_line(aes(y = release, group = unit), size = 0.3,
                  color = "#A6CEE3") +
        geom_point(data = dp.summary, aes(y = sim.mean), size = 1.5,
                   color = "#1F78B4") +
        geom_line(data = dp.summary, aes(y = sim.mean), size = 1.1,
                  color = "#1F78B4") +
        geom_line(data = dp.summary, aes(y = dp), size = 1.1, alpha = 0.7,
                  color = "#33A02C") +
        theme_bw() +
        scale_y_continuous(limits = c(0, y.tick.max + 5),
                           breaks = seq(0, y.tick.max + 5, 5)) +
        scale_x_continuous(breaks = dp.summary$time) +
        ylab("Dissolution (%)") +
        xlab(paste0("Time (", time.unit, ")")) +
        ggtitle(paste0("Summary of ", n.units, " Simulated Dissolution ",
                       "Profiles for ", product, "."),
                subtitle = paste0("Blue: mean simulated profile; ",
                                  "Green: targeted profile; ",
                                  "Thin lines: individual profiles."))
    } else {#use boxplots
      y.tick.max <- ceiling(max(c(dp.summary$sim.qt95, dp.summary$sim.bpmax)))

      p1 <- ggplot(dp.summary, aes(x = time)) +
        geom_hline(yintercept = 85, size = 1.2, color = "gray",
                   alpha = 0.7, linetype = "dashed") +
        geom_boxplot(aes(ymin = sim.bpmin, lower = sim.qt25,
                         middle = sim.median, upper = sim.qt75,
                         ymax = sim.bpmax, group = time),
                     fill = "#A6CEE3", color = "#A6CEE3",
                     stat = "identity", width = 0.5, alpha = 0.5) +
        geom_point(aes(y = sim.mean), size = 1.5, color = "#1F78B4") +
        geom_line(aes(y = sim.mean), size = 1.1, color = "#1F78B4") +
        geom_line(aes(y = dp), size = 1.1, alpha = 0.7, color = "#33A02C") +
        theme_bw() +
        scale_y_continuous(limits = c(0, y.tick.max + 5),
                           breaks = seq(0, y.tick.max + 5, 5)) +
        scale_x_continuous(breaks = dp.summary$time) +
        ylab("Dissolution (%)") +
        xlab(paste0("Time (", time.unit, ")")) +
        ggtitle(paste0("Summary of ", n.units, " Simulated Dissolution ",
                       "Profiles for ", product, "."),
                subtitle = paste0("Blue: mean simulated profile; ",
                                  "Green: targeted profile; ",
                                  "Boxplot: summary of individual profiles."))
    }
  }# end plot

  # output ---------------------------------------------------------------------
  if (isTRUE(use.model)) {
    if (isTRUE(plot)) {
      sim.out <- list(sim.summary = dp.summary[, -c(16, 17)],
                      sim.disso = dp.ind.out, sim.info = sim.info,
                      model.par.ind = model.par.ind, sim.plot = p1)
    } else {
      sim.out <- list(sim.summary = dp.summary[, -c(16, 17)],
                      sim.disso = dp.ind.out, sim.info = sim.info,
                      model.par.ind = model.par.ind)
    }
  } else {# no model.par.ind if not using model
    if (isTRUE(plot)) {
      sim.out <- list(sim.summary = dp.summary[, -c(16, 17)],
                      sim.disso = dp.ind.out, sim.info = sim.info,
                      sim.plot = p1)
    } else {
      sim.out <- list(sim.summary = dp.summary[, -c(16, 17)],
                      sim.disso = dp.ind.out, sim.info = sim.info)
    }
  }

  # if (!is.null(oldseed))
  #  .GlobalEnv$.Random.seed <- oldseed
  # else
  #   rm(".Random.seed", envir = .GlobalEnv)

  # for another function find.dp.by.f2
  class(sim.out) <- append(class(sim.out), "sim.dp.out")
  return(sim.out)
}# end function sim.dp
