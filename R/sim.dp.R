#' Simulate Dissolution Profiles
#'
#' Function to simulate dissolution profiles based on mathematical models or
#' multivariate normal distribution.
#'
#' @importFrom MASS mvrnorm
#' @usage
#' sim.dp(tp, dp, dp.cv, model = c("Weibull", "first-order"),
#'        model.par = NULL, seed = NULL, n.units = 12L, product,
#'        max.disso = 100, ascending = FALSE, message = FALSE,
#'        time.unit = c("min", "h"), plot = TRUE,
#'        plot.max.unit = 36L)
#'
#' @param tp *Numeric* vector of time points for the dissolution profiles.
#'   See Details.
#' @param dp,dp.cv *Numeric* vectors of the *target mean dissolution profile*
#'   (`dp`) and its respective CV at time points `tp` (`dp.cv`). See Details.
#' @param model *Character* strings of model names. Currently only `"Weibull"`
#'   and `"first-order"` models are supported.
#' @param model.par A *list* with model parameters. If missing, a list of
#'   random `model.par` will be generated. See Details.
#' @param seed *Integer* seed value for reproducibility. If missing, a random
#'   seed will be generated for reproducibility purpose.
#' @param n.units An *integer* indicating the number of individual profiles
#'   to be simulated.
#' @param product *Character* strings indicating the product name of the
#'   simulated profiles. If missing, a random name with 3 letters and 4 digits
#'   will be generated.
#' @param max.disso *Numeric* value for the maximum possible dissolution.
#'   In theory, the maximum dissolution is 100%, but in practice, it is not
#'   uncommon to see higher values, such as 105%, or much lower values,
#'   especially for products with low solubility.
#' @param ascending *Logical*. If `TRUE`, simulated profiles will always
#'   increase with time. Only applicable when the approach based on
#'   multivariate normal distribution is used. See Details.
#' @param message *Logical*. If `TRUE`, basic information of the simulation
#'   will be printed on screen.
#' @param time.unit *Character* strings indicating the unit of time. It should
#'   be either `"min"` for minute or `"h"` for hour. It is mainly used for
#'   and making plot and generating `model.par` and `dp.cv` when they are not
#'   provided by the user. @seealso [calcf2()].
#' @param plot *Logical*. If `TRUE`, a a dissolution versus time plot will be
#'   included in the output.
#' @param plot.max.unit *Integer*. If the number of individual units is no more
#'   than this value, the mean and all individual profiles will be plotted;
#'   otherwise, individual profiles will be represented by boxplots at each
#'   time point. Therefore, to avoid overplotting, this value should not be
#'   too large. @seealso [calcf2()].
#'
#' @return A list of 3 to 5 components:
#'   - `sim.summary`: A *data frame* with summary statistics of all
#'     individual dissolution profiles.
#'   - `sim.disso`: A *data frame* with all individual dissolution profiles.
#'   - `sim.info`: A *data frame* with information of the simulation such as
#'     the seed number and number of individual profiles. If modelling approach
#'     is used, the data frame will contain model parameters as well.
#'   - `model.par.ind`: A *data frame* of all individual model parameters
#'     that were used for the simulation of individual dissolution profiles.
#'     Available only if the modelling approach is used, i.e., when `dp`
#'     is missing.
#'   - `sim.dp`: A plot. Available only if the argument `plot` is `TRUE`.
#'
#' @details
#' ## Simulation approaches
#' The approach used to simulate individual dissolution profiles depends on if
#' the *target mean dissolution profile* `dp` is provided by the user or not.
#' - If `dp` is not provided, then it will be calculated using `tp`, `model`,
#'   and `model.par`. The parameters defined by `model.par` are considered as
#'   the *population parameter*; consequently, the calculated `dp` is
#'   considered as the *targeted population profile*. In addition, `n.units`
#'   sets of *individual model parameters* will be simulated based on
#'   exponential error model, and *individual dissolution profiles* will be
#'   generated using those individual parameters. The mean of all simulated
#'   individual profiles, `sim.mean`, included in one of the outputs data
#'   frames, `sim.summary`, will be *similar, but not identical, to `dp`*.
#'   The difference between `sim.mean` and `dp` is determined by the number of
#'   units and the variability of the model parameters. In general, the larger
#'   the number of units and the lower of the variability, the smaller the
#'   difference. Additional details on the mathematical models are given below.
#' - If `dp` is provided, then `n.units` of individual dissolution profiles
#'   will be simulated using multivariate normal distribution. The mean of all
#'   simulated individual profiles, `sim.mean`, will be *identical to `dp`*.
#'   In such case, it is recommended that `dp.cv`, the CV at time points `tp`,
#'   should also be provided by the user. If `dp.cv` is not provided, it will
#'   be generated automatically such that the CV is between 10% and 20% for time
#'   points up to 10 min, between 1% and 3% for time points where the
#'   dissolution is more than 95%, between 3% and 5% for time points where the
#'   dissolution is between 90% and 95%, and between 5% and 10% for the rest of
#'   time points. Whether the `dp.cv` is provided or generated automatically,
#'   the CV calculated using all individual profiles will be equal to `dp.cv`.
#'   Additional details on this approach are given below.
#'
#' ## Minimum required arguments that must be provided by the user
#' If `dp` is provided by the user, logically, `tp` must also be provided, and
#' the approach based on multivariate normal distribution is used, as explained
#' above. If `dp` is not provided, `tp` could be missing, i.e., a simple
#' function call such as `sim.dp()` will simulate dissolution profiles. In such
#' case, a default `tp` will be generated depending on the `time.unit`: if the
#' `time.unit` is `"min"`, then `tp` would be `c(5, 10, 15, 20, 30, 45, 60)`;
#' otherwise the `tp` would be `c(1, 2, 3, 4, 6, 8, 10, 12)`. The rest
#' arguments are either optional, or required by the function but default
#' values will be used.
#'
#' ## Notes on mathematical models
#' The first-order model is expressed as
#' \deqn{f_t = f_\mathrm{max} \left(1 - %
#'   e^{-k\left(t - t_\mathrm{lag}\right)}\right),}{f(t) = fmax%
#'   (1 - exp(-k(t - tlag))),}
#' and the Weibull model was expressed either as
#' \deqn{f_t = f_\mathrm{max} \left(1 - %
#'   e^{-\left(\frac{t - t_\mathrm{lag}}{\mathrm{MDT}}%
#'   \right)^\beta}\right)}{f(t) = fmax(1 - exp(-((t - tlag)/MDT)^\beta))}
#' or
#' \deqn{f_t = f_\mathrm{max} \left(1 - %
#'   e^{-\frac{(t - t_\mathrm{lag})^\beta}{\alpha}}\right)}{f(t) = fmax%
#'   (1 - exp(-(((t - tlag)^\beta)/\alpha))),}
#' where \eqn{f_\mathrm{max}}{fmax} is the maximum dissolution,
#' \eqn{\mathrm{MDT}}{MDT} is the mean dissolution time,
#' \eqn{t_\mathrm{lag}}{tlag} is the lag time, \eqn{\alpha}{\alpha} and
#' \eqn{\beta}{\beta} are the scale and shape factor in Weibull function,
#' and \eqn{k}{k} is the rate constant in the first-order model. Obviously,
#' The two Weibull models are mathematically equivalent by letting
#' \eqn{\alpha = \mathrm{MDT}^\beta}{\alpha = MDT^\beta}.
#'
#' Individual model parameter were simulated according to the exponential
#' error model
#' \deqn{P_i = P_\mu e^{N\left(0, \sigma^2\right)}}{P(i) = P(\mu)%
#'   exp(N(0, \sigma^2)),}
#' where \eqn{P_i}{P(i)} and \eqn{P_\mu}{P(\mu)} denote the individual and
#' population model parameters; \eqn{N\left(0, \sigma^2\right)}{N(0, \sigma^2)}
#' represents the normal distribution with mean zero and variance \eqn{\sigma^2}
#' (\eqn{\sigma = \mathrm{CV}/100}{\sigma = CV/100}).
#'
#' ## How to supply `model.par`
#' The argument `model.par` should be supplied as a *named list* of model
#' parameters as explained above, and their respective CV for simulation of
#' individual parameters. Therefore, for the first-order model, three pairs of
#' parameters should be specified: `fmax/fmax.cv`, `k/k.cv`, and `tlag/tlag.cv`;
#' and for Weibull model, four pairs: `fmax/fmax.cv`, `tlag/tlag.cv`,
#' `beta/beta.cv`, and either `alpha/alpha.cv` or `mdt/mdt.cv`, depending on
#' the mathematical formula used. CV should be given in percentage, e.g., if
#' CV for `beta` is 30%, it should be specified as `beta.cv = 30`, *not*
#' `beta.cv = 0.3`. The order of the parameters does not matter, but the name
#' of the parameters must be given *exactly same* as described above.
#' For example:
#'   - `model.par = list(fmax = 100, fmax.cv = 5, k = 0.6, k.cv = 25,
#'      tlag = 0, tlag.cv = 0)` for the first-order model.
#'   - `model.par = list(fmax = 100, fmax.cv = 5, tlag = 5, tlag.cv = 10,
#'      mdt = 15, mdt.cv = 20, beta = 1.5, beta.cv = 5)`, or
#'   - `model.par = list(fmax = 100, fmax.cv = 5, tlag = 5, tlag.cv = 10,
#'      alpha = 60, alpha.cv = 20, beta = 1.5, beta.cv = 5)` for Weibull models.
#'
#' ## Notes on multivariate normal distribution approach
#' When this approach is used, depending on `dp/dp.cv`, there is no guarantee
#' that all individual profiles increase with time; near the end of the time
#' points, some individual profiles may decrease, especially when the
#' dissolution is close to the plateau and the variability is high. This can
#' happen to real life data, especially for those products with drug substances
#' that are unstable in dissolution medium. To force increase for all profiles,
#' set `ascending = TRUE`. Depending on the data, the program may take long
#' time to run, or may even fail.
#'
#' @examples
#' # time points
#' tp <- c(5, 10, 15, 20, 30, 45, 60)
#'
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
#' # using Weibull model to simulate 100 profiles without lag time
#' mod.par <- list(fmax = 100, fmax.cv = 2, tlag = 0, tlag.cv = 0,
#'                 mdt = 20, mdt.cv = 5, beta = 2.2, beta.cv = 5)
#' d2 <- sim.dp(tp, model.par = mod.par, seed = 100, n.units = 100L,
#'              plot = FALSE)
#'
#' # using multivariate normal distribution approach
#' # target mean profile with same length as tp
#' dp <- c(39, 56, 67, 74, 83, 90, 94)
#'
#' # CV% at each time point
#' dp.cv <- c(19, 15, 10, 8, 8, 5, 3)
#'
#' # simulation
#' d3 <- sim.dp(tp, dp = dp, dp.cv = dp.cv, seed = 1234, plot = FALSE)
#'
#' @export
sim.dp <- function(tp, dp, dp.cv, model = c("Weibull", "first-order"),
                   model.par = NULL, seed = NULL, n.units = 12L, product,
                   max.disso = 100, ascending = FALSE, message = FALSE,
                   time.unit = c("min", "h"), plot = TRUE,
                   plot.max.unit = 36L) {
  # initial checking -----------------------------------------------------------
  model     <- match.arg(model)
  time.unit <- match.arg(time.unit)

  # default tp if missing
  if (missing(tp)) {
    miss.tp <- TRUE
    if (time.unit == "min") {
      tp <- c(5, 10, 15, 20, 30, 45, 60)
    } else {
      tp <- c(1:4, seq(6, 12, 2))
    }
  } else miss.tp <- FALSE

  if (missing(dp)) use.model <- TRUE else use.model <- FALSE

  # if no seed is provided, generate one
  if (is.null(seed)) {
    seed <- sample(1:(.Machine$integer.max - 1), 1)
  }
  set.seed(seed)

  # if model.par is given by user, check for error
  if (!is.null(model.par)) {
    if (model == "first-order") {
      if (length(setdiff(c("fmax", "fmax.cv", "k", "k.cv", "tlag", "tlag.cv"),
                         names(model.par))) != 0) {
        stop("Model parameters are incorrect. Three pairs of parameters are\n",
             "required: 'fmax/fmax.cv', 'k/k.cv', and 'tlag/tlag.cv'.\n")
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
  }# end !is.null(model.par)

  # mathematical modelling approach when dp is missing -------------------------
  if (isTRUE(use.model)) {
    # remove time 0 if any
    tp <- tp[tp > 0]

    ## generate default model parameters if missing ------------------
    if (is.null(model.par)) {
      if (model == "first-order") {
        if (time.unit == "min") {# typical IR dosage form if unit is 'min'
          model.par <- list(fmax    = round(100*exp(rnorm(1, 0, 0.025)), 2),
                            fmax.cv = 3,
                            k       = round(0.1*exp(rnorm(1, 0, 0.3)), 6),
                            k.cv    = 30,
                            tlag    = 0,
                            tlag.cv = 0)
        } else {# typical ER dosage form if unit is 'h'
          model.par <- list(fmax    = round(98*exp(rnorm(1, 0, 0.03)), 2),
                            fmax.cv = 3,
                            k       = round(0.3*exp(rnorm(1, 0, 0.3)), 6),
                            k.cv    = 30,
                            tlag    = round(0.3*exp(rnorm(1, 0, 0.25)), 2),
                            tlag.cv = 20)
        }# end missing model.par for first-order model
      } else {# for Weibull model
        if (time.unit == "min") {# typical IR dosage form if unit is 'min'
          model.par <- list(fmax    = round(100*exp(rnorm(1, 0, 0.025)), 2),
                            fmax.cv = 3,
                            tlag    = 0, tlag.cv = 0,
                            mdt     = round(10*exp(rnorm(1, 0, 0.3)), 2),
                            mdt.cv  = 30,
                            beta    = round(0.8*exp(rnorm(1, 0, 0.2)), 6),
                            beta.cv = 20)
        } else {# typical ER dosage form if unit is 'h'
          model.par <- list(fmax = round(98*exp(rnorm(1, 0, 0.05)), 2),
                            fmax.cv = 5,
                            tlag = round(0.5*exp(rnorm(1, 0, 0.3)), 2),
                            tlag.cv = 30,
                            mdt = round(4*exp(rnorm(1, 0, 0.2)), 2),
                            mdt.cv = 10,
                            beta = round(2*exp(rnorm(1, 0, 0.3)), 6),
                            beta.cv = 20)
        }# end ER dosage form
      }# end missing model.par for Weibull model
      if (model.par$fmax > max.disso) model.par$fmax <- max.disso
    }# end is.null(model.par)

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
      dp <- model.par$fmax*(1 - exp(-model.par$k*tp.adj))

      # individual model par
      k.ind <- model.par$k*exp(rnorm(n.units, 0, model.par$k.cv/100))

      # individual dissolution data. each column is one unit
      dp.ind <- matrix(NA, nrow = length(tp), ncol = n.units)

      for (i in seq_len(n.units)) {
        tp.adj.ind <- tp - tlag.ind[i]
        tp.adj.ind[tp.adj.ind < 0] <- 0
        dp.ind[, i] <- fmax.ind[i]*(1 - exp(-k.ind[i]*tp.adj.ind))
      }

      # save individual parameters for output
      model.par.ind <- data.frame(
        fmax.ind = fmax.ind,
        k.ind    = k.ind,
        tlag.ind = tlag.ind,
        stringsAsFactors = FALSE
      )

      # output simulation info
      sim.info <- data.frame(
        fmax      = model.par$fmax,
        fmax.cv   = model.par$fmax.cv,
        k         = model.par$k,
        k.cv      = model.par$k.cv,
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

       alpha.ind <- model.par$alpha*exp(rnorm(n.units,0,model.par$alpha.cv/100))

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
        cat("- k         = ", model.par$k, "\n", sep = "")
        cat("- k.cv      = ", model.par$k.cv, "%\n\n", sep = "")
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
      cat("Seed number used: ", seed, "\n")
      if (!missing(dp.cv)) {
        warning("'dp.cv' is given without 'dp', so it is ignored.\n")
      }
    } # end message # end missing dp -----------------------
  } else {# with dp, using multivariate normal distribution approach -----------
    if (isTRUE(miss.tp)) {
      stop("Time points 'tp' must be provided.")
    }

    # check length, tp = dp = dp.cv
    if (length(tp) != length(dp)) {
      stop("Length of 'tp' and 'dp' must be equal. Check your data.")
    }

    if (!missing(dp.cv)) {
      if (length(dp) != length(dp.cv)) {
        stop("Length of 'dp' and 'dp.cv' must be equal. Check your data.")
      }
    }

    # remove time 0 from tp if any, and corresponding point from dp,
    # and if dp.cv is not missing, from dp.cv as well
    if (isTRUE((length(which(tp == 0)) != 0))) {
      tp <- tp[-(which(tp == 0))]
      dp <- dp[-(which(tp == 0))]
      if (!missing(dp.cv)) dp.cv[-(which(tp == 0))]
    }

    ## generate dp.cv if missing -------------------------------------
    # CV 10--20% for time points up to 10 min, 3--5% for dissolution between
    # 90--95%, 1--3% for dissolution > 95%, 5--10% for the rest points,
    if (missing(dp.cv)) {
      # initialize cv within 5--10% for all points.
      dp.cv <- rep(sample(5:10, 1), length(tp))

      # check if there is any time point <= 10, if so, cv within 10--20%
      if (time.unit == "min") {
        if (length(which(tp <= 10)) != 0) {
          dp.cv[which(tp <= 10)] <- sample(10:20, 1)
        }
      } else {# convert h to min
        if (length(which(tp*60 <= 10)) != 0) {
          dp.cv[which(tp*60 <= 10)] <- sample(10:20, 1)
        }
      }# end checking time points <= 10 min

      # if dp > 90%, then cv about 5%, which is more realistic
      dp.cv <- ifelse(dp > 95, sample(1:3, 1), dp.cv)
      dp.cv <- ifelse(dp >= 90 & dp <= 95, sample(3:5, 1), dp.cv)

      # lastly, check if there's dp = 0, if so, cv = 0.
      dp.cv[dp <= 0] <- 0
    }# end missing dp.cv

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
           "distribution, 'n.units' should not be less than number of time\n",
           "points 'tp'. Please either decrease the number of time points\n",
           "or increase 'n.units'. Alternatively, you can use the approach\n",
           "based on mathematical models.")
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
  if (isTRUE(use.model)) {
    dp.cv <- rep(NA, length(tp))
  } else {
    dp.cv <- c(0, dp.cv)
  }

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
  # class(sim.out) <- append(class(sim.out), "sim.dp.out")
  attr(sim.out, "come.from") <- "sim.dp"
  return(sim.out)
}# end function sim.dp
