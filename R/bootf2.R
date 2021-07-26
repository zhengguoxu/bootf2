#' Estimate 90% Confidence Intervals of f2 Using Bootstrap Methodology
#'
#' Main function to estimate 90% confidence intervals of f2 using bootstrap
#' methodology.
#'
#' @importFrom stats qnorm pnorm glm
#' @importFrom utils packageVersion sessionInfo
#' @usage
#' bootf2(test, ref, path.in, file.in, path.out, file.out,
#'        n.boots = 10000L, seed = 306, digits = 2L, alpha = 0.05,
#'        regulation = c("EMA", "FDA", "WHO"), min.points = 1L,
#'        both.TR.85 = FALSE, print.report = TRUE,
#'        report.style = c("concise", "intermediate", "detailed"),
#'        f2.type = c("all", "est.f2", "exp.f2", "bc.f2",
#'                    "vc.exp.f2", "vc.bc.f2"),
#'        ci.type = c("all", "normal", "basic", "percentile",
#'                    "bca.jackknife", "bca.boot"),
#'        quantile.type = c("all", 1:9, "boot"),
#'        jackknife.type = c("nt+nr", "nt*nr", "nt=nr"),
#'        time.unit = c("min", "h"), output.to.screen = TRUE,
#'        sim.data.out = FALSE)
#'
#' @param test,ref Data frames of dissolution profiles of test and reference
#'   product if \code{path.in} and \code{file.in} are not specified; otherwise,
#'   they should be character strings indicating the worksheet names of the
#'   Excel file where the dissolution data is saved. Required format:
#'   the first column should be time and the rest columns are dissolution
#'   data of each unit. See Input/Output in Details.
#' @param path.in,file.in,path.out,file.out Character strings of input and
#'   output directories and file names. See Input/Output in Details.
#' @param n.boots An integer indicating the number of bootstrap samples.
#' @param seed Numeric seed value for reproducibility. If missing, a random
#'   seed will be generated for reproducibility purpose.
#' @param digits An integer indicating the decimal points for the output.
#' @param alpha A numeric value between 0 and 1.
#'   \eqn{(1-2\times \alpha)\times 100}{(1 - 2*alpha)*100}
#'   confidence interval will be estimated.
#' @param regulation Character strings indicating regulatory guidelines.
#'   See Regulation in Details.
#' @param min.points An integer indicating the minimum time points to be used
#'   to calculate f2. See Regulation in Details.
#' @param digits An integer indicating the decimal points for the output.
#' @param both.TR.85 Logical. If \code{TRUE}, the old (incorrect)
#'   interpretation (that both the test and reference should release more than
#'   85%) will be used for f2 calculation when \code{regulation = 'FDA'}.
#'   See Regulation in Details.
#' @param print.report Logical. If \code{TRUE}, a plain text report will be
#'   produced.
#' @param report.style \code{'concise'} style will produce the confidence
#'   intervals for the f2 estimators; \code{'intermediate'} style will add
#'   a list of individual f2s for all bootstrap samples; \code{detailed}
#'   style will further add individual bootstrap samples along with their f2.
#'   See Input/Output in Details.
#' @param f2.type Character strings indicating which f2 estimators should be
#'   calculated.
#' @param ci.type Character strings indicating which type of  confidence
#'   intervals should be estimated.
#' @param quantile.type Character strings indicating the type of percentile.
#' @param jackknife.type Character strings indicating the type for jackknife
#'   method.
#' @param time.unit Character strings indicating the unit of time.
#'   It should be either \code{"min"} for minute or \code{"h"} for hour.
#'   It is mainly used for checking CV rules and making plot.
#'   See Regulation in Details.
#' @param output.to.screen Logical. If \code{TRUE}, a concise style summary
#'   report will be printed to screen.
#' @param sim.data.out Logical. If \code{TRUE}, all individual bootstrap data
#'   sets will be included in the output.
#' @returns A list of 3 or 5 components.
#'
#'   - \code{boot.ci}: A data frame of bootstrap f2 confidence intervals.
#'   - \code{boot.f2}: A data frame of all individual f2 values.
#'   - \code{boot.info}: A data frame with detailed information of bootstrap,
#'       such as the time points used for calculation of f2 and the number of
#'       \code{NA}s.
#'   - \code{boot.t} and \code{boot.r}: Individual bootstrap samples if
#'      \code{sim.data.out = TRUE}.
#'
#' @examples
#' \dontrun{
#' # time points
#' tp <- c(5, 10, 15, 20, 30, 45, 60)
#' # model.par for reference with low variability
#' par.r <- list(fmax = 100, fmax.cv = 3, mdt = 15, mdt.cv = 14,
#'               tlag = 0, tlag.cv = 0, beta = 1.5, beta.cv = 8)
#' # simulate reference data
#' dr <- sim.dp(tp, model.par = par.r, seed = 100, plot = FALSE)
#' # model.par for test
#' par.t <- list(fmax = 100, fmax.cv = 3, mdt = 12.29, mdt.cv = 12,
#'               tlag = 0, tlag.cv = 0, beta = 1.727, beta.cv = 9)
#' # simulate test data with low variability
#' dt <- sim.dp(tp, model.par = par.t, seed = 100, plot = FALSE)
#'
#' # bootstrap
#' bootf2(dt$sim.disso, dr$sim.disso, n.boots = 100)
#' }
#'
#' @export
bootf2 <- function(test, ref, path.in, file.in, path.out, file.out,
                   n.boots = 10000L, seed = 306L, digits = 2L, alpha = 0.05,
                   regulation = c("EMA", "FDA", "WHO"), min.points = 1L,
                   both.TR.85 = FALSE, print.report = TRUE,
                   report.style = c("concise", "intermediate", "detailed"),
                   f2.type = c("all", "est.f2", "exp.f2", "bc.f2",
                               "vc.exp.f2", "vc.bc.f2"),
                   ci.type = c("all", "normal", "basic", "percentile",
                               "bca.jackknife", "bca.boot"),
                   quantile.type = c("all", 1:9, "boot"),
                   jackknife.type = c("nt+nr", "nt*nr", "nt=nr"),
                   time.unit = c("min", "h"), output.to.screen = TRUE,
                   sim.data.out = FALSE) {
  # for output info.
  dt.name <- noquote(deparse1(substitute(test)))
  dr.name <- noquote(deparse1(substitute(ref)))

  # initial check --------------------------------------------------------------
  regulation     <- match.arg(regulation)
  report.style   <- match.arg(report.style)
  f2.type        <- match.arg(f2.type)
  ci.type        <- match.arg(ci.type)
  quantile.type  <- match.arg(quantile.type)
  jackknife.type <- match.arg(jackknife.type)
  time.unit      <- match.arg(time.unit)

  if (all(isTRUE(both.TR.85), regulation != "FDA")) {
    stop("'both.TR.85 = TRUE' is only valid when 'regulation = FDA'.")
  }

  if (any(missing(test), missing(ref))) {
    stop("Both 'test' and 'ref' have to be specified.")
  }

  if (quantile.type == "all") {
    q.type <- 1:9
  } else if (quantile.type == "boot") {
    q.type <- NULL
  } else {
    q.type <- as.numeric(quantile.type)
  }

  # read data ------------------------------------------------------------------
  if (all(missing(path.in), missing(file.in))) {
    data.t <- as.matrix(test, rownames.force = FALSE)
    data.r <- as.matrix(ref, rownames.force = FALSE)
    path.in <- file.in1 <- NA # for output name
  } else if (all(missing(path.in), !missing(file.in))) {
    stop("Please provide the directory 'path.in' where the file is stored.")
  } else {# for path.in not missing
    if (missing(file.in)) {
      stop("Please provide the name of the data file.")
    }

    # if path.in specified incorrectly
    if (!dir.exists(path.in)) {
      stop("The directory you specified does not exist. Check your spelling.")
    }

    path.in <- normalizePath(path.in, winslash = "/")
    path.in <- ifelse(regmatches(path.in, regexpr(".$", path.in)) == "/",
                      path.in, paste0(path.in, "/"))
    file.in1 <- file.in# for output name
    file.in <- paste0(path.in, file.in)
    if (!file.exists(file.in)) {
      stop(paste0("The file you specified does not exist. Don't forget to ",
                  "include\nthe extension 'xlsx' or 'xls' in the file name."))
    }

    # package readxl::read_excel
    data.t <- as.matrix(read_excel(file.in, test, col_types = "numeric"))
    data.r <- as.matrix(read_excel(file.in, ref, col_types = "numeric"))
  }# end read data

  nt <- NCOL(data.t) - 1
  nr <- NCOL(data.r) - 1

  if (is.null(seed)) {
    seed <- sample(1:(.Machine$integer.max - 1), 1)
  }
  set.seed(seed)

  # bootstrap data -------------------------------------------------------------
  # initialize result ------------------------------------------------
  boot.t <- boot.r  <- vector(mode = "list", length = n.boots)

  if (f2.type == "all") { # est.f2, exp.f2, bc.f2, vc.exp.f2, vc.bc.f2, tp
    boot.f2 <- matrix(0, nrow = n.boots, ncol = 6, byrow = TRUE,
                      dimnames = list(rep("", n.boots),
                                      c("est.f2", "exp.f2", "bc.f2",
                                        "vc.exp.f2", "vc.bc.f2", "tp")))
  } else { # f2, tp
    boot.f2 <- matrix(0, nrow = n.boots, ncol = 2, byrow = TRUE,
                      dimnames = list(rep("", n.boots), c(f2.type, "tp")))
  }

  # bootstrap index ------------------------------------------------------------
  # implement the same bootstrap algorithm as boot package to compare results
  # use safer version according to manual ?sample
  resample <- function(x, ...) x[sample.int(length(x), replace = TRUE, ...)]

  # function in boot need 1 data and 1 index option. each row is
  # considered as one multivariate observation so need to transpose it
  prod <- c(rep(1, nt), rep(2, nr))
  bt.array <- matrix(NA, nrow = n.boots, ncol = nt + nr)

  for (i in 1:2) {# fill by columns
    index <- seq_len(nt + nr)[prod == i]
    bt.array[, index] <- resample(index, n.boots*length(index))
  }

  for (i in 1:n.boots) {
    boot.t[[i]]  <- data.t[, c(1, bt.array[i, 1:nt] + 1)]
    boot.r[[i]]  <- data.r[, c(1, bt.array[i, (nt + 1):(nt + nr)] - nt + 1)]
    boot.f2[i, ] <- calcf2(test = boot.t[[i]], ref = boot.r[[i]],
                           regulation = regulation, digits = digits,
                           cv.rule = FALSE, min.points = min.points,
                           both.TR.85 = both.TR.85, message = FALSE,
                           f2.type = f2.type, plot = FALSE,
                           time.unit = time.unit)
  }

  # f2 original ------------------------------------------------------
  # calculate f2s with original data without regard to variability
  f2o <- calcf2(test = data.t, ref = data.r, regulation = regulation,
                digits = digits, cv.rule = FALSE, min.points = min.points,
                both.TR.85 = both.TR.85, message = FALSE, f2.type = f2.type,
                plot = FALSE, time.unit = time.unit)

  # interpolation function -------------------------------------------
  # function for interpolation of quantile. Davison, Ch5, Eq. 5.8.
  # modified from R boot package internal function 'norm.inter()'.
  normal.inter <- function(boot.f2, alpha) {#
    btf2  <- as.vector(boot.f2)
    btf2  <- btf2[is.finite(btf2)]
    n.f2  <- length(btf2)
    rk    <- (n.f2 + 1)*alpha
    k     <- trunc(rk)
    inds  <- seq_along(k)
    out   <- inds
    tstar <- sort(btf2)
    ints  <- (k == rk)
    if (any(ints)) out[inds[ints]] <- tstar[k[inds[ints]]]
    out[k == 0]    <- tstar[1L]
    out[k == n.f2] <- tstar[n.f2]
    not   <- function(v) xor(rep(TRUE, length(v)), v)
    temp  <- inds[not(ints) & k != 0 & k != n.f2]
    temp1 <- qnorm(alpha[temp])
    temp2 <- qnorm(k[temp]/(n.f2 + 1))
    temp3 <- qnorm((k[temp] + 1)/(n.f2 + 1))
    tk    <- tstar[k[temp]]
    tk1   <- tstar[k[temp] + 1L]
    out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2)*(tk1 - tk)
    return(out)
  }

  # confidence intervals -------------------------------------------------------

  # normal interval --------------------------------------------------
  # normal approximation: L,U = f2o - bias -/+ mean.err*z(1-alpha)
  # Davison, Bootstrap Methods and Their Application, CUP, 1997. Ch5,
  if (any(ci.type == "all", ci.type == "normal")) {
    normal.ci <- function(boot.f2, f2o, alpha) {
      btf2  <- as.vector(boot.f2)
      btf2  <- btf2[is.finite(btf2)]
      btf2.mean <- mean(btf2, na.rm = TRUE)
      btf2.var  <- var(btf2, na.rm = TRUE)
      2*f2o - btf2.mean + c(-1, 1)*sqrt(btf2.var)*qnorm(1 - alpha)
    }

    ci.normal <- data.frame(f2.type = NA, ci.type = NA,
                            ci.lower = NA, ci.upper = NA)

    for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
      ci.normal[i, 1]   <- dimnames(boot.f2)[[2]][i]
      ci.normal[i, 2]   <- "Normal"
      ci.normal[i, 3:4] <- normal.ci(boot.f2 = boot.f2[, i], f2o = f2o[[i]],
                                     alpha = alpha)
    }
  }# end normal CI

  # basic interval ---------------------------------------------------
  # Basic CI. Davidson, Ch5, Eq 5.6
  if (any(ci.type == "all", ci.type == "basic")) {
    basic.ci <- function(boot.f2, f2o, alpha) {
      btf2  <- as.vector(boot.f2)
      btf2  <- btf2[is.finite(btf2)]
      2*f2o - normal.inter(boot.f2 = btf2, alpha = alpha)
    }

    ci.basic <- data.frame(f2.type = NA, ci.type = NA,
                           ci.lower = NA, ci.upper = NA)

    for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
      ci.basic[i, 1]   <- dimnames(boot.f2)[[2]][i]
      ci.basic[i, 2]   <- "Basic"
      ci.basic[i, 3:4] <- basic.ci(boot.f2 = boot.f2[, i], f2o = f2o[i],
                                   alpha = c(1 - alpha, alpha))
    }
  } # end basic interval

  # percentile interval ----------------------------------------------
  if (any(ci.type == "all", ci.type == "percentile")) {
    ci.percentile <- data.frame(f2.type = NA, ci.type = NA,
                                ci.lower = NA, ci.upper = NA)

    if (quantile.type == "boot") {# same as boot package
      for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
        ci.percentile[i, 1]   <- dimnames(boot.f2)[[2]][[i]]
        ci.percentile[i, 2]   <- "Percentile (boot)  "
        ci.percentile[i, 3:4] <- normal.inter(boot.f2 = boot.f2[, i],
                                              alpha = c(alpha, 1 - alpha))
      }
    } else if (quantile.type == "all") {
      k <- 0
      for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
        for (j in seq_along(q.type)) {
          k <- k + 1
          ci.percentile[k, 1]   <- dimnames(boot.f2)[[2]][[i]]
          ci.percentile[k, 2]   <-
            ifelse(quantile.type == "all",
                   paste0("Percentile (Type ", j, ")"),
                   paste0("Percentile (Type ", q.type, ")"))

          ci.percentile[k, 3:4] <- quantile(boot.f2[, i],
                                            probs = c(alpha, 1 - alpha),
                                            na.rm = TRUE, names = FALSE,
                                            type = q.type[[j]])
        } # end quantile type loop j

        # Davison, Ch5, same as boot package percentile ci
        k <- k + 1
        ci.percentile[k, 1]   <- dimnames(boot.f2)[[2]][[i]]
        ci.percentile[k, 2]   <- "Percentile (boot)  "
        ci.percentile[k, 3:4] <- normal.inter(boot.f2 = boot.f2[, i],
                                              alpha = c(alpha, 1 - alpha))
      } # end f2 type loop i
    } else {# type = 1, 2, ..., 9
      for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
        ci.percentile[i, 1]   <- dimnames(boot.f2)[[2]][[i]]
        ci.percentile[i, 2]   <- paste0("Percentile (Type ", q.type, ")")
        ci.percentile[i, 3:4] <- quantile(boot.f2[, i],
                                          probs = c(alpha, 1 - alpha),
                                          na.rm = TRUE, names = FALSE,
                                          type = q.type)
      } # end quantile type loop j
    }
  } # end percentile CI

  # BCa, jackknife ---------------------------------------------------
  # bca CI, acceleration a by jackknife.
  # An introduction to the Bootstrap, by Efron B. and Tibshirani, R., 1993
  if (any(ci.type == "all", ci.type == "bca.jackknife")) {
    # function to obtain acceleration a by jackkinfe
    jackf2 <- function(data.t, data.r) {
      # remove time 0 point if any
      data.t <- data.t[data.t[, 1] != 0, ]
      data.r <- data.r[data.r[, 1] != 0, ]

      # initialize result
      # est.f2.jk <- exp.f2.jk <- bc.f2.jk <- rep(0, nt + nr)
      # f2.jk <- data.frame(est.f2.jk = NA, exp.f2.jk = NA,
      #                     bc.f2.jk = NA, tp.jk = NA)
      if (f2.type == "all") {
        n.col <- 6
        jk.f2.name <- c("est.f2", "exp.f2", "bc.f2",
                        "vc.exp.f2", "vc.bc.f2", "tp")
      } else {
        n.col <- 2
        jk.f2.name <- c(f2.type, "tp")
      }

      if (jackknife.type == "nt+nr") {
        jk.f2 <- matrix(NA, nrow = nt + nr, ncol = n.col,
                        dimnames = list(rep("", nt + nr), jk.f2.name))
        # jackknife with test data
        for (i in 1:nt) {
          jk.f2[i, 1:n.col] <-
            calcf2(test = data.t[, -(i + 1)], ref = data.r,
                   regulation = regulation, digits = digits, cv.rule = FALSE,
                   min.points = min.points, both.TR.85 = both.TR.85,
                   message = FALSE, f2.type = f2.type, plot = FALSE,
                   time.unit = time.unit)
        }
        # jackknife with reference data
        for (j in 1:nr) {
          jk.f2[i + j, 1:n.col] <-
            calcf2(test = data.t, ref = data.r[, -(j + 1)],
                   regulation = regulation, digits = digits, cv.rule = FALSE,
                   min.points = min.points, both.TR.85 = both.TR.85,
                   message = FALSE, f2.type = f2.type, plot = FALSE,
                   time.unit = time.unit)
        }
      } else if (jackknife.type == "nt*nr") {
        jk.f2 <- matrix(NA, nrow = nt*nr, ncol = n.col,
                        dimnames = list(rep("", nt*nr), jk.f2.name))

        k <- 0
        for (i in 1:nt) {
          for (j in 1:nr) {
            k <- k + 1
            jk.f2[k, 1:n.col] <-
              calcf2(test = data.t[, -(i + 1)], ref = data.r[, -(j + 1)],
                     regulation = regulation, digits = digits, cv.rule = FALSE,
                     min.points = min.points, both.TR.85 = both.TR.85,
                     message = FALSE, f2.type = f2.type, plot = FALSE,
                     time.unit = time.unit)
          }
        }
      } else {# in this case, need nt = nr
        if (!all.equal(nt, nr)) {# usu. nt = nr = 12, not a problem
          stop("To use this type of jackknife, the number of test and ",
               "reference data should be equal.")
        } else {
          jk.f2 <- matrix(NA, nrow = nt, ncol = n.col,
                          dimnames = list(rep("", nt), jk.f2.name))
          for (i in 1:nt) {
            jk.f2[i, 1:n.col] <-
              calcf2(test = data.t[, -(i + 1)], ref = data.r[, -(i + 1)],
                     regulation = regulation, digits = digits, cv.rule = FALSE,
                     min.points = min.points, both.TR.85 = both.TR.85,
                     message = FALSE, f2.type = f2.type, plot = FALSE,
                     time.unit = time.unit)
          }
        }
      }# end jackknife.type

      jk.f2.mean <- colMeans(jk.f2, na.rm = TRUE)[-n.col]

      # accelerated alpha and jackknife mean
      if (f2.type == "all") {
        a <- c(est.f2.a = NA, exp.f2.a = NA, bc.f2.a = NA,
               vc.exp.f2.a = NA, vc.bc.f2.a = NA)
      } else {
        aname <- paste0(f2.type, ".a")
        a <- c(NA)
        names(a) <- aname
      }

      for (i in 1:length(a)) {
        a[i] <- sum((jk.f2.mean[[i]] - jk.f2[, i])^3, na.rm = TRUE)/
          (6*(sum((jk.f2.mean[[i]] - jk.f2[, i])^2, na.rm = TRUE))^1.5)
      }
      jk.f2.mean <- jk.f2.mean[is.finite(a)]
      a <- a[is.finite(a)]
      return(c(a, jk.f2.mean))
    }# end fun jackf2

    # now get value a
    a.jack <- jackf2(data.t, data.r)

    bca.ci.jack <- function(boot.f2, f2o, a, alpha) {
      btf2 <- as.vector(boot.f2)
      btf2 <- btf2[is.finite(btf2)]
      n.f2 <- length(boot.f2)
      z0 <- qnorm(sum(btf2 < f2o, na.rm = TRUE)/n.f2)
      if (is.finite(z0)) {
        a1 <- pnorm(z0 + (z0 + qnorm(alpha))/(1 - a*(z0 + qnorm(alpha))))
        a2 <- pnorm(z0 + (z0 + qnorm(1 - alpha))/(1 - a*(z0 + qnorm(1-alpha))))
        return(normal.inter(boot.f2 = btf2, alpha = c(a1, a2)))
      } else return(NA)
    }

    ci.bca.jackknife <- data.frame(f2.type = NA, ci.type = NA,
                                   ci.lower = NA, ci.upper = NA)


    for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
      ci.bca.jackknife[i, 1]   <- dimnames(boot.f2)[[2]][i]
      ci.bca.jackknife[i, 2]   <- "BCa (Jackknife)"
      if (all(dimnames(boot.f2)[[2]][i] %in% names(a.jack),
              !is.na(f2o[[i]]))) {
        ci.bca.jackknife[i, 3:4] <-
          bca.ci.jack(boot.f2 = boot.f2[, i], f2o = f2o[[i]],
                      a = a.jack[[i]], alpha = alpha)
      } else {
        ci.bca.jackknife[i, 3:4] <- rep(NA, 2)
      }
    } # end f2 type loop i
  } # end BCa CI by jackknife

  # BCa, boot --------------------------------------------------------
  # BCa, by empirical regression, Davison, Sec 2.7.4
  # also, ref boot package, empinf.reg function. same as Mendyk's bootf2BCA
  if (any(ci.type == "all", ci.type == "bca.boot")) {
    bca.ci.boot <- function(boot.f2, bt.array, f2o, nt, nr, alpha) {
      btf2 <- as.vector(boot.f2)
      inds <- is.finite(boot.f2)
      btf2 <- btf2[inds]
      n.f2 <- length(boot.f2)
      # construct frequency table from bt.array
      table.f <- matrix(NA, nrow = n.f2, ncol = nt + nr)
      for (i in 1:n.f2) {
        tmp1 <- table(bt.array[i, ])
        tmp2 <- as.numeric(names(tmp1))
        tmp3 <- 1:(nt + nr) %in% tmp2
        tmp3[tmp2] <- tmp1
        table.f[i, ] <- tmp3
      }

      table.n <- matrix(c(rep(nt, nt), rep(nr, nr)), nrow = n.f2,
                        ncol = nt + nr, byrow = TRUE)

      # X is the dependent, response variable is the bootstrapped f2
      X <- (table.f/table.n)[, -c(1, nt + 1)]
      X <- X[inds, ]
      beta <- coef(glm(btf2 ~ X))[-1L]

      L <- rep(0, nt + nr)
      L[-c(1, nt + 1)] <- beta
      prod <- c(rep(1, nt), rep(2, nr))
      L <- L - tapply(L, prod, mean)[prod]

      # accelerate a
      a <- sum(L^3)/(6*sum(L^2)^1.5)

      z0 <- qnorm(sum(btf2 < f2o, na.rm = TRUE)/n.f2)
      if (is.finite(z0)) {
        a1 <- pnorm(z0 + (z0 + qnorm(alpha))/(1 - a*(z0 + qnorm(alpha))))
        a2 <- pnorm(z0 + (z0 + qnorm(1 - alpha))/(1 - a*(z0 + qnorm(1-alpha))))
        return(normal.inter(boot.f2 = btf2, c(a1, a2)))
      } else return(NA)
    } # end function bca.ci.boot

    ci.bca.boot <- data.frame(f2.type = NA, ci.type = NA,
                              ci.lower = NA, ci.upper = NA)

    for (i in seq_along(1:(NCOL(boot.f2) - 1))) {
      ci.bca.boot[i, 1] <- dimnames(boot.f2)[[2]][i]
      ci.bca.boot[i, 2] <- "BCa (boot)"
      if (is.na(f2o[[i]])) {
        ci.bca.boot[i, 3:4] <- rep(NA, 2)
      } else {
        ci.bca.boot[i, 3:4] <-
          bca.ci.boot(boot.f2 = boot.f2[, i], bt.array = bt.array,
                      f2o = f2o[[i]], nt = nt, nr = nr, alpha = alpha)
      }
    }
  } # end BCa CI by regression

  ## TODO: add mode CIs later

  # output results -------------------------------------------------------------
  # prepare output file ----------------------------------------------
  if (ci.type == "all") {
    boot.f2.ci <- rbind(ci.normal, ci.basic, ci.percentile,
                        ci.bca.jackknife, ci.bca.boot)
  } else if (ci.type == "normal") {
    boot.f2.ci <- ci.normal
  } else if (ci.type == "basic") {
    boot.f2.ci <- ci.basic
  } else if (ci.type == "percentile") {
    boot.f2.ci <- ci.percentile
  } else if (ci.type == "bca.jackknife") {
    boot.f2.ci <- ci.bca.jackknife
  } else {
    boot.f2.ci <- ci.bca.boot
  }

  # check output -----------------------------------------------------
  # if output file is missing, auto-generate one for user as
  # 'data.t_vs_data.r_time.zone_Date_HHMMSS.txt'
  if (isTRUE(print.report)) {
    if (missing(path.out)) {
      path.out <- getwd()
    } else {
      if(!dir.exists(path.out)) {
        stop("The directory you specified does not exist. Check your spelling.")
      }
    }
    path.out <- normalizePath(path.out, winslash = "/")
    path.out <- ifelse(regmatches(path.out, regexpr(".$", path.out)) == "/",
                       path.out, paste0(path.out, "/"))

    if (missing(file.out)) {
      file.out1 <- paste0(dt.name, "_vs_", dr.name, "_",
                          format(Sys.time(), "%Z"), "_",
                          format(Sys.Date(), "%F"), "_",
                          format(Sys.time(), "%H%M%S"), ".txt")
    } else {# if file.out provided with wrong extension
      file.out1 <-
        ifelse(regmatches(file.out, regexpr(".{3}$", file.out)) == "txt",
                          file.out, paste0(file.out, ".txt"))
    }
    file.out <- paste0(path.out, file.out1)
  } else {# no report
    path.out <- file.out <- NA
  }

  # output additional information for transparency -------------------
  boot.info <- data.frame(
    test                = dt.name,
    ref                 = dr.name,
    n.boots             = n.boots,
    seed                = seed,
    regulation          = regulation,
    min.points          = min.points,
    alpha               = alpha,
    digits              = digits,
    both.TR.85          = both.TR.85,
    f2.type             = f2.type,
    ci.type             = ci.type,
    quantile.type       = ifelse(ci.type %in% c("all", "percentile"),
                                 quantile.type, NA),
    jackknife.type      = ifelse(ci.type %in% c("all", "bca.jackknife"),
                                 jackknife.type, NA),
    time.unit           = time.unit,
    print.report        = print.report,
    report.style        = ifelse(isTRUE(print.report), report.style, NA),
    sim.data.out        = sim.data.out,
    tp1                 = sum(boot.f2[, "tp"] == 1),
    tp2                 = sum(boot.f2[, "tp"] == 2),
    t(f2o),
    path.in             = path.in,
    file.in             = file.in1,
    path.out            = path.out,
    file.out            = file.out1,
    stringsAsFactors = FALSE
  )

  if (f2.type == "all") {
    boot.info[["est.f2.na"]]    <- sum(is.na(boot.f2[, "est.f2"]))
    boot.info[["exp.f2.na"]]    <- sum(is.na(boot.f2[, "exp.f2"]))
    boot.info[["bc.f2.na"]]     <- sum(is.na(boot.f2[, "bc.f2"]))
    boot.info[["vc.exp.f2.na"]] <- sum(is.na(boot.f2[, "vc.exp.f2"]))
    boot.info[["vc.bc.f2.na"]]  <- sum(is.na(boot.f2[, "vc.bc.f2"]))
  } else if (f2.type == "est.f2") {
    boot.info[["est.f2.na"]] <- sum(is.na(boot.f2[, "est.f2"]))
  } else if (f2.type == "exp.f2") {
    boot.info[["exp.f2.na"]] <- sum(is.na(boot.f2[, "exp.f2"]))
  } else if (f2.type == "bc.f2") {
    boot.info[["bc.f2.na"]]  <- sum(is.na(boot.f2[, "bc.f2"]))
  } else if (f2.type == "vc.exp.f2") {
    boot.info[["vc.exp.f2.na"]] <- sum(is.na(boot.f2[, "vc.exp.f2"]))
  } else {
    boot.info[["vc.bc.f2.na"]] <- sum(is.na(boot.f2[, "vc.bc.f2"]))
  }

  # write report ---------------------------------------------------------------
  if (isTRUE(print.report)) {
    if (report.style == "concise") {
      sink(file.out, split = output.to.screen)
    } else {
      sink(file.out, split = FALSE)
    }

    # header information. for all styles -----------------------------
    cat("=================================================================\n")
    cat("|                                                               |\n")
    cat("|  Comparison of Dissolution Profiles by Bootstrap f2 Method.   |\n")
    cat("|_______________________________________________________________|\n")
    cat("|                                                               |\n")
    cat("| Smimilarity Criterion:                                        |\n")
    cat("| ----------------------                                        |\n")
    cat("|     To conclude similarity, the lower limit of 90% confidence |\n")
    cat("| interval should be greater than or equal to 50.               |\n")
    cat("|                                                               |\n")
    cat("=================================================================\n\n")
    # width for numeric values of ci
    cilo <- if (digits <= 2) {
      paste0("%-5", ".", digits, "f")
    } else {
      paste0("%-", digits + 3, ".", digits, "f")
    }

    ciup <- if (digits <= 2) {
      paste0("%5", ".", digits, "f")
    } else {
      paste0("%", digits + 3, ".", digits, "f")
    }

    # f2 and confidence intervals ------------------------------------
    # for all styles
    cat("============================================\n")
    cat("|              Main Results:               |\n")
    cat("|  f2 and Its Confidence Intervals (CI)s   |\n")
    cat("============================================\n\n")
    # for estimated f2 ----
    if (f2.type %in% c("all", "est.f2")) {
      ci.estf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "est.f2", ]
      ci.estf2.n <- NROW(ci.estf2)
      cat("----------------------\n")
      cat("* Estimated f2 Values\n")
      cat("----------------------\n")
      cat("  - with original data                :   ",
          round(f2o[["est.f2"]], digits), "\n", sep = "")
      if (ci.type %in% c("all", "bca.jackknife")) {
        cat("  - with original data (by Jackknife) :   ",
            round(a.jack[["est.f2"]], digits), "\n", sep = "")
      }
      cat("  - with bootstrapped data (Mean)     :   ",
          round(mean(boot.f2[, "est.f2"], na.rm = TRUE), digits),
          "\n", sep = "")
      cat("  - with bootstrapped data (Median)   :   ",
          round(median(boot.f2[, "est.f2"], na.rm = TRUE), digits),
          "\n\n", sep = "")

      cat("-----------------------\n")
      cat("* Confidence Intervals\n")
      cat("-----------------------\n")

      ci.header(digits = digits)

      for (i in 1:ci.estf2.n) {
        cat(sprintf(paste0("%23s", "%3s", cilo, "%3s", ciup),
                    ci.estf2[i, 2], "", ci.estf2[i, 3], "", ci.estf2[i, 4]),
            "\n", sep = "")
      }
      cat("  ----------------------------------------------------\n")
      cat("  Out of", n.boots, "bootstrapped data sets,\n")
      cat("  - Number of f2 calculated with 1 time point :   ",
          boot.info[["tp1"]], "\n", sep = "")
      cat("  - Number of f2 calculated with 2 time point :   ",
          boot.info[["tp2"]], "\n", sep = "")
      cat("  - Number of f2 cannot be calculated (NA)    :   ",
          boot.info[["est.f2.na"]], "\n", sep = "")
      cat("  ----------------------------------------------------\n")
      if (f2.type == "all") {
        cat("_____________________________________________________",
            "_________________\n\n", sep = "")
      } else {
        cat("\n\n")
      }
    } # end est.f2

    # for expected f2 ------------------------------------------------
    if (f2.type %in% c("all", "exp.f2")) {
      ci.expf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "exp.f2", ]
      ci.expf2.n <- NROW(ci.expf2)
      cat("---------------------\n")
      cat("* Expected f2 Values\n")
      cat("---------------------\n")
      cat("  - with original data                :   ",
          round(f2o[["exp.f2"]], digits), "\n", sep = "")

      if (ci.type %in% c("all", "bca.jackknife")) {
        cat("  - with original data (by Jackknife) :   ",
            round(a.jack[["exp.f2"]], digits), "\n", sep = "")
      }

      cat("  - with bootstrapped data (Mean)     :   ",
          round(mean(boot.f2[, "exp.f2"], na.rm = TRUE), digits),
          "\n", sep = "")
      cat("  - with bootstrapped data (Median)   :   ",
          round(median(boot.f2[, "exp.f2"], na.rm = TRUE), digits),
          "\n\n", sep = "")

      cat("-----------------------\n")
      cat("* Confidence Intervals\n")
      cat("-----------------------\n")

      ci.header(digits = digits)

      for (i in 1:ci.expf2.n) {
        cat(sprintf(paste0("%23s", "%3s", cilo, "%3s", ciup),
                    ci.expf2[i, 2], "", ci.expf2[i, 3], "", ci.expf2[i, 4]),
            "\n", sep = "")
      }

      cat("  ----------------------------------------------------\n")
      cat("  Out of", n.boots, "bootstrapped data sets,\n")
      cat("  - Number of f2 calculated with 1 time point :   ",
          boot.info[["tp1"]], "\n", sep = "")
      cat("  - Number of f2 calculated with 2 time point :   ",
          boot.info[["tp2"]], "\n", sep = "")
      cat("  - Number of f2 cannot be calculated (NA)    :   ",
          boot.info[["exp.f2.na"]], "\n", sep = "")
      cat("  ----------------------------------------------------\n")
      if (f2.type == "all") {
        cat("_____________________________________________________",
            "_________________\n\n", sep = "")
      } else {
        cat("\n\n")
      }
    } # end exp.f2

    # for bias-corrected f2 ------------------------------------------
    if (f2.type %in% c("all", "bc.f2")) {
      ci.bcf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "bc.f2", ]
      ci.bcf2.n <- NROW(ci.bcf2)

      cat("---------------------------\n")
      cat("* Bias-Corrected f2 Values\n")
      cat("---------------------------\n")
      cat("  - with original data                :   ",
          round(f2o[["bc.f2"]], digits), "\n", sep = "")

      if (ci.type %in% c("all", "bca.jackknife")) {
        cat("  - with original data (by Jackknife) :   ",
            ifelse ("bc.f2" %in% names(a.jack),
                    round(a.jack[["bc.f2"]], digits), NA), "\n", sep = "")
      }

      cat("  - with bootstrapped data (Mean)     :   ",
          round(mean(boot.f2[, "bc.f2"], na.rm = TRUE), digits),
          "\n", sep = "")
      cat("  - with bootstrapped data (Median)   :   ",
          round(median(boot.f2[, "bc.f2"], na.rm = TRUE), digits),
          "\n\n", sep = "")

      cat("-----------------------\n")
      cat("* Confidence Intervals\n")
      cat("-----------------------\n")

      ci.header(digits = digits)

      for (i in 1:ci.bcf2.n) {
        cat(sprintf(paste0("%23s", "%3s", cilo, "%3s", ciup),
                    ci.bcf2[i, 2], "", ci.bcf2[i, 3], "", ci.bcf2[i, 4]),
            "\n", sep = "")
      }
      cat("  ----------------------------------------------------\n")
      cat("  Out of", n.boots, "bootstrapped data sets,\n")
      cat("  - Number of f2 calculated with 1 time point :   ",
          boot.info[["tp1"]], "\n", sep = "")
      cat("  - Number of f2 calculated with 2 time point :   ",
          boot.info[["tp2"]], "\n", sep = "")
      cat("  - Number of f2 cannot be calculated (NA)    :   ",
          boot.info[["bc.f2.na"]], "\n", sep = "")
      cat("  ----------------------------------------------------\n")
      if (f2.type == "all") {
        cat("_____________________________________________________",
            "_________________\n\n", sep = "")
      } else {
        cat("\n\n")
      }
    } # end bc.f2

    # for variance corrected expected f2 ----
    if (f2.type %in% c("all", "vc.exp.f2")) {
      ci.vcexpf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "vc.exp.f2", ]
      ci.vcexpf2.n <- NROW(ci.vcexpf2)
      cat("----------------------------------------\n")
      cat("* Variance-corrected Expected f2 Values\n")
      cat("----------------------------------------\n")
      cat("  - with original data                :   ",
          round(f2o[["vc.exp.f2"]], digits), "\n", sep = "")

      if (ci.type %in% c("all", "bca.jackknife")) {
        cat("  - with original data (by Jackknife) :   ",
            round(a.jack[["vc.exp.f2"]], digits), "\n", sep = "")
      }


      cat("  - with bootstrapped data (Mean)     :   ",
          round(mean(boot.f2[, "vc.exp.f2"], na.rm = TRUE), digits),
          "\n", sep = "")
      cat("  - with bootstrapped data (Median)   :   ",
          round(median(boot.f2[, "vc.exp.f2"], na.rm = TRUE), digits),
          "\n\n", sep = "")

      cat("-----------------------\n")
      cat("* Confidence Intervals\n")
      cat("-----------------------\n")

      ci.header(digits = digits)

      for (i in 1:ci.vcexpf2.n) {
        cat(sprintf(paste0("%23s", "%3s", cilo, "%3s", ciup), ci.vcexpf2[i, 2],
                    "", ci.vcexpf2[i, 3], "", ci.vcexpf2[i, 4]),
            "\n", sep = "")
      }

      cat("  ----------------------------------------------------\n")
      cat("  Out of", n.boots, "bootstrapped data sets,\n")
      cat("  - Number of f2 calculated with 1 time point :   ",
          boot.info[["tp1"]], "\n", sep = "")
      cat("  - Number of f2 calculated with 2 time point :   ",
          boot.info[["tp2"]], "\n", sep = "")
      cat("  - Number of f2 cannot be calculated (NA)    :   ",
          boot.info[["vc.exp.f2.na"]], "\n", sep = "")
      cat("  ----------------------------------------------------\n")
      if (f2.type == "all") {
        cat("_____________________________________________________",
            "_________________\n\n", sep = "")
      } else {
        cat("\n\n")
      }
    } # end vc.exp.f2

    # for variance- and bias-corrected f2 ----
    if (f2.type %in% c("all", "vc.bc.f2")) {
      ci.vcbcf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "vc.bc.f2", ]
      ci.vcbcf2.n <- NROW(ci.vcbcf2)

      cat("-----------------------------------------\n")
      cat("* Variance- and Bias-Corrected f2 Values\n")
      cat("-----------------------------------------\n")
      cat("  - with original data                :   ",
          round(f2o[["vc.bc.f2"]], digits), "\n", sep = "")

      if (ci.type %in% c("all", "bca.jackknife")) {
        cat("  - with original data (by Jackknife) :   ",
            ifelse ("vc.bc.f2" %in% names(a.jack),
                    round(a.jack[["vc.bc.f2"]], digits), NA), "\n", sep = "")
      }

      cat("  - with bootstrapped data (Mean)     :   ",
          round(mean(boot.f2[, "vc.bc.f2"], na.rm = TRUE), digits),
          "\n", sep = "")
      cat("  - with bootstrapped data (Median)   :   ",
          round(median(boot.f2[, "vc.bc.f2"], na.rm = TRUE), digits),
          "\n\n", sep = "")

      cat("-----------------------\n")
      cat("* Confidence Intervals\n")
      cat("-----------------------\n")

      ci.header(digits = digits)

      for (i in 1:ci.vcbcf2.n) {
        cat(sprintf(paste0("%23s", "%3s", cilo, "%3s", ciup),
                    ci.vcbcf2[i, 2], "", ci.vcbcf2[i, 3], "", ci.vcbcf2[i, 4]),
            "\n", sep = "")
      }
      cat("  ----------------------------------------------------\n")
      cat("  Out of", n.boots, "bootstrapped data sets,\n")
      cat("  - Number of f2 calculated with 1 time point :   ",
          boot.info[["tp1"]], "\n", sep = "")
      cat("  - Number of f2 calculated with 2 time point :   ",
          boot.info[["tp2"]], "\n", sep = "")
      cat("  - Number of f2 cannot be calculated (NA)    :   ",
          boot.info[["vc.bc.f2.na"]], "\n", sep = "")
      cat("  ----------------------------------------------------\n")
      if (f2.type == "all") {
        cat("_____________________________________________________",
            "_________________\n\n", sep = "")
      } else {
        cat("\n\n")
      }
    } # end vc.bc.f2

    # system and program info ----------------------------------------
    # output function settings so anyone who read the results can reproduce it
    cat("============================================\n")
    cat("| Function Settings and System Information |\n")
    cat("============================================\n\n")
    cat("---------------------\n")
    cat("* Function Settings\n")
    cat("---------------------\n")
    cat("  - test              :   ", dt.name, "\n", sep = "")
    cat("  - ref               :   ", dr.name, "\n", sep = "")
    cat("  - n.boots           :   ", n.boots, "\n", sep = "")
    cat("  - seed              :   ", seed, "\n", sep = "")
    cat("  - digits            :   ", digits, "\n", sep = "")
    cat("  - alpha             :   ", alpha, " (", (1 - 2*alpha)*100, "% CI)\n",
        sep = "")
    cat("  - regulation        :   ", regulation, "\n", sep = "")
    cat("  - min.points        :   ", min.points, "\n", sep = "")
    cat("  - both.TR.85        :   ", both.TR.85, "\n", sep = "")
    cat("  - print.report      :   ", print.report, "\n", sep = "")
    cat("  - report.style      :   ", report.style, "\n", sep = "")
    cat("  - f2.type           :   ", f2.type, "\n", sep = "")
    cat("  - ci.type           :   ", ci.type, "\n", sep = "")
    cat("  - quantile.type     :   ", quantile.type, "\n", sep = "")
    cat("  - jackknife.type    :   ", jackknife.type, "\n", sep = "")
    cat("  - time.unit         :   ", time.unit, "\n", sep = "")
    cat("  - output.to.screen  :   ", output.to.screen, "\n", sep = "")
    cat("  - sim.data.out      :   ", sim.data.out, "\n", sep = "")
    cat("  - path.in           :   ", path.in, "\n", sep = "")
    cat("  - file.in           :   ", file.in1, "\n", sep = "")
    cat("  - path.out          :   ", path.out, "\n", sep = "")
    cat("  - file.out          :   ", file.out1, "\n\n", sep = "")

    cat("---------------------\n")
    cat("* System Information\n")
    cat("---------------------\n")

    rinfo   <- strsplit(sessionInfo()$R.version$version.string, " ")
    sysinfo <- Sys.info()
    bootf2v <- as.character(packageVersion("bootf2"))
    cat("  - Operating System Name     :   ",
        sysinfo["sysname"], " ", sysinfo["release"], "\n", sep = "")
    cat("  - Operating System Version  :   ",
        sysinfo["version"], "\n", sep = "")
    cat("  - Machine Node Name         :   ",
        sysinfo["nodename"], "\n", sep = "")
    cat("  - User Name                 :   ",
        sysinfo["user"], "\n", sep = "")
    cat("  - Time Zone                 :   ",
        Sys.timezone(), "\n", sep = "")
    cat("  - R Version                 :   ",
        rinfo[[1]][3], " ", rinfo[[1]][4], "\n", sep = "")
    cat("  - Package bootf2 Version    :   ", bootf2v, "\n", sep = "")
    cat("_____________________________________________________",
        "_________________\n\n", sep = "")
    cat("The current report was generated at ", format(Sys.time(), "%H:%M:%S"),
        " on ", format(Sys.time(), "%F"), " ", format(Sys.time(), "%Z"), " by",
        "\nuser '", sysinfo["user"], "' on machine '", sysinfo["nodename"],
        "', and saved as text file\n'", file.out1, "' at the location:\n'",
        path.out, "'.\n", sep = "")
    cat("=============================================================",
        "=========\n\n", sep = "")
    # end concise report style.

    if (report.style %in%  c("intermediate", "detailed")) {
      cat("Individual f2 for each bootstrapped data set\n")
      cat("___________________________________________________________",
          "___________\n\n", sep = "")
      if (f2.type == "all") {
        if (digits <= 2) {
          cat("Bootstrap  est.f2  exp.f2  bc.f2   vc.exp.f2  vc.bc.f2  ",
              "Time Points\n", sep = "")
          for (i in 1:n.boots) {
            cat(sprintf(paste0("%2s", "%-7s", "%2s",
                               paste0("%-6", ".", digits, "f"), "%2s",
                               paste0("%-6", ".", digits, "f"), "%2s",
                               paste0("%-6", ".", digits, "f"), "%2s",
                               paste0("%-9", ".", digits, "f"), "%2s",
                               paste0("%-8", ".", digits, "f"), "%7s", "%-s"),
                        "", i, "", boot.f2[i, 1], "", boot.f2[i, 2], "",
                        boot.f2[i, 3], "", boot.f2[i, 4], "", boot.f2[i, 5],
                        "", boot.f2[i, 6]),
                "\n", sep = "")
          }
        } else { # digits > 2
          cat("Bootstrap  est.f2", rep(" ", digits - 1), "exp.f2",
              rep(" ", digits - 1), "bc.f2 ", rep(" ", digits - 1),
              "vc.exp.f2", rep(" ", digits - 1), "vc.bc.f2",
              rep(" ", digits - 1), "Time Points\n", sep = "")
          for (i in 1:n.boots) {
            cat(sprintf(paste0("%2s", "%-7s", "%2s",
                               paste0("%-", digits+3, ".", digits, "f"), "%2s",
                               paste0("%-", digits+3, ".", digits, "f"), "%2s",
                               paste0("%-", digits+3, ".", digits, "f"), "%2s",
                               paste0("%-", digits+6, ".", digits, "f"), "%2s",
                               paste0("%-", digits+5, ".", digits, "f"), "%7s",
                               "%-s"), "", i, "", boot.f2[i, 1], "",
                        boot.f2[i, 2], "", boot.f2[i, 3], "", boot.f2[i, 4], "",
                        boot.f2[i, 5], "", boot.f2[i, 6]), "\n", sep = "")
          }
        } # end digits > 2
      } else { # f2.type individual
        if (digits <= 2) {
          cat("Bootstrap No.   ", f2.type, "   No. Time Points\n", sep = "")
          for (i in 1:n.boots) {
            cat(sprintf(paste0("%4s", "%-9s", "%3s",
                               paste0("%-6", ".", digits, "f"), "%3s", "%-s"),
                        "", i, "", boot.f2[i, 1], "", boot.f2[i, 2]),
                "\n", sep = "")
          }
        } else { # digits > 2
          cat("Bootstrap No.   ", f2.type, rep(" ", digits + 1),
              "No. Time Points\n", sep = "")
          for (i in 1:n.boots) {
            cat(sprintf(paste0("%4s", "%-9s", "%3s",
                               paste0("%-", digits + 4, ".", digits, "f"),
                               "%7s", "%-s"), "", i, "", boot.f2[i, 1], "",
                        boot.f2[i, 4]),
                "\n", sep = "")
          }
        } # end digits > 2
      } # end f2.type
      cat("___________________________________________________________",
          "___________\n\n", sep = "")
    } # end intermediate

    # this part only for detailed report. long process time. not recommended
    if (report.style == "detailed") {

      cat("\n=====================================================",
          "===========================\n\n", sep = "")
      cat("Individual bootstrapped data set and its f2s\n")
      cat("_______________________________________________________",
          "_________________________\n\n", sep = "")
      cat("Original data set for test\n--------------------------\n")
      print(as.data.frame(data.t), row.names = FALSE)
      cat("\nOriginal data set for reference\n")
      cat("-------------------------------\n")
      print(as.data.frame(data.r), row.names = FALSE)

      tpo <- data.r[, 1][data.r[, 1] > 0][1:f2o[["tp"]]]

      if (f2.type %in% c("all", "est.f2")) {
        cat("\n\nEstimated f2                    :   ", f2o[["est.f2"]],
            "\n", sep = "")
      }

      if (f2.type %in% c("all", "exp.f2")) {
        cat("Expected f2                     :   ", f2o[["exp.f2"]],
            "\n", sep = "")
      }

      if (f2.type %in% c("all", "bc.f2")) {
        cat("Bias-corrected f2               :   ", f2o[["bc.f2"]],
            "\n", sep = "")
      }

      if (f2.type %in% c("all", "vc.exp.f2")) {
        cat("Variance-corrected expected f2  :   ", f2o[["vc.exp.f2"]],
            "\n", sep = "")
      }

      if (f2.type %in% c("all", "vc.bc.f2")) {
        cat("Variance- and bias-corrected f2 :   ", f2o[["vc.bc.f2"]],
            "\n", sep = "")
      }

      cat("Time Points Used                :   ", f2o[["tp"]], " (",
          paste0(tpo[1:length(tpo) - 1], rep(", ", length(tpo) - 1)), "and ",
          tpo[length(tpo)], " ", time.unit, ")\n", sep = "")
      cat("-------------------------------------------------------",
          "-------------------------\n\n", sep = "")

      for (i in 1:n.boots) {
        dimnames(boot.t[[i]])[[2]] <-
          c("time", paste0("U", 1:(NCOL(boot.t[[i]]) - 1)))

        dimnames(boot.r[[i]])[[2]] <-
          c("time", paste0("U", 1:(NCOL(boot.r[[i]]) - 1)))

        tpb <- boot.t[[i]][, 1][boot.t[[i]][, 1] > 0][1:boot.f2[i, "tp"]]

        cat("Bootstrap data set for test       :       ", i, "\n", sep = "")
        cat("-----------------------------------\n")
        print(boot.t[[i]], row.names = FALSE)
        cat("\n")
        cat("Bootstrap data set for reference  :       ", i, "\n", sep = "")
        cat("-----------------------------------\n")
        print(boot.r[[i]], row.names = FALSE)

        if (f2.type %in% c("all", "est.f2")) {
          cat("\n\nEstimated f2                    :   ",
              boot.f2[[i, "est.f2"]], "\n", sep = "")
        }

        if (f2.type %in% c("all", "exp.f2")) {
          cat("Expected f2                     :   ",
              boot.f2[[i, "exp.f2"]], "\n", sep = "")
        }

        if (f2.type %in% c("all", "bc.f2")) {
          cat("Bias-corrected f2               :   ",
              boot.f2[[i, "bc.f2"]], "\n", sep = "")
        }

        if (f2.type %in% c("all", "vc.exp.f2")) {
          cat("Variance-corrected expected f2  :   ",
              boot.f2[[i, "vc.exp.f2"]], "\n", sep = "")
        }

        if (f2.type %in% c("all", "vc.bc.f2")) {
          cat("Variance- and bias-corrected f2 :   ",
              boot.f2[[i, "vc.bc.f2"]], "\n", sep = "")
        }

        cat("Time Points Used                :   ", boot.f2[[i, "tp"]], " (",
            paste0(tpo[1:length(tpo) - 1], rep(", ", length(tpo) - 1)), "and ",
            tpo[length(tpo)], " ", time.unit, ")\n", sep = "")
        cat("-------------------------------------------------------",
            "-------------------------\n\n", sep = "")
      } # end loop for i data set
      cat("\n=====================================================",
          "===========================\n\n", sep = "")
    } # end detailed report

    sink()
  } # end isTRUE(print.report)

  if (isTRUE(sim.data.out)) {
    invisible(list(boot.ci = boot.f2.ci[order(boot.f2.ci[, 1]), ],
                   boot.f2 = as.data.frame(boot.f2, row.names = ""),
                   boot.info = boot.info, boot.t = boot.t, boot.r = boot.r))
  } else {# save some memory
    invisible(list(boot.ci = boot.f2.ci[order(boot.f2.ci[, 1]), ],
                   boot.f2 = as.data.frame(boot.f2, row.names = ""),
                   boot.info = boot.info))
  }
}
