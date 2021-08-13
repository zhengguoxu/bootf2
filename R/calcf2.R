#' Calculate Similarity Factor \eqn{f_2}{f2}
#'
#' Main function to calculate \eqn{f_2}{f2} according to different regulatory
#' guidelines.
#'
#' @importFrom readxl read_excel excel_sheets
#' @importFrom stats quantile reshape var rnorm median quantile lm coef
#' @import ggplot2
#' @usage
#' calcf2(test, ref, path.in, file.in, path.out, file.out,
#'        regulation = c("EMA", "FDA", "WHO", "Canada", "ANVISA"),
#'        cv.rule = TRUE, message = FALSE, min.points = 3L,
#'        f2.type = c("est.f2", "exp.f2", "bc.f2", "vc.exp.f2",
#'                    "vc.bc.f2", "all"), both.TR.85 = FALSE,
#'        digits = 2L, time.unit = c("min", "h"),  plot = TRUE,
#'        plot.start.time = 0, plot.max.unit = 24L)
#'
#' @param test,ref *Data frames* of dissolution profiles of test and reference
#'   product if `path.in` and `file.in` are not specified; otherwise, they
#'   should be *character* strings indicating the worksheet names of the Excel
#'   file where the dissolution data is saved. See Input/Output in Details.
#' @param path.in,file.in,path.out,file.out *Character* strings of input and
#'   output directories and file names. See Input/Output in Details.
#' @param regulation *Character* strings indicating regulatory guidelines. See
#'   Regulation in Details.
#' @param cv.rule *Logical*. If `TRUE`, CV rule will be checked according
#'   to regulatory guidelines. See Regulation in Details.
#' @param message *Logical*. If `TRUE`, the results and messages will be
#'   printed on screen. Users are recommended to set it to `TRUE`.
#' @param min.points An *integer* indicating the minimum time points to be used
#'   to calculate \eqn{f_2}{f2}. The default value 3 should be used for
#'   conventional \eqn{f_2}{f2} calculation. This parameter is mainly used for
#'   bootstrap \eqn{f_2}{f2} method. See Regulation in Details.
#'   @seealso [bootf2()].
#' @param f2.type *Character* strings indicating which \eqn{f_2}{f2} estimators
#'   should be calculated. For conventional \eqn{f_2}{f2} calculation, the
#'   default `"est.f2"` should be used. Other estimators are mainly for the
#'   bootstrap method. @seealso [bootf2()].
#' @param both.TR.85 *Logical*. If `TRUE`, and if `regulation = "FDA"`, all
#'   measurements up to the time points at which both test and reference
#'   products dissolve more than 85% will be used to calculate \eqn{f_2}{f2}.
#'   This is the conventional, but incorrect, interpretation of the US FDA rule.
#'   Therefore, the argument should only be set to `TRUE` for validation purpose
#'   such as comparing the results from old literature that use the wrong
#'   interpretation to calculate \eqn{f_2}{f2}. See Regulation in Details.
#' @param digits An *integer* indicating the decimal points for the output.
#' @param time.unit *Character* strings indicating the unit of time. It should
#'   be either `"min"` for minute or `"h"` for hour. It is mainly used for
#'   checking CV rules and making plot. See Regulation in Details.
#' @param plot *Logical*. If `TRUE`, a dissolution versus time plot will be
#'   printed.
#' @param plot.start.time *Numeric* value indicating the starting time for the
#'   plot.
#' @param plot.max.unit *Integer*. If the number of individual units is no more
#'   than this value, the mean and all individual profiles will be plotted;
#'   otherwise, individual profiles will be represented by boxplots at each
#'   time point. Therefore, to avoid overplotting, this value should not be
#'   too large. @seealso [calcf2()].
#'
#' @return A *data frame* of \eqn{f_2}{f2} type and \eqn{f_2}{f2} value, the
#'   number of time points used for the calculation (`f2.tp`), indication if
#'   both test and reference dissolve more than 85% at 15 min (`d85at15`), and
#'   other information used for the calculation.
#'
#' @details
#' ## Minimum required arguments that must be provided by the user
#' Arguments `test` and `ref` must be provided by the user. They should be `R`
#' `data frames`, with *time as the first column*, and all individual profiles
#' profiles as the rest columns, or mean profile as the second column if only
#' mean profile is available. The actual names of the columns do not matter
#' since they will be renamed internally.
#'
#' ## Input/Output
#' The dissolution data of test and reference product can either be provided as
#' *data frames* for `test` and `ref`, as explained above, or be read from an
#' *Excel file* with data of test and reference stored in *separate worksheets*.
#' In the latter case, the argument `path.in`, the directory where the Excel
#' file is, and `file.in`, the name of the Excel file *including the file
#' extension `.xls` or `.xlsx`*, must be provided. In such case, the argument
#' `test` and `ref` must be *the names of the worksheets in quotation marks*.
#' The first column of each Excel worksheet must be time, and the rest columns
#' are individual dissolution profiles, or the second column must be mean
#' profile if only mean data is available. The first row should be column names,
#' such as time, unit01, unit02, ... The actual names of the columns do not
#' matter as they will be renamed internally.
#'
#' Arguments `path.out` and `file.out` are the names of the output directory
#' and file. It is an overkill to output such simple calculations; therefore,
#' unless these two arguments are specified by the user, results are printed
#' on screen by default.
#'
#' ## Regulation
#'
#' To apply \eqn{f_2}{f2} method, different regulatory guidelines have slightly
#' different requirements. Some requirements are almost universal, such as same
#' time points for the test and reference product, minimum 3 time points
#' (excluding time zero), and twelve individual profiles for each formulation.
#' Other requirements are slightly different among different regulatory
#' guidelines, or at least interpreted differently. Two main issues are the
#' rules for the variability (CV Rule) and time points where dissolution is more
#' than 85% (85% Rule).
#'
#' ### CV rule
#' - `EMA`, `Canada`, and `ANVISA`: The CV of the *first time point* should not
#'   be greater than 20%, and the CV of the rest time points should not be
#'   greater than 10%.
#' - `WHO`: The CV should not be greater than 20% for *time points up to
#'   10 min*, and not greater than 10% for the rest time points.
#' - `FDA`: US FDA is more flexible. The CV for the *early time points* should
#'   not be greater than 20%, and for the rest time points, not greater than
#'   10%.
#'
#' The phrase *the first time point* in `EMA` rule was later interpreted as all
#' time points up to 10 min, according to an unofficial communication with an
#' European regulator. This makes the *`EMA` rule the same as `WHO` rule*. For
#' example, if there are 5 min and 10 min time points in the dissolution
#' profiles, the CV for both 5 min and 10 min should not be greater than 20%.
#'
#' The *first time point* in `ANVISA` rule corresponds to *40% of the total
#' collected points*. For example, for a dissolution profile with five
#' collection times, the first two collection times are considered first points.
#'
#' The phrase *early time points* in `FDA` rule is typically interpreted as
#' those points up to 15 min, sometimes even up to 20 min according to
#' an unofficial communication with FDA staff. In the function `calcf2()`, the
#' cutting point for FDA rule is 15 min.
#'
#' ### 85% Rule
#' This rule is implemented as follows:
#' - `EMA`, `FDA`, `Canada`, and `ANVISA`: Only one measurement is considered
#'   after 85% of dissolution for any product.
#' - `WHO`: Dissolution profiles should be 'cut' at the time point where
#'   the reference release more than 85%. Therefore, `WHO` rule only differs
#'   from rule of `EMA`, `FDA`, `Canada`, and `ANVISA` when test product
#'   dissolve faster than reference. If reference product dissolve faster, then
#'   rules of all five regulatory bodies are same in this regard.
#'
#' ### Notes on conventional FDA rule
#' The exact phrase in the guidance of US FDA regarding this rule is that
#' "*Only one measurement should be considered after 85% dissolution of both
#' the products*." Due to the ambiguous word "both" used in the sentence, the
#' conventional interpretation was that all measurements up to the time point
#' at which both test and reference dissolved more than 85% should be included
#' in the calculation of \eqn{f_2}{f2}. However, this is only true when both
#' test and reference dissolve more than 85% at the same time points.
#'
#' Consider the following example:
#'
#' |time |test |reference|
#' |---: |---: |---:     |
#' | 5   | 7   |10       |
#' |10   |15   |20       |
#' |15   |50   |55       |
#' |20   |69   |86       |
#' |30   |82   |90       |
#' |45   |84   |95       |
#' |60   |86   |97       |
#'
#' According to conventional interpretation, all measurements up to 60 min
#' should be included to calculate \eqn{f_2}{f2} because both test and reference
#' dissolved more than 85% only at 60 min, not at any earlier time point.
#' However, in such case, there would be 4 measurement of reference (20, 30, 45,
#' and 60 min) included in the calculation, which would be a direct
#' contradictory to the phrase "Only *one measurement* should be considered
#' after 85% ..." in the same statement in the guidance!
#'
#' In an unofficial communication using this example, an FDA staff confirmed
#' that only the first 4 time points (up to 20 min) would be used. In other
#' words, *FDA rule in this regard is the same as EMA rule*.
#'
#' The statement in `ANVISA` guideline also uses the word "ambos" (means both),
#' which could also lead to the similar confusion. Follow the same logic as
#' demonstrated above, it should also be interpreted as the same rule in EMA
#' guideline.
#'
#' Read vignette *Introduction to bootf2* for more details.
#'
#' @examples
#' tp <- c(5, 10, 15, 20, 30, 45, 60)
#'
#' mod.par.t <- list(fmax = 100, fmax.cv = 2, tlag = 0, tlag.cv = 0,
#'                   mdt = 20, mdt.cv = 5, beta = 2.2, beta.cv = 5)
#'
#' d.t <- sim.dp(tp, model.par = mod.par.t, seed = 100, n.units = 120L,
#'               plot = FALSE)$sim.disso
#'
#' mod.par.r <- list(fmax = 100, fmax.cv = 2, tlag = 0, tlag.cv = 0,
#'                   mdt = 25, mdt.cv = 4, beta = 2.1, beta.cv = 3)
#'
#' d.r <- sim.dp(tp, model.par = mod.par.r, seed = 100, n.units = 120L,
#'               plot = FALSE)$sim.disso
#'
#' # set `message = TRUE` to view the compliance of the regulatory guidelines.
#' calcf2(d.t, d.r, plot = FALSE)
#'
#' @export
calcf2 <- function(test, ref, path.in, file.in, path.out, file.out,
                   regulation = c("EMA", "FDA", "WHO", "Canada", "ANVISA"),
                   cv.rule = TRUE, message = FALSE, min.points = 3L,
                   f2.type = c("est.f2", "exp.f2", "bc.f2", "vc.exp.f2",
                               "vc.bc.f2", "all"), both.TR.85 = FALSE,
                   digits = 2L, time.unit = c("min", "h"),  plot = TRUE,
                   plot.start.time = 0, plot.max.unit = 24L) {
  # initial check --------------------------------------------------------------
  regulation <- match.arg(regulation)
  f2.type    <- match.arg(f2.type)
  time.unit  <- match.arg(time.unit)

  if (any(missing(test), missing(ref))) {
    stop("Both 'test' and 'ref' have to be specified.")
  }

  if (all(isTRUE(both.TR.85), regulation != "FDA")) {
    stop("'both.TR.85 = TRUE' is only applicable when 'regulation = FDA'.")
  }

  if (any(all(!missing(path.out), missing(file.out)),
          all(missing(path.out), !missing(file.out)))) {
    stop("You should provided both 'path.out' and 'file.out'.")
  }

  # read data ------------------------------------------------------------------
  if (all(missing(path.in), missing(file.in))) {
    data.t <- as.matrix(test, rownames.force = FALSE)
    data.r <- as.matrix(ref, rownames.force = FALSE)
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
    file.in <- paste0(path.in, file.in)
    if (!file.exists(file.in)) {
      stop(paste0("The file you specified does not exist. Don't forget to ",
                  "include\nthe extension 'xlsx' or 'xls' in the file name."))
    }

    sheet.names <- excel_sheets(file.in)
    if (!(test %in% sheet.names)) (
      stop("The name of the work sheet 'test' is wrong. Check your spelling.")
    )

    if (!(ref %in% sheet.names)) (
      stop("The name of the work sheet 'ref' is wrong. Check your spelling.")
    )

    # package readxl::read_excel
    data.t <- as.matrix(read_excel(file.in, test, col_types = "numeric"))
    data.r <- as.matrix(read_excel(file.in, ref, col_types = "numeric"))
  }# end read data

  # remove time 0 point if any
  data.t <- data.t[data.t[, 1] != 0, ]
  data.r <- data.r[data.r[, 1] != 0, ]

  # time points for T and R have to be same
  if (!all(data.t[, 1] == data.r[, 1])) {
    stop("Time points of two data sets should be same! ",
         "Please check your data.\n\n")
  }

  # simple way to determine if mean/individual data was used:
  # if data has > 2 columns (time, mean) ==> individual data used
  if (all(dim(data.t)[[2]] == 2, dim(data.r)[[2]] == 2)) {# for mean data
    cv.rule  <- FALSE
    use.mean <- TRUE
    mean.t   <- data.t[, 2]
    mean.r   <- data.r[, 2]
  } else {# for individual data
    use.mean <- FALSE
    nt       <- NCOL(data.t) - 1
    nr       <- NCOL(data.r) - 1
    mean.t   <- rowMeans(data.t[, 2:NCOL(data.t)], na.rm = TRUE)
    mean.r   <- rowMeans(data.r[, 2:NCOL(data.r)], na.rm = TRUE)
  }

  # where to cut time points to calculate f2 -----------------------------------
  # put time, mean T and R into same data frame
  mean.tr <- cbind(time = data.t[, 1], mean.t, mean.r)

  # determine if mean dissolution >= 85% at 15 min -----------------------------
  if (time.unit == "min") {
    if (length(which(mean.tr[, 1] <= 15)) != 0) {
      if (all(mean.tr[max(which(mean.tr[, 1] <= 15)), 2] > 85,
              mean.tr[max(which(mean.tr[, 1] <= 15)), 3] > 85)) {
        d85at15 <- "yes"
      } else d85at15 <- "no"
    } else d85at15 <- "no"
  } else {# time unit h
    if (length(which(mean.tr[, 1]*60 <= 15)) != 0) {
      if (all(mean.tr[max(which(mean.tr[, 1]*60 <= 15)), 2] > 85,
              mean.tr[max(which(mean.tr[, 1]*60 <= 15)), 3] > 85)) {
        d85at15 <- "yes"
      } else d85at15 <- "no"
    } else d85at15 <- "no"
  }


  # give all time points as initial value in case no dissolution > 85%
  tp85 <- NROW(mean.tr)

  # according to unofficial communication with FDA staff, EMA = FDA,
  # previous 'consensus' for both T and R > 85% is incorrect
  if (all(regulation != "WHO", isFALSE(both.TR.85))) {
    for (i in 1:NROW(mean.tr)) {
      if (any(mean.tr[i, 2] > 85, mean.tr[i, 3] > 85)) {
        tp85 <- i
        break
      }
    }# end EMA/FDA/Canada/ANVISA
  } else {# WHO: cut where reference reach 85% release
    for (i in 1:NROW(mean.tr)) {# 3rd column ref
      if (mean.tr[i, 3] > 85) {
        tp85 <- i
        break
      }
    } # end WHO
  }

  # to compare results from other programs with 'old convention',
  # still needs 'old consensus'.
  if (all(regulation == "FDA", isTRUE(both.TR.85))) {
    for (i in 1:NROW(mean.tr)) {
      if (all(mean.tr[i, 2] > 85, mean.tr[i, 3] > 85)) {
        tp85   <- i
        break
      }
    }
    msg1 <- paste0("*Note: Argument 'both.TR.85' is 'TRUE', which is the ",
                   "wrong interpretation\nof the guidance. This should only ",
                   "be used for cases such as checking the\ncalculation ",
                   "published in the old literature.\n\n")
  } else msg1 <- NULL

  # mean data points used for f2 calculation
  mt.tp85 <- mean.tr[1:tp85, 2]
  mr.tp85 <- mean.tr[1:tp85, 3]

  # calculate f2s if usable points >= minimum points ---------------------------
  if (tp85 >= min.points) {# whether mean or individual data, always get est.f2
    est.f2 <- 100 - 25*log10(1 + sum((mt.tp85 - mr.tp85)^2)/tp85)
  } else {
    est.f2 <- NA
  }

  # get cv/variance etc only when individual data is used ----------------------
  if (isFALSE(use.mean)) {
    # get var/cv
    var.t <- apply(data.t[, -1], 1, var, na.rm = TRUE)
    var.r <- apply(data.r[, -1], 1, var, na.rm = TRUE)

    cv.t <- sqrt(var.t)/mean.t*100
    cv.r <- sqrt(var.r)/mean.r*100

    # NaN could be generated for mean = 0, better safe than sorry.
    cv.t[is.na(cv.t)] <- 0
    cv.r[is.na(cv.r)] <- 0

    # output mean and cv later so user can see
    dtr <- cbind(mean.tr, cv.t, cv.r)

    # check cv criteria when cv.rule = TRUE -------------------------------
    if (isTRUE(cv.rule)) {
      # CV cutting point, default 10 min for EMA/WHO, FDA rule more flexible,
      # could be 15 min or even 20 min according to unofficial communication
      if (regulation %in% c("EMA", "WHO", "Canada")) {
        cv.cut.tp <- 10.0
      } else if (regulation == "FDA") {# TO DO: need better work for FDA here
        cv.cut.tp <- 15.0
      } else {# ANVISA: 1st 40% of total time points can have CV 20%.
        cv.cut.tp <- ifelse(time.unit == "min",
                            mean.tr[floor(NROW(mean.tr)*0.4), 1],
                            mean.tr[floor(NROW(mean.tr)*0.4), 1]*60)
      }

      # check if there is any time point <= cv.cut.tp
      if (time.unit == "min") {
        if (length(which(dtr[1:tp85, 1] <= cv.cut.tp)) != 0) {
          # exists point <= cv.cut.tp
          cv.tp <- max(which(dtr[1:tp85, 1] <= cv.cut.tp))
        } else {# all time points > cv.cut.tp
          cv.tp <- 0
        }
      } else {# h ==> min
        if (length(which(dtr[1:tp85, 1]*60 <= cv.cut.tp)) != 0) {
          cv.tp <- max(which(dtr[1:tp85, 1]*60 <= cv.cut.tp))
        } else {
          cv.tp <- 0
        }
      }

      # now check if CV is OK. better rounded CV first then compare
      # to 20%/10% rule. express CV precision with digits
      if (any(all(cv.tp != 0, round(dtr[1:cv.tp, 4:5], digits) <= 20,
                  round(dtr[(cv.tp + 1):tp85, 4:5], digits) <= 10),
              all(cv.tp == 0, round(dtr[1:tp85, 4:5], digits) <= 10))) {
        cv.ok <- TRUE
      } else {
        cv.ok <- FALSE
      }
    } else {# when mean profiles is used
      cv.ok <- FALSE
    }

    # determine if individual dissolution > 85% at 15 min for Canadian rule
    # determine if mean dissolution >= 85% at 15 min
    if (regulation == "Canada") {
      if (time.unit == "min") {
        if (length(which(mean.tr[, 1] <= 15)) != 0) {
          if(all(data.t[max(which(mean.tr[, 1] <= 15)), -1] > 85,
                 data.r[max(which(mean.tr[, 1] <= 15)), -1] > 85)) {
            d85at15 <- "yes"
          } else  d85at15 <- "no"
        } else  d85at15 <- "no"
      } else {# time unit h
        if (length(which(mean.tr[, 1]*60 <= 15)) != 0) {
          if (all(data.t[max(which(mean.tr[, 1]*60 <= 15)), -1] > 85,
                  data.r[max(which(mean.tr[, 1]*60 <= 15)), -1] > 85)) {
            d85at15 <- "yes"
          } else d85at15 <- "no"
        } else d85at15 <- "no"
      }
    }# end Canadian rule for 85% at 15 min

    # for exp.f2, bc.f2, vc.exp.f2, vc.bc.f2 -----------------------------------
    # next section mainly for bootstrap f2, article Shah et al, 1998
    if (f2.type != "est.f2") {# get all estimators
      # get variance
      vt.tp85 <- var.t[1:tp85]
      vr.tp85 <- var.r[1:tp85]

      # table 7 in Shah's article 1998.
      sum.mean.sqr <- sum((mt.tp85 - mr.tp85)^2)
      sum.mean.var <- sum(vt.tp85/nt + vr.tp85/nr)

      # weighted variance
      wvt <- (0.5 + vt.tp85/(vt.tp85 + vr.tp85))*vt.tp85
      wvr <- (0.5 + vr.tp85/(vt.tp85 + vr.tp85))*vr.tp85
      wvt[is.na(wvt)] <- 0
      wvr[is.na(wvr)] <- 0
      sum.mean.var.w <- sum(wvt/nt + wvr/nr)

      # now calculate other f2 estimators
      if (tp85 >= min.points) {
        exp.f2    <- 100 - 25*log10(1 + (sum.mean.sqr + sum.mean.var)/tp85)
        vc.exp.f2 <- 100 - 25*log10(1 + (sum.mean.sqr + sum.mean.var.w)/tp85)
        if (sum.mean.sqr - sum.mean.var > -tp85) {
          bc.f2 <- 100 - 25*log10(1 + (sum.mean.sqr - sum.mean.var)/tp85)
        } else bc.f2 <- NA
        if (sum.mean.sqr - sum.mean.var.w > -tp85) {
          vc.bc.f2 <- 100 - 25*log10(1 + (sum.mean.sqr - sum.mean.var.w)/tp85)
        } else vc.bc.f2 <- NA
      } else {
        exp.f2    <- NA
        bc.f2     <- NA
        vc.exp.f2 <- NA
        vc.bc.f2  <- NA
      }
    }# end f2.type != "est.f2"
  }# end isFALSE(use.mean)

  # overkill for simple f2 calculation but just in case someone need it
  if (all(!missing(path.out), !missing(file.out))) {
    if (!dir.exists(path.out)) {
      stop("The directory you specified does not exist. Check your spelling.")
    }

    path.out <- normalizePath(path.out, winslash = "/")
    path.out <- ifelse(regmatches(path.out, regexpr(".$", path.out)) == "/",
                      path.out, paste0(path.out, "/"))
    file.out <- ifelse(regmatches(file.out, regexpr(".{3}$", file.out))=="txt",
                       file.out, paste0(file.out, ".txt"))
    file.out <- paste0(path.out, file.out)

    sink(file = file.out, split = TRUE)
  }

  # print message --------------------------------------------------------------
  if (isTRUE(message)) {# set to FALSE for bootstrap f2
    # output regulation guideline
    cat("The f2 method was applied according to ", regulation, "'s BE ",
        ifelse(regulation == "FDA", "guidance", "guideline"), ".\n\n", sep = "")

    # print format depending on user set decimal width for time column
    wt  <- paste0("%-", max(nchar(mean.tr[, 1])) + 2, "s")

    # width for other column headers
    ws  <- paste0("%", digits + 8, "s")

    # width for numeric values of dissolution
    wd  <- paste0("%", digits + 8, ".", digits, "f")

    if (isFALSE(use.mean)) {## message for using individual data -----
      if (isTRUE(cv.rule)) {# cv.rule = TRUE
        cat("Individual data was provided with option 'cv.rule = TRUE',\n",
            "therefore, CV has been calculated and checked accordingly.\n\n",
            sep = "")
      } else {# cv.rule = FALSE
        cat("Individual data was provided with option 'cv.rule = FALSE',\n",
            "therefore, CV has been calculated but has not been checked.\n",
            "You should really consider setting 'cv.rule = TRUE' to comply\n",
            "regulatory requirements.\n\n", sep = "")
      }# end cv.rule = FALSE

      # output mean and CV so user can check
      cat("Calculated mean and CV as follows:\n")
      cat(sprintf(paste(wt, ws, ws, ws, ws), "Time", "Mean (T)",
                  "Mean (R)", "CV (T)", "CV (R)"), "\n")
      for (i in 1:tp85) {
        cat(sprintf(paste(wt, wd, wd, wd, wd), dtr[i, 1],
                    round(dtr[i, 2], digits), round(dtr[i, 3], digits),
                    round(dtr[i, 4], digits), round(dtr[i, 5], digits)), "\n")
      }

      # add a line to separate used data points and unused data points
      cat(paste0(rep("-", max(nchar(mean.tr[, 1])) + 2 + (digits + 8)*4 + 4),
                 collapse = ""), "\n")

      # print the rest data
      if (tp85 < NROW(dtr)) {
        for (i in (tp85 + 1):NROW(dtr)) {
          cat(sprintf(paste(wt, wd, wd, wd, wd), dtr[i, 1],
                      round(dtr[i, 2], digits), round(dtr[i, 3], digits),
                      round(dtr[i, 4], digits), round(dtr[i, 5], digits)), "\n")
        }
      }

      cat("==================================\n")
      cat("Number of units for test is      : nt = ", nt, "\n", sep = "")
      cat("Number of units for reference is : nr = ", nr, "\n\n", sep = "")

      if (all(isTRUE(cv.rule), isTRUE(cv.ok))) {
        cat("CV criteria fulfilled; therefore, f2 method can be applied.\n\n")
      } else if (all(isTRUE(cv.rule), isFALSE(cv.ok))) {# cv not ok
        if (regulation != "FDA") {
          cat("CV criteria not fulfilled; therefore, f2 method cannot",
              "be applied.\n\n")
          stop(paste0("You should consider alternative methods such as ",
                      "bootstrap f2.\n\n"))
        } else {
          warning(
            paste0("f2 was calculated while CV criterion is not strictly ",
                   "fulfilled; you \nmight want to consider an alternative ",
                   "method such as bootstrap f2.\n\n")
          )
        }
      } else if (isFALSE(cv.rule)) {
        cat("CV has not been checked.\n\n")
      }
    } else {## use mean data -----------------------------------------
      cat("Only mean data was provided, so CV cannot been checked.\n\n")

      cat(sprintf(paste(wt, ws, ws), "Time", "Mean (T)", "Mean (R)"), "\n")
      for (i in 1:tp85) {
        cat(sprintf(paste(wt, wd, wd), mean.tr[i, 1],
                    round(mean.tr[i, 2], digits),
                    round(mean.tr[i, 3], digits)), "\n")
      }

      # add a line to separate used data points and unused data points
      cat(paste0(rep("-", max(nchar(mean.tr[, 1])) + 2 + (digits + 8)*2 + 2),
                 collapse = ""), "\n")

      # print the rest data
      if (tp85 < NROW(mean.tr)) {
        for (i in (tp85 + 1):NROW(mean.tr)) {
          cat(sprintf(paste(wt, wd, wd), mean.tr[i, 1],
                      round(mean.tr[i, 2], digits),
                      round(mean.tr[i, 3], digits)), "\n")
        }
      }
      cat("\n")
    }# end use mean data

    # print f2 -----------------------------------------------------------------
    # only necessary for est.f2
    if (f2.type == "est.f2") {
      if (tp85 >= 3L) {
        cat("The time points above the dashed line are used in f2",
            "calculation.\n")
        cat("\nEstimated f2 = ", round(est.f2, digits),
            ifelse(is.null(msg1), "\n\n", "*\n\n"), sep = "")
        cat(msg1)
      } else {# tp85 < 3
        if (d85at15 == "yes") {
          cat("At least 3 time points are necessary to calculate f2, ",
              "but you have\n", tp85, ifelse(tp85 > 1, " points ", " point "),
              "only. However, ",
              ifelse(regulation != "Canada", "mean ", "individual "),
              "dissolution profiles of test and \nreference are more than ",
              "85% at 15 min. If the products are immediate-release \n",
              "formulations, then the profiles are considered similar ",
              "without statistical \nevaluation.", sep = "")
        } else {
          cat("Estimated f2 = ", round(est.f2, digits),
              ifelse(is.null(msg1), "\n\n", "*\n\n"), sep = "")
          cat(msg1)
          warning("Warning: f2 was calculated with less than 3 time points ",
                  "for information\npurpose only. You should add more ",
                  "earlier points in your dissolution method.")
        }
      }# end tp85 < 3
    } else {# not really useful. just for fun
      cat("The following f2 estimators are applicable for bootstrap method",
          "only.\n")
      cat("The time points above the dashed line are used in f2 calculation.\n")
      # exp.f2 --------------------------------------------------
      if (all(f2.type == "exp.f2", isFALSE(use.mean))) {
        cat("\nExpected f2 = ", round(exp.f2, digits),
            ifelse(is.null(msg1), "\n\n", "*\n\n"), sep = "")
        cat(msg1)
      }

      # vc.exp.f2 -----------------------------------------------
      if (all(f2.type == "vc.exp.f2", isFALSE(use.mean))) {
        cat("\nVariance-corrected expected f2 = ", round(vc.exp.f2, digits),
            ifelse(is.null(msg1), "\n\n", "*\n\n"), sep = "")
        cat(msg1)
      }

      # bc.f2 and vc.bc.f2 could be NA depending on variance ----
      if (all(f2.type == "bc.f2", isFALSE(use.mean))) {
        cat("\nBias-corrected f2 = ", round(bc.f2, digits), sep = "")
        if (sum.mean.sqr - sum.mean.var <= -tp85) {
          cat("*\nBias-corrected f2 cannot be calculated. See Shah's article\n",
              "for details (Pharm. Res., 1998, 15(6), 889-896).\n\n", sep = "")
        } else {
          cat(ifelse(is.null(msg1), "\n\n", "*\n\n"))
          cat(msg1)
        }
      }

      if (all(f2.type == "vc.bc.f2", isFALSE(use.mean))) {
        cat("\nVariance- and bias-corrected f2 = ", round(vc.bc.f2, digits),
            sep = "")
        if (sum.mean.sqr - sum.mean.var.w <= -tp85) {
          cat("*\nVariance- and bias-corrected f2 cannot be calculated. See\n",
              "Shah's article for details (Pharm. Res., 1998, 15(6), 889-896).",
              "\n\n", sep = "")
        } else {
          cat(ifelse(is.null(msg1), "\n\n", "*\n\n"))
          cat(msg1)
        }
      }

      if (all(f2.type == "all", isFALSE(use.mean))) {
        cat("\n                   Estimated f2 = ", round(est.f2, digits),
            ifelse(is.null(msg1), "\n", "*\n"), sep = "")

        cat("                    Expected f2 = ", round(exp.f2, digits),
            ifelse(is.null(msg1), "\n", "*\n"), sep = "")

        cat(" Variance-corrected expected f2 = ", round(vc.exp.f2, digits),
            ifelse(is.null(msg1), "\n", "*\n"), sep = "")

        if (is.null(msg1)) {
          cat("              Bias-corrected f2 = ", round(bc.f2, digits),
              ifelse(sum.mean.sqr - sum.mean.var <= -tp85, "*\n", "\n"),
              sep = "")
        } else {
          cat("              Bias-corrected f2 = ", round(bc.f2, digits),
              ifelse(sum.mean.sqr - sum.mean.var <= -tp85, "**\n", "*\n"),
              sep = "")
        }

        if (is.null(msg1)) {
          cat("Variance- and bias-corrected f2 = ", round(vc.bc.f2, digits),
              ifelse(sum.mean.sqr - sum.mean.var.w <= -tp85, "*\n", "\n"),
              sep = "")
        } else {
          cat("Variance- and bias-corrected f2 = ", round(vc.bc.f2, digits),
              ifelse(sum.mean.sqr - sum.mean.var.w <= -tp85, "**\n", "*\n"),
              sep = "")
        }

        if (sum.mean.sqr - sum.mean.var <= -tp85) {
          cat("---------------------------------\n")
          cat("*Bias-corrected f2 cannot be calculated. See Shah's article\n",
              "for details (Pharm. Res., 1998, 15(6), 889-896).\n\n", sep = "")
        }

        if (sum.mean.sqr - sum.mean.var.w <= -tp85) {
          cat("---------------------------------\n")
          cat("*Variance- and bias-corrected f2 cannot be calculated. See\n",
              "Shah's article for details (Pharm. Res., 1998, 15(6), 889-896).",
              "\n\n", sep = "")
        }
        cat(msg1)
      }# end f2.type = all
    }
  }# end message

  # return output back to default screen
  if (all(!missing(path.out), !missing(file.out))) {
    sink()
  }

  # plot -----------------------------------------------------------------------
  if (isTRUE(plot)) {
    # need this to get rid of "no visible binding for global variable" notes:
    time <- release <- bp.min <- qt25 <- qt75 <- bp.max <- NULL

    # add time 0 for better visual
    tmp.t <- as.data.frame(rbind(rep(0, NCOL(data.t)), data.t),
                           row.names = "", stringsAsFactors = FALSE)
    tmp.r <- as.data.frame(rbind(rep(0, NCOL(data.r)), data.r),
                           row.names = "", stringsAsFactors = FALSE)

    # RColorBrewer::brewer.pal(12, "Paired")
    # 2 "#1F78B4" and 6 "#E31A1C" for T and R
    cols <- c("R" = "#1F78B4", "T" = "#E31A1C")

    if (isFALSE(use.mean)) {# print individual profiles if suitable
      if (all(NCOL(data.t)-1 < plot.max.unit, NCOL(data.r)-1 < plot.max.unit)){
        names(tmp.t) <- c("time", paste0("unit.",1:(NCOL(tmp.t)-1)))
        names(tmp.r) <- c("time", paste0("unit.",1:(NCOL(tmp.r)-1)))
        # result data.frame: time unit release
        p.dt <- reshape(tmp.t, varying = 2:NCOL(tmp.t), direction = "long",
                        timevar = "unit", v.names = "release", idvar = "time")
        #p.dt$form <- "T"

        p.dr <- reshape(tmp.r, varying = 2:NCOL(tmp.r), direction = "long",
                        timevar = "unit", v.names = "release", idvar = "time")
        #p.dr$form <- "R"

        #p.dtr <- as.data.frame(rbind(p.dt, p.dr), stringsAsFactors = FALSE)

        p.mt <- data.frame(time = c(0, data.t[, 1]), release = c(0, mean.t),
                           stringsAsFactors = FALSE)
        p.mr <- data.frame(time = c(0, data.r[, 1]), release = c(0, mean.r),
                           stringsAsFactors = FALSE)

        #p.mtr <- as.data.frame(rbind(p.mt, p.mr), stringsAsFactors = FALSE)
        y.tick.max <- ceiling(max(c(p.dt$release, p.dr$release)))

        p.tr <- ggplot(p.dt, aes(x = time, y = release)) +
          geom_hline(yintercept = 85, size = 1.2, color = "gray",
                     alpha = 0.7, linetype = "dashed") +
          geom_line(aes(group = unit, color = "T"), alpha = 0.5, size = 0.3)+
          geom_line(data = p.dr, aes(group = unit, color = "R"),
                    alpha = 0.5, size = 0.3)+
          geom_point(data = p.mt, size = 1.5, aes(color = "T")) +
          geom_line(data = p.mt, size = 1.1, aes(color = "T")) +
          geom_point(data = p.mr, size = 1.5, aes(color = "R")) +
          geom_line(data = p.mr, size = 1.1, aes(color = "R")) +
          theme_bw()+
          theme(legend.title = element_blank(), legend.justification = c(1, 0),
                legend.position = c(0.98, 0.02)) +
          scale_colour_manual(values = cols) +
          scale_y_continuous(breaks = seq(0, y.tick.max + 5, 5),
                             limits = c(0, y.tick.max + 5)) +
          scale_x_continuous(breaks = p.mt$time) +
          ylab("Dissolution (%)") +
          xlab(paste0("Time (", time.unit, ")")) +
          ggtitle(paste0("Mean Dissolution Profiles of Test and Reference ",
                         "with Estimated f2 = ", round(est.f2, digits)),
                  subtitle = "Thin lines represent individual profiles")
      } else {# use boxplot
        p.dt <- data.frame(
          time   = tmp.t[, 1],
          mean   = rowMeans(tmp.t[, -1], na.rm = TRUE),
          median = apply(tmp.t[, -1], 1, median, na.rm = TRUE),
          qt25   = apply(tmp.t[, -1], 1, quantile, probs = 0.25, na.rm = TRUE),
          qt75   = apply(tmp.t[, -1], 1, quantile, probs = 0.75, na.rm = TRUE),
          bp.min = apply(tmp.t[, -1], 1, bpwhisker.l),
          bp.max = apply(tmp.t[, -1], 1, bpwhisker.u),
        # qt05   = apply(tmp.t[, -1], 1, quantile, probs = 0.05, na.rm = TRUE),
        # qt95   = apply(tmp.t[, -1], 1, quantile, probs = 0.95, na.rm = TRUE),
        #  form   = "T",
          stringsAsFactors = FALSE
        )

        p.dr <- data.frame(
          time   = tmp.r[, 1],
          mean   = rowMeans(tmp.r[, -1], na.rm = TRUE),
          median = apply(tmp.r[, -1], 1, median, na.rm = TRUE),
          qt25   = apply(tmp.r[, -1], 1, quantile, probs = 0.25, na.rm = TRUE),
          qt75   = apply(tmp.r[, -1], 1, quantile, probs = 0.75, na.rm = TRUE),
          bp.min = apply(tmp.r[, -1], 1, bpwhisker.l),
          bp.max = apply(tmp.r[, -1], 1, bpwhisker.u),
        # qt05   = apply(tmp.r[, -1], 1, quantile, probs = 0.05, na.rm = TRUE),
        # qt95   = apply(tmp.r[, -1], 1, quantile, probs = 0.95, na.rm = TRUE),
        #  form   = "R",
          stringsAsFactors = FALSE
        )

        #p.dtr <- as.data.frame(rbind(p.dt, p.dr), stringsAsFactors = FALSE)
        y.tick.max <- ceiling(max(c(p.dt$bp.max, p.dr$bp.max)))

        p.tr <- ggplot(p.dt, aes(x = time, y = mean)) +
          geom_hline(yintercept = 85, size = 1.2, color = "gray",
                     alpha = 0.7, linetype = "dashed") +
          geom_point(size = 1.5, aes(color = "T")) +
          geom_line(size = 1.1, aes(color = "T")) +
          geom_boxplot(aes(ymin = bp.min, lower = qt25, middle = median,
                           upper = qt75, ymax = bp.max, group = time),
                           fill = cols[[2]], color = cols[[2]],
                       alpha = 0.5, stat = "identity", width = 0.7) +
          geom_point(data = p.dr, size = 1.5, aes(color = "R")) +
          geom_line(data = p.dr, size = 1.1, aes(color = "R")) +
          geom_boxplot(aes(ymin = bp.min, lower = qt25, middle = median,
                           upper = qt75, ymax = bp.max, group = time),
                           fill = cols[[1]], color = cols[[1]],
                       data = p.dr, alpha = 0.5, stat = "identity",
                       width = 0.7) +
          theme_bw()+
          theme(legend.title = element_blank(), legend.justification = c(1, 0),
                legend.position = c(0.98, 0.02)) +
          scale_color_manual(values = cols) +
          scale_y_continuous(breaks = seq(0, y.tick.max + 5, 5),
                             limits = c(0, y.tick.max + 5)) +
          scale_x_continuous(breaks = p.dt$time) +
          ylab("Dissolution (%)") +
          xlab(paste0("Time (", time.unit, ")")) +
          ggtitle(paste0("Mean Dissolution Profiles of Test and Reference ",
                         "with Estimated f2 = ", round(est.f2, digits)),
                  subtitle = "Individual prifiles are summarized by boxplots")
      }
    } else {# with mean data
      names(tmp.t)[1:2] <- names(tmp.r)[1:2] <- c("time", "mean")

      y.tick.max <- ceiling(max(tmp.t$mean, tmp.r$mean))

      p.tr <- ggplot(tmp.t, aes(x = time, y = mean)) +
        geom_hline(yintercept = 85, size = 1.2, color = "gray",
                   alpha = 0.7, linetype = "dashed") +
        geom_point(size = 1.5, aes(color = "T")) +
        geom_line(size = 1.1, aes(color = "T")) +
        geom_point(data = tmp.r, size = 1.5, aes(color = "R")) +
        geom_line(data = tmp.r, size = 1.1, aes(color = "R")) +
        theme_bw()+
        theme(legend.title = element_blank(), legend.justification = c(1, 0),
              legend.position = c(0.98, 0.02)) +
        scale_color_manual(values = cols) +
        scale_y_continuous(breaks = seq(0, y.tick.max + 5, 5),
                           limits = c(0, y.tick.max + 5)) +
        scale_x_continuous(breaks = tmp.t$time) +
        ylab("Dissolution (%)") +
        xlab(paste0("Time (", time.unit, ")")) +
        ggtitle(paste0("Mean Dissolution Profiles of Test and Reference ",
                       "with Estimated f2 = ", round(est.f2, digits)))
    }
    print(p.tr)
  }# end plot

  # return non-rounded result for further calculation if necessary--------------
  # same for bootstrap or non-bootstrap
  if (f2.type == "est.f2") {
    # invisible(c(est.f2 = est.f2, tp = tp85))
    invisible(
      data.frame(f2.type = "est.f2", f2.value = est.f2, f2.tp = tp85,
                 d85at15 = d85at15, regulation = regulation,
                 cv.rule = cv.rule, min.points = min.points,
                 stringsAsFactors = FALSE)
      )
  } else if (all(f2.type == "exp.f2", isFALSE(use.mean))) {
    # invisible(c(exp.f2 = exp.f2, tp = tp85))
    invisible(
      data.frame(f2.type = "exp.f2", f2.value = exp.f2, f2.tp = tp85,
                 d85at15 = d85at15, regulation = regulation,
                 cv.rule = cv.rule, min.points = min.points,
                 stringsAsFactors = FALSE)
    )
  } else if (all(f2.type == "bc.f2", isFALSE(use.mean))) {
    # invisible(c(bc.f2 = bc.f2, tp = tp85))
    invisible(
      data.frame(f2.type = "bc.f2", f2.value = bc.f2, f2.tp = tp85,
                 d85at15 = d85at15, regulation = regulation,
                 cv.rule = cv.rule, min.points = min.points,
                 stringsAsFactors = FALSE)
    )
  } else if (all(f2.type == "vc.exp.f2", isFALSE(use.mean))) {
    # invisible(c(vc.exp.f2 = vc.exp.f2, tp = tp85))
    invisible(
      data.frame(f2.type = "vc.exp.f2", f2.value = vc.exp.f2, f2.tp = tp85,
                 d85at15 = d85at15, regulation = regulation,
                 cv.rule = cv.rule, min.points = min.points,
                 stringsAsFactors = FALSE)
    )
  } else if (all(f2.type == "vc.bc.f2", isFALSE(use.mean))) {
    # invisible(c(vc.bc.f2 = vc.bc.f2, tp = tp85))
    invisible(
      data.frame(f2.type = "vc.bc.f2", f2.value = vc.bc.f2, f2.tp = tp85,
                 d85at15 = d85at15, regulation = regulation,
                 cv.rule = cv.rule, min.points = min.points,
                 stringsAsFactors = FALSE)
    )
  } else if (all(f2.type == "all", isFALSE(use.mean))) {
    # invisible(c(est.f2 = est.f2, exp.f2 = exp.f2, bc.f2 = bc.f2,
    #             vc.exp.f2 = vc.exp.f2, vc.bc.f2 = vc.bc.f2, tp = tp85))
    invisible(
      data.frame(
        f2.type = c("est.f2", "exp.f2", "bc.f2", "vc.exp.f2", "vc.bc.f2"),
        f2.value = c(est.f2, exp.f2, bc.f2, vc.exp.f2, vc.bc.f2),
        f2.tp = tp85, d85at15 = d85at15, regulation = regulation,
        cv.rule = cv.rule, min.points = min.points, stringsAsFactors = FALSE)
    )
  }
}
