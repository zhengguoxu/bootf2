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
## - rpt.ci
## - rpt.f2
## - rpt.concise
## - rpt.detailed
## - rpt.intermediate
## - rpt.screen

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

#' @param boot.info A data frame of bootstrap information from `bootf2`
#'   function.
#' @rdname helper
ci.header <- function(boot.info) {
  # need some trick for print format depending on user set decimal
  # for CI type: "%24s", rest "%3s", lower, "%3s", upper
  # header ----
  if (boot.info$digits <= 2) {
    cat("          Types of         Lower   Upper\n",
        "    Confidence Intervals   <----------->\n", sep = "")
  } else {
    cat("          Types of         Lower", rep(" ", 2*boot.info$digits - 1),
        "Upper\n", "    Confidence Intervals   <----",
        rep("-", 2*boot.info$digits - 1), "---->\n", sep = "")
  }
}

#' @param tp,ref.dp Numeric vector of time points \code{tp} and their
#'   corresponding mean dissolution profiles \code{ref.dp}.
#' @param digits An integer indicating the decimal points for the output.
#' @param model Strings of model names. Currently only 'Weibull' and
#'   'first-order' models are supported.
#' @param max.disso Numeric value indicating the maximum dissolution.
#' @param time.unit Character strings indicating the unit of time. It should
#'   be either \code{"min"} for minute or \code{"h"} for hour. It is mainly
#'   used for checking CV rules and making plot. @seealso [calcf2()].
#' @rdname helper
mod.ref <- function(tp, ref.dp, digits = 4, model, max.disso, time.unit) {
  # add random noise in case mean values comes from simulation
  ref.dp <- ref.dp + rnorm(length(ref.dp), sd = 0.01)

  fmax.ini <- max(ref.dp) + 1
  if (fmax.ini > max.disso) fmax.ini <- max.disso

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

    k.ini <- -coef(mod.ini)[[1]]

    ctr <- nls.control(maxiter = 500)

    mod <- minpack.lm::nlsLM(
      ref.dp ~ fmax*(1 - exp(-k*(tp - tlag))), data = tmp1, control = ctr,
      start = c(fmax = fmax.ini, tlag = tlag.ini, k = k.ini)
    )

    # coef(mod)[[2]] --> tlag,
    if (any(all(time.unit == "min", coef(mod)[[2]] <= 5),
            coef(mod)[[2]] <= 0)) {
      mod.par <- data.frame(model = "first-order",
                            fmax = coef(mod)[[1]], fmax.cv = 5,
                            tlag = 0, tlag.cv = 0,
                            k = coef(mod)[[3]], k.cv = 30,
                            time.unit = time.unit, stringsAsFactors = FALSE)
    } else {
      mod.par <- data.frame(model = "first-order",
                            fmax  = coef(mod)[[1]], fmax.cv = 5,
                            tlag = coef(mod)[[2]], tlag.cv = 30,
                            k = coef(mod)[[3]], k.cv = 30,
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


#' @param f2type Character strings indicating the f2 type.
#' @param btsum A data frame of descriptive statistics of the bootstrap data
#'   set.
#' @rdname helper
rpt.ci <- function(f2type, btsum, boot.info) {
  cat("------------------------------------------------------------\n")
  cat("  Out of", boot.info$n.boots, "bootstrapped data sets,\n")
  cat("  - Number of f2 calculated with 1 time point             :   ",
      btsum[btsum$f2.type == f2type, ]$f2.tp1, "\n", sep = "")
  cat("  - Number of f2 calculated with 2 time point             :   ",
      btsum[btsum$f2.type == f2type, ]$f2.tp2, "\n", sep = "")
  cat("  - Number of f2 cannot be calculated (NA)                :   ",
      btsum[btsum$f2.type == f2type, ]$boot.na, "\n", sep = "")
  cat("  - Number of T and R dissolving more than 85% at 15 min  :   ",
      btsum[btsum$f2.type == f2type, ]$d85at15, "\n", sep = "")
  cat("------------------------------------------------------------\n")
  if (boot.info$f2.type == "all") {
    cat("_____________________________________________________",
        "_________________\n\n", sep = "")
  } else {
    cat("\n\n")
  }
}

#' @param f2type Character strings indicating the f2 type.
#' @param f2o Vector of f2 values calculated with the original data set
#' @param a.jack Data frame of acceleration from `jackf2` function
#' @rdname helper
rpt.f2 <- function(f2type, f2o, boot.info, a.jack, btsum) {
  cat("  - with original data                        :   ",
      round(f2o[[f2type]], boot.info$digits), "\n", sep = "")
  if (boot.info$ci.type %in% c("all", "bca.jackknife")) {
    for (i in seq_len(NROW(a.jack))) {
      cat("  - with original data (by Jackknife, ", a.jack[i, "type"],
          ")  :   ", round(a.jack[i, f2type], boot.info$digits), "\n", sep = "")
    }
  }
  cat("  - with bootstrapped data (Mean)             :   ",
      round(btsum[btsum$f2.type == f2type, ]$boot.mean, boot.info$digits),
      "\n", sep = "")
  cat("  - with bootstrapped data (Median)           :   ",
      round(btsum[btsum$f2.type == f2type, ]$boot.median, boot.info$digits),
      "\n\n", sep = "")

  cat("-----------------------\n")
  cat("* Confidence Intervals\n")
  cat("-----------------------\n")

  ci.header(boot.info)
}

#' @param boot.f2.ci Matrix of f2 values from bootstrap data sets
#' @rdname helper
rpt.concise <- function(boot.f2.ci, boot.info, f2o, a.jack, btsum) {
  rpt.screen(boot.f2.ci, boot.info, f2o, a.jack, btsum)
  sysinfo <- Sys.info()
  cat("The current report was generated at ",
      format(Sys.time(), "%H:%M:%S"), " on ", format(Sys.time(), "%F"),
      " ", format(Sys.time(), "%Z"), " by\nuser '", sysinfo["user"],
      "' on machine '", sysinfo["nodename"], "', and saved as text file\n'",
      boot.info$file.out, "' at the location:\n'", boot.info$path.out, "'.\n",
      sep = "")
  cat("=============================================================",
      "=========\n\n", sep = "")
  # end concise report style.
}

#' @param data.t,data.r Input data sets for test and reference.
#' @param boot.t,boot.r List of bootstrap data sets for test and reference.
#' @param boot.f2 Matrix of f2 calculated from bootstrap data sets.
#' @rdname helper
rpt.detailed <- function(data.t, data.r, boot.t, boot.r, boot.f2, boot.info,
                         f2o) {
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

  tpo <- data.r[, 1][data.r[, 1] > 0][1:f2o[["f2.tp"]]]
  cat("\n-----------------------------------\n\n")
  if (boot.info$f2.type %in% c("all", "est.f2")) {
    cat("- Estimated f2                            :   ", f2o[["est.f2"]],
        "\n", sep = "")
  }

  if (boot.info$f2.type %in% c("all", "exp.f2")) {
    cat("- Expected f2                             :   ", f2o[["exp.f2"]],
        "\n", sep = "")
  }

  if (boot.info$f2.type %in% c("all", "bc.f2")) {
    cat("- Bias-corrected f2                       :   ", f2o[["bc.f2"]],
        "\n", sep = "")
  }

  if (boot.info$f2.type %in% c("all", "vc.exp.f2")) {
    cat("- Variance-corrected expected f2          :   ", f2o[["vc.exp.f2"]],
        "\n", sep = "")
  }

  if (boot.info$f2.type %in% c("all", "vc.bc.f2")) {
    cat("- Variance- and bias-corrected f2         :   ", f2o[["vc.bc.f2"]],
        "\n", sep = "")
  }

  cat("- Time Points Used                        :   ", f2o[["f2.tp"]], " (",
      paste0(tpo[1:length(tpo) - 1], rep(", ", length(tpo) - 1)), "and ",
      tpo[length(tpo)], " ", boot.info$time.unit, ")\n", sep = "")

  cat("- T & R dissolved more than 85% at 15 min :   ",
      ifelse(unique(f2o[["d85at15"]]) == 1, "Yes", "no"), "\n", sep = "")
  cat("-------------------------------------------------------",
      "-------------------------\n\n\n", sep = "")

  for (i in 1:boot.info$n.boots) {
    dimnames(boot.t[[i]])[[2]] <-
      c("time", paste0("U", 1:(NCOL(boot.t[[i]]) - 1)))

    dimnames(boot.r[[i]])[[2]] <-
      c("time", paste0("U", 1:(NCOL(boot.r[[i]]) - 1)))

    tpb <- boot.t[[i]][, 1][boot.t[[i]][, 1] > 0][1:boot.f2[i, "f2.tp"]]

    cat("Bootstrap data set for test       :       ", i, "\n", sep = "")
    cat("-----------------------------------\n")
    print(boot.t[[i]], row.names = FALSE)
    cat("\n")
    cat("Bootstrap data set for reference  :       ", i, "\n", sep = "")
    cat("-----------------------------------\n")
    print(boot.r[[i]], row.names = FALSE)
    cat("\n-----------------------------------\n\n")
    if (boot.info$f2.type %in% c("all", "est.f2")) {
      cat("- Estimated f2                            :   ",
          boot.f2[[i, "est.f2"]], "\n", sep = "")
    }

    if (boot.info$f2.type %in% c("all", "exp.f2")) {
      cat("- Expected f2                             :   ",
          boot.f2[[i, "exp.f2"]], "\n", sep = "")
    }

    if (boot.info$f2.type %in% c("all", "bc.f2")) {
      cat("- Bias-corrected f2                       :   ",
          boot.f2[[i, "bc.f2"]], "\n", sep = "")
    }

    if (boot.info$f2.type %in% c("all", "vc.exp.f2")) {
      cat("- Variance-corrected expected f2          :   ",
          boot.f2[[i, "vc.exp.f2"]], "\n", sep = "")
    }

    if (boot.info$f2.type %in% c("all", "vc.bc.f2")) {
      cat("- Variance- and bias-corrected f2         :   ",
          boot.f2[[i, "vc.bc.f2"]], "\n", sep = "")
    }

    cat("- Time Points Used                        :   ",
        boot.f2[[i, "f2.tp"]], " (",
        paste0(tpb[1:(length(tpb) - 1)], rep(", ", length(tpo) - 1)), "and ",
        tpb[length(tpb)], " ", boot.info$time.unit, ")\n", sep = "")
    cat("- T & R dissolved more than 85% at 15 min :   ",
        ifelse(boot.f2[i, "d85at15"] == 1, "Yes", "no"), "\n", sep = "")
    cat("-------------------------------------------------------",
        "-------------------------\n\n", sep = "")
  } # end loop for i data set
  cat("\n=====================================================",
      "===========================\n\n", sep = "")
}


#' @rdname helper
rpt.info <- function(boot.info) {
  # system and program info ----------------------------------------
  # output function settings so anyone who read the results can reproduce it
  cat("============================================\n")
  cat("| Function Settings and System Information |\n")
  cat("============================================\n\n")
  cat("---------------------\n")
  cat("* Function Settings\n")
  cat("---------------------\n")
  cat("  - test              :   ", boot.info$test, "\n", sep = "")
  cat("  - ref               :   ", boot.info$ref, "\n", sep = "")
  cat("  - n.boots           :   ", boot.info$n.boots, "\n", sep = "")
  cat("  - seed              :   ", boot.info$seed, "\n", sep = "")
  cat("  - digits            :   ", boot.info$digits, "\n", sep = "")
  cat("  - alpha             :   ", boot.info$alpha, " (",
      (1 - 2*boot.info$alpha)*100, "% CI)\n", sep = "")
  cat("  - regulation        :   ", boot.info$regulation, "\n", sep = "")
  cat("  - min.points        :   ", boot.info$min.points, "\n", sep = "")
  cat("  - both.TR.85        :   ", boot.info$both.TR.85, "\n", sep = "")
  cat("  - print.report      :   ", boot.info$print.report, "\n", sep = "")
  cat("  - report.style      :   ", boot.info$report.style, "\n", sep = "")
  cat("  - f2.type           :   ", boot.info$f2.type, "\n", sep = "")
  cat("  - ci.type           :   ", boot.info$ci.type, "\n", sep = "")
  cat("  - quantile.type     :   ", boot.info$quantile.type, "\n", sep = "")
  cat("  - jackknife.type    :   ", boot.info$jackknife.type, "\n", sep = "")
  cat("  - time.unit         :   ", boot.info$time.unit, "\n", sep = "")
  cat("  - output.to.screen  :   ", boot.info$output.to.screen, "\n", sep = "")
  cat("  - sim.data.out      :   ", boot.info$sim.data.out, "\n", sep = "")
  cat("  - path.in           :   ", boot.info$path.in, "\n", sep = "")
  cat("  - file.in           :   ", boot.info$file.in, "\n", sep = "")
  cat("  - path.out          :   ", boot.info$path.out, "\n", sep = "")
  cat("  - file.out          :   ", boot.info$file.out, "\n\n", sep = "")

  cat("---------------------\n")
  cat("* System Information\n")
  cat("---------------------\n")
  sysinfo0 <- sessionInfo()
  rinfo   <- strsplit(sessionInfo()$R.version$version.string, " ")
  sysinfo1 <- Sys.info()
  bootf2v <- as.character(packageVersion("bootf2"))
  cat("  - OS Platform             :   ", sysinfo0$platform, "\n", sep = "")
  cat("  - OS Name and Release     :   ", sysinfo0$running, "\n", sep = "")
  cat("  - Machine Node Name       :   ",
      sysinfo1["nodename"], "\n", sep = "")
  cat("  - User Name               :   ",
      sysinfo1["user"], "\n", sep = "")
  cat("  - Time Zone               :   ",
      Sys.timezone(), "\n", sep = "")
  cat("  - R Version               :   ",
      rinfo[[1]][3], " ", rinfo[[1]][4], "\n", sep = "")
  cat("  - Package bootf2 Version  :   ", bootf2v, "\n", sep = "")
  cat("_____________________________________________________",
      "_________________\n\n", sep = "")
}# end function rpt.info


#' @rdname helper
rpt.intermediate <- function(boot.info, boot.f2) {
  cat("Individual f2 for each bootstrapped data set\n")
  cat("___________________________________________________________",
      "___________\n\n", sep = "")
  if (boot.info$f2.type == "all") {
    if (boot.info$digits <= 2) {
      cat("Bootstrap  est.f2  exp.f2  bc.f2   vc.exp.f2  vc.bc.f2  ",
          "f2.tp    TR>85%\n", sep = "")
      for (i in 1:boot.info$n.boots) {
        cat(sprintf(paste0("%2s", "%-7s", "%2s",
                           paste0("%-6", ".", boot.info$digits, "f"), "%2s",
                           paste0("%-6", ".", boot.info$digits, "f"), "%2s",
                           paste0("%-6", ".", boot.info$digits, "f"), "%2s",
                           paste0("%-9", ".", boot.info$digits, "f"), "%2s",
                           paste0("%-8", ".", boot.info$digits, "f"), "%2s",
                           "%-9s", "%-9s"),
                    "", i, "", boot.f2[i, 1], "", boot.f2[i, 2], "",
                    boot.f2[i, 3], "", boot.f2[i, 4], "", boot.f2[i, 5],
                    "", boot.f2[i, 6], boot.f2[i, 7]),
            "\n", sep = "")
      }
    } else { # digits > 2
      cat("Bootstrap  est.f2", rep(" ", boot.info$digits - 1), "exp.f2",
          rep(" ", boot.info$digits - 1), "bc.f2 ",
          rep(" ", boot.info$digits - 1), "vc.exp.f2",
          rep(" ", boot.info$digits - 1), "vc.bc.f2",
          rep(" ", boot.info$digits - 1), "f2.tp    TR>85%\n", sep = "")
      for (i in 1:boot.info$n.boots) {
        cat(sprintf(paste0("%2s", "%-7s", "%2s",
                           paste0("%-", boot.info$digits + 3, ".",
                                  boot.info$digits, "f"), "%2s",
                           paste0("%-", boot.info$digits + 3, ".",
                                  boot.info$digits, "f"), "%2s",
                           paste0("%-", boot.info$digits + 3, ".",
                                  boot.info$digits, "f"), "%2s",
                           paste0("%-", boot.info$digits + 6, ".",
                                  boot.info$digits, "f"), "%2s",
                           paste0("%-", boot.info$digits+5, ".",
                                  boot.info$digits, "f"), "%2s",
                           "%-9s", "%-9s"), "", i, "",
                    boot.f2[i, 1], "", boot.f2[i, 2], "", boot.f2[i, 3], "",
                    boot.f2[i, 4], "", boot.f2[i, 5], "", boot.f2[i, 6],
                    boot.f2[i, 7]), "\n", sep = "")
      }
    } # end digits > 2
  } else { # f2.type individual
    cat("Bootstrap No.   ", boot.info$f2.type,
        rep(" ", boot.info$digits + 13 - nchar(boot.info$f2.type)),
        "f2.tp    TR>85%\n", sep = "")
    for (i in 1:boot.info$n.boots) {
      cat(sprintf(paste0("%4s", "%-9s", "%3s",
                         paste0("%-", boot.info$digits + 13, ".",
                                boot.info$digits, "f"), "%-9s", "%-9s"),
                  "", i, "", boot.f2[i, 1], boot.f2[i, 2], boot.f2[i, 3]),
          "\n", sep = "")
    }
  } # end f2.type
  cat("___________________________________________________________",
      "___________\n\n", sep = "")
}


#' @rdname helper
rpt.screen <- function(boot.f2.ci, boot.info, f2o, a.jack, btsum) {

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

  cat("=================================================================\n")
  cat("* T & R dissolved more than 85% at 15 min in original data : ",
      ifelse(f2o[["d85at15"]] == 1, "YES|", "NO |"), "\n", sep = "")
  cat("=================================================================\n\n")

  # f2 and confidence intervals ------------------------------------
  cat("============================================\n")
  cat("|              Main Results:               |\n")
  cat("|  f2 and Its Confidence Intervals (CI)s   |\n")
  cat("============================================\n\n")

  # width for numeric values of ci
  cilo <- if (boot.info$digits <= 2) {
    paste0("%-5", ".", boot.info$digits, "f")
  } else {
    paste0("%-", boot.info$digits + 3, ".", boot.info$digits, "f")
  }

  ciup <- if (boot.info$digits <= 2) {
    paste0("%5", ".", boot.info$digits, "f")
  } else {
    paste0("%", boot.info$digits + 3, ".", boot.info$digits, "f")
  }

  # for estimated f2 ------------------------------------------
  if (boot.info$f2.type %in% c("all", "est.f2")) {
    ci.estf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "est.f2", ]
    ci.estf2.n <- NROW(ci.estf2)
    cat("-------------------------------\n")
    cat("* Estimated f2 Values (est.f2)\n")
    cat("-------------------------------\n")

    rpt.f2("est.f2", f2o, boot.info, a.jack, btsum)

    for (i in 1:ci.estf2.n) {
      cat(sprintf(paste0("%24s", "%3s", cilo, "%3s", ciup),
                  ci.estf2[i, 2], "", ci.estf2[i, 3], "", ci.estf2[i, 4]),
          "\n", sep = "")
    }

    rpt.ci("est.f2", btsum, boot.info)
  } # end est.f2

  # for expected f2 -------------------------------------------
  if (boot.info$f2.type %in% c("all", "exp.f2")) {
    ci.expf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "exp.f2", ]
    ci.expf2.n <- NROW(ci.expf2)
    cat("------------------------------\n")
    cat("* Expected f2 Values (exp.f2)\n")
    cat("------------------------------\n")

    rpt.f2("exp.f2",f2o, boot.info, a.jack, btsum)

    for (i in 1:ci.expf2.n) {
      cat(sprintf(paste0("%24s", "%3s", cilo, "%3s", ciup),
                  ci.expf2[i, 2], "", ci.expf2[i, 3], "", ci.expf2[i, 4]),
          "\n", sep = "")
    }

    rpt.ci("exp.f2", btsum, boot.info)
  } # end exp.f2

  # for bias-corrected f2 -------------------------------------
  if (boot.info$f2.type %in% c("all", "bc.f2")) {
    ci.bcf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "bc.f2", ]
    ci.bcf2.n <- NROW(ci.bcf2)
    cat("-----------------------------------\n")
    cat("* Bias-Corrected f2 Values (bc.f2)\n")
    cat("-----------------------------------\n")

    rpt.f2("bc.f2", f2o, boot.info, a.jack, btsum)

    for (i in 1:ci.bcf2.n) {
      cat(sprintf(paste0("%24s", "%3s", cilo, "%3s", ciup),
                  ci.bcf2[i, 2], "", ci.bcf2[i, 3], "", ci.bcf2[i, 4]),
          "\n", sep = "")
    }

    rpt.ci("bc.f2", btsum, boot.info)
  } # end bc.f2

  # for variance corrected expected f2 ------------------------
  if (boot.info$f2.type %in% c("all", "vc.exp.f2")) {
    ci.vcexpf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "vc.exp.f2", ]
    ci.vcexpf2.n <- NROW(ci.vcexpf2)
    cat("----------------------------------------------------\n")
    cat("* Variance-corrected Expected f2 Values (vc.exp.f2)\n")
    cat("----------------------------------------------------\n")

    rpt.f2("vc.exp.f2", f2o, boot.info, a.jack, btsum)

    for (i in 1:ci.vcexpf2.n) {
      cat(sprintf(paste0("%24s", "%3s", cilo, "%3s", ciup), ci.vcexpf2[i, 2],
                  "", ci.vcexpf2[i, 3], "", ci.vcexpf2[i, 4]),
          "\n", sep = "")
    }

    rpt.ci("vc.exp.f2", btsum, boot.info)
  } # end vc.exp.f2

  # for variance- and bias-corrected f2 -----------------------
  if (boot.info$f2.type %in% c("all", "vc.bc.f2")) {
    ci.vcbcf2   <- boot.f2.ci[boot.f2.ci[["f2.type"]] == "vc.bc.f2", ]
    ci.vcbcf2.n <- NROW(ci.vcbcf2)
    cat("----------------------------------------------------\n")
    cat("* Variance- and Bias-Corrected f2 Values (vc.bc.f2)\n")
    cat("----------------------------------------------------\n")

    rpt.f2("vc.bc.f2", f2o, boot.info, a.jack, btsum)

    for (i in 1:ci.vcbcf2.n) {
      cat(sprintf(paste0("%24s", "%3s", cilo, "%3s", ciup),
                  ci.vcbcf2[i, 2], "", ci.vcbcf2[i, 3], "", ci.vcbcf2[i, 4]),
          "\n", sep = "")
    }

    rpt.ci("vc.bc.f2", btsum, boot.info)
  } # end vc.bc.f2

  rpt.info(boot.info)
}# end rpt.screen function
