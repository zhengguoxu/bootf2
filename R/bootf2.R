#' Estimate 90% Confidence Intervals of \eqn{f_2}{f2} with Bootstrap Methodology
#'
#' Main function to estimate 90% confidence intervals of \eqn{f_2}{f2} using
#' bootstrap methodology.
#'
#' @importFrom stats qnorm pnorm glm
#' @importFrom utils packageVersion sessionInfo
#' @usage
#' bootf2(test, ref, path.in, file.in, path.out, file.out,
#'        n.boots = 10000L, seed = 306L, digits = 2L, alpha = 0.05,
#'        regulation = c("EMA", "FDA", "WHO","Canada", "ANVISA"),
#'        min.points = 1L, both.TR.85 = FALSE, print.report = TRUE,
#'        report.style = c("concise", "intermediate", "detailed"),
#'        f2.type = c("all", "est.f2", "exp.f2", "bc.f2",
#'                    "vc.exp.f2", "vc.bc.f2"),
#'        ci.type = c("all", "normal", "basic", "percentile",
#'                    "bca.jackknife", "bca.boot"),
#'        quantile.type = c("all", as.character(1:9), "boot"),
#'        jackknife.type = c("all", "nt+nr", "nt*nr", "nt=nr"),
#'        time.unit = c("min", "h"), output.to.screen = TRUE,
#'        sim.data.out = FALSE)
#'
#' @param test,ref *Data frames* of dissolution profiles of test and reference
#'   product if `path.in` and `file.in` are not specified; otherwise, they
#'   should be *character* strings indicating the worksheet names of the Excel
#'   file where the dissolution data is saved. See Input/Output in Details.
#' @param path.in,file.in,path.out,file.out *Character* strings of input and
#'   output directories and file names. See Input/Output in Details.
#' @param n.boots An *integer* indicating the number of bootstrap samples.
#' @param seed *Integer* seed value for reproducibility. If missing, a random
#'   seed will be generated for reproducibility purpose.
#' @param digits An *integer* indicating the decimal points for the output.
#' @param alpha A *numeric* value between 0 and 1 to estimate
#'   \eqn{(1-2\times \alpha)\times 100}{(1 - 2*alpha)*100} confidence interval.
#' @param regulation *Character* strings indicating regulatory guidelines.
#'   @seealso [calcf2()] for details on regulation rules.
#' @param min.points An *integer* indicating the minimum time points to be used
#'   to calculate \eqn{f_2}{f2}. For conventional \eqn{f_2}{f2} calculation, the
#'   default is 3, however, for bootstrap \eqn{f_2}{f2}, the value should be
#'   lower as there might be less time points available in certain bootstrap
#'   samples. The default is 1. @seealso [calcf2()].
#' @param both.TR.85 *Logical*. If `TRUE`, and if `regulation = "FDA"`, all
#'   measurements up to the time points at which both test and reference
#'   products dissolve more than 85% will be used to calculate \eqn{f_2}{f2}.
#'   This is the conventional, but incorrect, interpretation of the US FDA rule.
#'   Therefore, the argument should only be set to `TRUE` for validation purpose
#'   such as comparing the results from old literatures that use the wrong
#'   interpretation to calculate \eqn{f_2}{f2}. @seealso [calcf2()] for details
#'   on regulation rules.
#' @param print.report *Logical*. If `TRUE`, a plain text report will be
#'   produced. See Input/Output in Details.
#' @param report.style `"concise"` style produces the estimators and their
#'   confidence intervals; `"intermediate"` style adds a list of individual
#'   \eqn{f_2}{f2}s for all bootstrap samples in the end of `"concise"` report;
#'   `"detailed"` style further adds individual bootstrap samples along with
#'   their \eqn{f_2}{f2}s in the end of `"intermediate"` report. See
#'   Input/Output in Details.
#' @param f2.type *Character* strings indicating which type of \eqn{f_2}{f2}
#'   estimator should be calculated. See Types of estimators in Details.
#' @param ci.type *Character* strings indicating which type of confidence
#'   interval should be estimated. See Types of confidence intervals in
#'   Details.
#' @param quantile.type *Character* strings indicating the type of percentile.
#' @param jackknife.type *Character* strings indicating the type of jackknife
#'   method. See Details.
#' @param time.unit *Character* strings indicating the unit of time. It should
#'   be either `"min"` for minute or `"h"` for hour. It is mainly used for
#'   checking CV rules and making plot. @seealso [calcf2()].
#' @param output.to.screen *Logical*. If `TRUE`, a `"concise"` style summary
#'   report will be printed on screen. See Input/Output in Details.
#' @param sim.data.out *Logical*. If `TRUE`, all individual bootstrap data
#'   sets will be included in the output.
#'
#' @returns A list of 3 or 5 components.
#' - `boot.ci`: A *data frame* of bootstrap \eqn{f_2}{f2} confidence intervals.
#' - `boot.f2`: A *data frame* of all individual \eqn{f_2}{f2} values for all
#'   bootstrap data set.
#' - `boot.info`: A *data frame* with detailed information of bootstrap for
#'   reproducibility purpose, such as all arguments used in the function, time
#'   points used for calculation of \eqn{f_2}{f2}, and the number of `NA`s.
#' - `boot.summary`: A *data frame* with descriptive statistics of the
#'    bootstrap \eqn{f_2}{f2}.
#' - `boot.t` and `boot.r`: *Lists* of individual bootstrap samples for test
#'   and reference product if `sim.data.out = TRUE`.
#'
#' @details
#' ## Minimum required arguments that must be provided by the user
#' Arguments `test` and `ref` must be provided by the user. They should be `R`
#' `data frames`, with *time as the first column*, and all individual profiles
#' profiles as the rest columns. The actual names of the columns do not matter
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
#' are individual dissolution profiles. The first row should be column names,
#' such as time, unit01, unit02, ... The actual names of the columns do not
#' matter as they will be renamed internally.
#'
#' Arguments `path.out` and `file.out` are the names of the output directory
#' and file. If they are not provided, but argument `print.report` is `TRUE`,
#' then a plain text report will be generated automatically in the current
#' working directory with file name `test_vs_ref_TZ_YYYY-MM-DD_HHMMSS.txt`,
#' where `test` and `ref` are data set names of test and reference, `TZ` is the
#' time zone such as `CEST`, `YYYY-MM-DD` is the numeric date format and
#' `HHMMSS` is the numeric time format for hour, minute, and second.
#'
#' For a quick check, set argument `output.to.screen = TRUE`, a summary report
#' very similar to `concise` style report will be printed on screen.
#'
#' ## Types of Estimators
#' According to Shah et al, the population \eqn{f_2}{f2} for the inference is
#' \deqn{f_2 = 100-25\log\left(1 + \frac{1}{P}\sum_{i=1}^P%
#'   \left(\mu_{\mathrm{T},i} - \mu_{\mathrm{R},i}\right)^2 \right)\,,}{%
#'   f2 = 100 - 25 log(1 + 1/P(\sum(\mu(Ti) - \mu(Ri))^2)),}
#' where \eqn{P} is the number of time points; \eqn{\mu_{\mathrm{T},i}}{\mu(Ti)}
#' and \eqn{\mu_{\mathrm{R},i}}{\mu(Ri)} are *population mean* of test and
#' reference product at time point \eqn{i}, respectively; \eqn{\sum_{i=1}^P}{%
#' \sum} is the summation from \eqn{i = 1} to \eqn{P}.
#'
#' Five estimators for \eqn{f_2}{f2} are included in the function:
#' 1. The estimated \eqn{f_2}{f2}, denoted by \eqn{\hat{f}_2}{est.f2}, is the
#'    one written in various regulatory guidelines. It is expressed differently,
#'    but mathematically equivalently, as
#'    \deqn{\hat{f}_2 = 100-25\log\left(1 + \frac{1}{P}\sum_{i=1}^P\left(%
#'      \bar{X}_{\mathrm{T},i} - \bar{X}_{\mathrm{R},i}\right)^2 \right)\:,}{%
#'      est.f2 = 100 - 25 log(1 + 1/P(\sum(X(Ti) - X(Ri))^2)),}
#'    where \eqn{P} is the number of time points;
#'    \eqn{\bar{X}_{\mathrm{T},i}}{X(Ti)} and
#'    \eqn{\bar{X}_{\mathrm{R},i}}{X(Ri)} are mean dissolution data at the
#'    \eqn{i}th time point of *random samples* chosen from the test and the
#'    reference population, respectively. Compared to the equation of population
#'    \eqn{f_2}{f2} above, the only difference is that in the equation of
#'    \eqn{\hat{f}_2}{est.f2} the *sample means* of dissolution profiles replace
#'    the *population means* for the approximation. *In other words, a point
#'    estimate is used for the statistical inference in practice*.
#' 2. The Bias-corrected \eqn{f_2}{f2}, denoted by
#'    \eqn{\hat{f}_{2,\mathrm{bc}}}{bc.f2}, was described in the article of
#'    Shah et al, as
#'    \deqn{\hat{f}_{2,\mathrm{bc}} = 100-25\log\left(1 + \frac{1}{P}%
#'      \left(\sum_{i=1}^P\left(\bar{X}_{\mathrm{T},i} - %
#'      \bar{X}_{\mathrm{R},i}\right)^2 - \frac{1}{n}\sum_{i=1}^P%
#'      \left(S_{\mathrm{T},i}^2 + S_{\mathrm{R},i}^2\right)\right)\right)\,,}{%
#'      bc.f2 = 100 - 25 log(1 + 1/P(\sum(X(Ti) - X(Ri))^2 - 1/n\sum(S(Ti)^2 + %
#'      S(Ri)^2))),}
#'    where \eqn{S_{\mathrm{T},i}^2}{S(Ti)^2} and
#'    \eqn{S_{\mathrm{R},i}^2}{S(Ri)^2} are unbiased estimates of variance at
#'    the \eqn{i}th time points of random samples chosen from test and reference
#'    population, respectively; and \eqn{n} is the sample size.
#' 3. The variance- and bias-corrected \eqn{f_2}{f2}, denoted by
#'    \eqn{\hat{f}_{2,\mathrm{vcbc}}}{vc.bc.f2}, does not assume equal weight of
#'    variance as \eqn{\hat{f}_{2,\mathrm{bc}}}{bc.f2} does.
#'    \deqn{\hat{f}_{2, \mathrm{vcbc}} = 100-25\log\left(1 +%
#'      \frac{1}{P}\left(\sum_{i=1}^P \left(\bar{X}_{\mathrm{T},i} -%
#'        \bar{X}_{\mathrm{R},i}\right)^2 - \frac{1}{n}\sum_{i=1}^P%
#'        \left(w_{\mathrm{T},i}\cdot S_{\mathrm{T},i}^2 +%
#'        w_{\mathrm{R},i}\cdot S_{\mathrm{R},i}^2\right)\right)\right)\,,}{%
#'        vc.bc.f2 = 100 -25 log(1 + 1/P(\sum(X(Ti) - X(Ri))^2 - 1/n\sum(w(Ti)%
#'        S(Ti)^2 + w(Ri)S(Ri)^2))),}
#'    where \eqn{w_{\mathrm{T},i}}{w(Ti)} and \eqn{w_{\mathrm{R},i}}{w(Ri)} are
#'    weighting factors for variance of test and reference products,
#'    respectively, which can be calculated as follows:
#'    \deqn{w_{\mathrm{T},i} = 0.5 + \frac{S_{\mathrm{T},i}^2}%
#'      {S_{\mathrm{T},i}^2 + S_{\mathrm{R},i}^2}\,,}{w(Ti) = 0.5 + %
#'     S(Ti)^2/(S(Ti)^2 + S(Ri)^2),} and
#'    \deqn{w_{\mathrm{R},i} = 0.5 + \frac{S_{\mathrm{R},i}^2}%
#'      {S_{\mathrm{T},i}^2 + S_{\mathrm{R},i}^2}\,.}{w(Ri) = 0.5 + %
#'      S(Ri)^2/(S(Ti)^2 + S(Ri)^2).}
#' 4. The expected \eqn{f_2}{f2}, denoted by \eqn{\hat{f}_{2, \mathrm{exp}}}{%
#'    exp.f2}, is calculated based on the mathematical expectation of estimated
#'    \eqn{f_2}{f2},
#'    \deqn{\hat{f}_{2, \mathrm{exp}} = 100-25\log\left(1 + \frac{1}{P}%
#'      \left(\sum_{i=1}^P\left(\bar{X}_{\mathrm{T},i} - %
#'      \bar{X}_{\mathrm{R},i}\right)^2 + \frac{1}{n}\sum_{i=1}^P%
#'      \left(S_{\mathrm{T},i}^2 + S_{\mathrm{R},i}^2\right)\right)\right)\,,}{%
#'      exp.f2 = 100 - 25 log(1 + 1/P(\sum(X(Ti) - X(Ri))^2 + 1/n\sum(%
#'      S(Ti)^2 + S(Ri)^2))),}
#'    using mean dissolution profiles and variance from samples for the
#'    approximation of population values.
#' 5. The variance-corrected expected \eqn{f_2}{f2}, denoted by
#'    \eqn{\hat{f}_{2, \mathrm{vcexp}}}{vc.exp.f2}, is calculated as
#'    \deqn{\hat{f}_{2, \mathrm{vcexp}} = 100-25\log\left(1 +%
#'      \frac{1}{P}\left(\sum_{i=1}^P \left(\bar{X}_{\mathrm{T},i} -%
#'        \bar{X}_{\mathrm{R},i}\right)^2 + \frac{1}{n}\sum_{i=1}^P%
#'        \left(w_{\mathrm{T},i}\cdot S_{\mathrm{T},i}^2 +%
#'        w_{\mathrm{R},i}\cdot S_{\mathrm{R},i}^2\right)\right)\right)\,.}{%
#'       vc.exp.f2 = 100 - 25 log(1 + 1/P(\sum(X(Ti) - X(Ri))^2 + 1/n\sum(w(Ti)%
#'        S(Ti)^2 + w(Ri)S(Ri)^2))).}
#'
#' ## Types of Confidence Interval
#' The following confidence intervals are included:
#' 1. The Normal interval with bias correction, denoted by `normal` in the
#'    function, is estimated according to Davison and Hinkley,
#'    \deqn{\hat{f}_{2, \mathrm{L,U}} = \hat{f}_{2, \mathrm{S}} - E_B \mp %
#'      \sqrt{V_B}\cdot Z_{1-\alpha}\,,}{f2(L,U) = f2(S) - E(B) -/+ %
#'      sqrt(V(B))Z(1-\alpha)),}
#'    where \eqn{\hat{f}_{2, \mathrm{L,U}}}{f2(L,U)} are the lower and upper
#'    limit of the confidence interval estimated from bootstrap samples;
#'    \eqn{\hat{f}_{2, \mathrm{S}}}{f2(S)} denotes the estimators described
#'    above; \eqn{Z_{1-\alpha}}{Z(1-\alpha)} represents the inverse of standard
#'    normal cumulative distribution function with type I error \eqn{\alpha};
#'    \eqn{E_B}{E(B)} and \eqn{V_B}{V(B)} are the *resampling estimates* of bias
#'    and variance calculated as
#'    \deqn{E_B = \frac{1}{B}\sum_{b=1}^{B}\hat{f}_{2,b}^\star - %
#'      \hat{f}_{2, \mathrm{S}} = \bar{f}_2^\star - \hat{f}_{2,\mathrm{S}}\,,}{%
#'      E(B) = 1/B\sum(f2(b)) - f2(S) = f2(b,m) - f2(S),}
#'    and
#'    \deqn{V_B = \frac{1}{B-1}\sum_{b=1}^{B} \left(\hat{f}_{2,b}^\star
#'      -\bar{f}_2^\star\right)^2\,,}{V(B) = 1/(B-1)\sum(f2(b) - f2(b,m))^2,}
#'    where \eqn{B} is the number of bootstrap samples;
#'    \eqn{\hat{f}_{2,b}^\star}{f2(b)} is the \eqn{f_2}{f2} estimate with
#'    \eqn{b}th bootstrap sample, and \eqn{\bar{f}_2^\star}{f2(b,m)} is the
#'    mean value.
#' 2. The basic interval, denoted by `basic` in the function, is estimated
#'    according to Davison and Hinkley,
#'    \deqn{\hat{f}_{2, \mathrm{L}} = 2\hat{f}_{2, \mathrm{S}} -%
#'      \hat{f}_{2,(B+1)(1-\alpha)}^\star\,,}{f2(L) = 2*f2(S) - %
#'      f2((B+1)(1-\alpha)),}
#'    and
#'    \deqn{\hat{f}_{2, \mathrm{U}} = 2\hat{f}_{2, \mathrm{S}} -%
#'      \hat{f}_{2,(B+1)\alpha}^\star\,,}{f2(U) = 2*f2(S) - f2((B+1)\alpha),}
#'    where \eqn{\hat{f}_{2,(B+1)\alpha}^\star}{f2((B+1)\alpha)} and
#'    \eqn{\hat{f}_{2,(B+1)(1-\alpha)}^\star}{f2((B+1)(1-\alpha))} are the
#'    \eqn{(B+1)\alpha}th and \eqn{(B+1)(1-\alpha)}th *ordered resampling
#'    estimates* of \eqn{f_2}{f2}, respectively. When \eqn{(B+1)\alpha} is not
#'    an integer, the following equation is used for interpolation,
#'    \deqn{\hat{f}_{2,(B+1)\alpha}^\star = \hat{f}_{2,k}^\star + %
#'      \frac{\Phi^{-1}\left(\alpha\right)-\Phi^{-1}\left(\frac{k}{B+1}\right)}%
#'      {\Phi^{-1}\left(\frac{k+1}{B+1}\right)-\Phi^{-1}%
#'      \left(\frac{k}{B+1}\right)}\left(\hat{f}_{2,k+1}^\star-%
#'      \hat{f}_{2,k}^\star\right),}{f2((B+1)\alpha) = f2(k) + %
#'      (\Phi^(-1)(\alpha) - \Phi^(-1)(k/(B+1)))/(\Phi^(-1)((k+1)/(B+1)) - %
#'      \Phi^(-1)(k/(B+1)))*(f2(k+1) - f2(k)),}
#'    where \eqn{k} is the *integer part* of \eqn{(B+1)\alpha},
#'    \eqn{\hat{f}_{2,k+1}^\star}{f2(k+1)} and \eqn{\hat{f}_{2,k}^\star}{%
#'    f2(k)} are the \eqn{(k+1)}th and the \eqn{k}th ordered resampling
#'    estimates of \eqn{f_2}{f2}, respectively.
#' 3. The percentile intervals, denoted by `percentile` in the function, are
#'    estimated using nine different types of quantiles, Type 1 to Type 9, as
#'    summarized in Hyndman and Fan's article and implemented in `R`'s
#'    `quantile` function. Using `R`'s `boot` package, program `bootf2BCA`
#'    outputs a percentile interval using the equation above for interpolation.
#'    To be able to compare the results among different programs, the same
#'    interval, denoted by `Percentile (boot)` in the function, is also
#'    included in the function.
#' 4. The bias-corrected and accelerated (BCa) intervals are estimated according
#'    to Efron and Tibshirani,
#'    \deqn{\hat{f}_{2, \mathrm{L}} = \hat{f}_{2, \alpha_1}^\star\,,}{%
#'    f2(L) = f2(\alpha1),}
#'    \deqn{\hat{f}_{2, \mathrm{U}} = \hat{f}_{2, \alpha_2}^\star\,,}{%
#'    f2(L) = f2(\alpha2),}
#'    where \eqn{\hat{f}_{2, \alpha_1}^\star}{f2(\alpha1)} and
#'    \eqn{\hat{f}_{2, \alpha_2}^\star}{f2(\alpha2)} are the \eqn{100\alpha_1}{%
#'    100\alpha1}th and the \eqn{100\alpha_2}{100\alpha2}th percentile of the
#'    resampling estimates of \eqn{f_2}{f2}, respectively. Type I errors
#'    \eqn{\alpha_1}{\alpha1} and \eqn{\alpha_2}{\alpha2} are obtained as
#'    \deqn{\alpha_1 = \Phi\left(\hat{z}_0 + \frac{\hat{z}_0 + \hat{z}_\alpha}%
#'      {1-\hat{a}\left(\hat{z}_0 + \hat{z}_\alpha\right)}\right),}{\alpha1 = %
#'      \Phi(z0 + (z0 + za)/(1 - a(z0 + za))),}
#'    and
#'    \deqn{\alpha_2 = \Phi\left(\hat{z}_0 + \frac{\hat{z}_0 + %
#'      \hat{z}_{1-\alpha}}{1-\hat{a}\left(\hat{z}_0 + %
#'      \hat{z}_{1-\alpha}\right)}\right),}{\alpha2 = \Phi(z0 +
#'      (z0 + z(1-\alpha))/(1 - a(z0 + z(1 - \alpha)))),}
#'    where \eqn{\hat{z}_0}{z0} and \eqn{\hat{a}}{a} are called
#'    *bias-correction* and *acceleration*, respectively.
#'
#'    There are different methods to estimate \eqn{\hat{z}_0}{z0} and
#'    \eqn{\hat{a}}{a}. Shah et al. used jackknife technique, denoted by
#'    `bca.jackknife` in the function,
#'    \deqn{\hat{z}_0 = \Phi^{-1}\left(\frac{N\left\{\hat{f}_{2,b}^\star <%
#'      \hat{f}_{2,\mathrm{S}} \right\}}{B}\right),}{z0 = %
#'      \Phi^(-1)(N(f2(b) < f2(S))/B)}
#'    and
#'    \deqn{\hat{a} = \frac{\sum_{i=1}^{n}\left(\hat{f}_{2,\mathrm{m}} -%
#'      \hat{f}_{2, i}\right)^3}{6\left(\sum_{i=1}^{n}\left(%
#'      \hat{f}_{2,\mathrm{m}} - \hat{f}_{2, i}\right)^2\right)^{3/2}}\,,}{%
#'      a = (\sum(f2(m) - f2(i)))^3/(6(\sum(f2(m) - f2(i))^2)^(3/2)),}
#'    where \eqn{N\left\{\cdot\right\}}{N(f2(b) < f2(S))} denotes the number of
#'    element in the set, \eqn{\hat{f}_{2, i}}{f2(i)} is the \eqn{i}th
#'    jackknife statistic, \eqn{\hat{f}_{2,\mathrm{m}}}{f2(m)} is the mean of
#'    the jackknife statistics, and \eqn{\sum} is the summation from 1 to
#'    sample size \eqn{n}.
#'
#'    Program `bootf2BCA` gives a slightly different BCa interval with `R`'s
#'    `boot` package. This approach, denoted by `bca.boot` in the function, is
#'    also implemented in the function for estimating the interval.
#'
#' ## Notes on the argument `jackknife.type`
#' For any sample with size $n$, the jackknife estimator is obtained by
#' estimating the parameter for each subsample omitting the \eqn{i}th
#' observation. However, when two samples (e.g., test and reference) are
#' involved, there are several possible ways to do it. Assuming sample size
#' of test and reference are \eqn{n_\mathrm{T}}{nt} and \eqn{n_\mathrm{R}}{nr},
#' the following three possibility are considered:
#' - Estimated by removing one observation from both test and reference samples.
#'   In this case, the prerequisite is \eqn{n_\mathrm{T} = n_\mathrm{R}}{nt=nr},
#'   denoted by `nt=nr` in the function. So if there are 12 units in test and
#'   reference data sets, there will be 12 jackknife estimators.
#' - Estimate the jackknife for test sample while keeping the reference data
#'   unchanged; and then estimate jackknife for reference sample while keeping
#'   the test sample unchanged. This is denoted by `nt+nr` in the function.
#'   This is the default method. So if there are 12 units in test and reference
#'   data sets, there will be \eqn{12 + 12 = 24} jackknife estimators.
#' - For each observation deleted from test sample, estimate jackknife for
#'   reference sample. This is denoted by `nt*nr` in the function. So if there
#'   are 12 units in test and reference data sets, there will be \eqn{12 \times
#'   12 = 144}{12*12 = 144} jackknife estimators.
#'
#'
#' @references Shah, V. P.; Tsong, Y.; Sathe, P.; Liu, J.-P. In Vitro
#'   Dissolution Profile Comparison---Statistics and Analysis of the
#'   Similarity Factor, \eqn{f_2}{f2}. *Pharmaceutical Research* 1998,
#'   **15** (6), 889--896. DOI: 10.1023/A:1011976615750.
#' @references Davison, A. C.; Hinkley, D. V. Bootstrap Methods and Their
#'   Application. Cambridge University Press, 1997.
#' @references Hyndman, R. J.; Fan, Y. Sample Quantiles in Statistical Packages.
#'   *The American Statistician* 1996, **50** (4), 361--365. DOI:
#'   /10.1080/00031305.1996.10473566.
#' @references Efron, B.; Tibshirani, R. An Introduction to the Bootstrap.
#'   Chapman & Hall, 1993.
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
#' bootf2(dt$sim.disso, dr$sim.disso, n.boots = 100, print.report = FALSE)
#' }
#'
#' @export
bootf2 <- function(test, ref, path.in, file.in, path.out, file.out,
                   n.boots = 10000L, seed = 306L, digits = 2L, alpha = 0.05,
                   regulation = c("EMA", "FDA", "WHO", "Canada", "ANVISA"),
                   min.points = 1L, both.TR.85 = FALSE, print.report = TRUE,
                   report.style = c("concise",  "intermediate", "detailed"),
                   f2.type = c("all", "est.f2", "exp.f2", "bc.f2",
                               "vc.exp.f2", "vc.bc.f2"),
                   ci.type = c("all", "normal", "basic", "percentile",
                               "bca.jackknife", "bca.boot"),
                   quantile.type = c("all", as.character(1:9), "boot"),
                   jackknife.type = c("all", "nt+nr", "nt*nr", "nt=nr"),
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

  if (any(missing(test), missing(ref))) {
    stop("Both 'test' and 'ref' have to be specified.")
  }

  if (all(isTRUE(both.TR.85), regulation != "FDA")) {
    stop("'both.TR.85 = TRUE' is only applicable when 'regulation = FDA'.")
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

  nt <- NCOL(data.t) - 1
  nr <- NCOL(data.r) - 1

  if (is.null(seed)) {
    seed <- sample(1:(.Machine$integer.max - 1), 1)
  }
  set.seed(seed)

  # f2 original ------------------------------------------------------
  # calculate f2s with original data without regard to variability
  tmp.f2o <- calcf2(test = data.t, ref = data.r, regulation = regulation,
                    digits = digits, cv.rule = FALSE, min.points = min.points,
                    both.TR.85 = both.TR.85, message = FALSE, f2.type = f2.type,
                    plot = FALSE, time.unit = time.unit)
  f2o <- c(tmp.f2o$f2.value, unique(tmp.f2o$f2.tp),
           unique(ifelse(tmp.f2o$d85at15 == "yes", 1, 0)))
  names(f2o) <- c(tmp.f2o$f2.type, "f2.tp", "d85at15")

  # bootstrap data -------------------------------------------------------------
  # bootstrap index --------------------------------------------------
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

  # initialize result ------------------------------------------------
  boot.t <- boot.r  <- vector(mode = "list", length = n.boots)
  boot.f2 <- matrix(NA, nrow = n.boots, ncol = length(f2o), byrow = TRUE,
                    dimnames = list(rep("", n.boots), names(f2o)))

  for (i in 1:n.boots) {
    boot.t[[i]]  <- data.t[, c(1, bt.array[i, 1:nt] + 1)]
    boot.r[[i]]  <- data.r[, c(1, bt.array[i, (nt + 1):(nt + nr)] - nt + 1)]
    tmp.f2       <- calcf2(test = boot.t[[i]], ref = boot.r[[i]],
                           regulation = regulation, digits = digits,
                           cv.rule = FALSE, min.points = min.points,
                           both.TR.85 = both.TR.85, message = FALSE,
                           f2.type = f2.type, plot = FALSE,
                           time.unit = time.unit)
    boot.f2[i, ] <- c(tmp.f2$f2.value, unique(tmp.f2$f2.tp),
                      unique(ifelse(tmp.f2$d85at15 == "yes", 1, 0)))
  }

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
  if (ci.type %in% c("all", "normal")) {
    normal.ci <- function(boot.f2, f2o, alpha) {
      btf2  <- as.vector(boot.f2)
      btf2  <- btf2[is.finite(btf2)]
      btf2.mean <- mean(btf2, na.rm = TRUE)
      btf2.var  <- var(btf2, na.rm = TRUE)
      2*f2o - btf2.mean + c(-1, 1)*sqrt(btf2.var)*qnorm(1 - alpha)
    }

    ci.normal <- data.frame(f2.type = NA, ci.type = NA,
                            ci.lower = NA, ci.upper = NA)

    for (i in seq_len(NCOL(boot.f2) - 2)) {
      ci.normal[i, 1]   <- dimnames(boot.f2)[[2]][i]
      ci.normal[i, 2]   <- "Normal"
      ci.normal[i, 3:4] <- normal.ci(boot.f2 = boot.f2[, i], f2o = f2o[[i]],
                                     alpha = alpha)
    }
  }# end normal CI

  # basic interval ---------------------------------------------------
  # Basic CI. Davidson, Ch5, Eq 5.6
  if (ci.type %in% c("all", "basic")) {
    basic.ci <- function(boot.f2, f2o, alpha) {
      btf2  <- as.vector(boot.f2)
      btf2  <- btf2[is.finite(btf2)]
      2*f2o - normal.inter(boot.f2 = btf2, alpha = alpha)
    }

    ci.basic <- data.frame(f2.type = NA, ci.type = NA,
                           ci.lower = NA, ci.upper = NA)

    for (i in seq_len(NCOL(boot.f2) - 2)) {
      ci.basic[i, 1]   <- dimnames(boot.f2)[[2]][i]
      ci.basic[i, 2]   <- "Basic"
      ci.basic[i, 3:4] <- basic.ci(boot.f2 = boot.f2[, i], f2o = f2o[i],
                                   alpha = c(1 - alpha, alpha))
    }
  } # end basic interval

  # percentile interval ----------------------------------------------
  if (ci.type %in% c("all", "percentile")) {
    ci.percentile <- data.frame(f2.type = NA, ci.type = NA,
                                ci.lower = NA, ci.upper = NA)

    if (quantile.type == "boot") {# same as boot package
      for (i in seq_len(NCOL(boot.f2) - 2)) {
        ci.percentile[i, 1]   <- dimnames(boot.f2)[[2]][[i]]
        ci.percentile[i, 2]   <- "Percentile (boot)"
        ci.percentile[i, 3:4] <- normal.inter(boot.f2 = boot.f2[, i],
                                              alpha = c(alpha, 1 - alpha))
      }
    } else if (quantile.type == "all") {
      k <- 0
      for (i in seq_len(NCOL(boot.f2) - 2)) {
        for (j in seq_along(q.type)) {
          k <- k + 1
          ci.percentile[k, 1]   <- dimnames(boot.f2)[[2]][[i]]
          ci.percentile[k, 2]   <- paste0("Percentile (Type ", j, ")")
            # ifelse(quantile.type == "all",
            #        paste0("Percentile (Type ", j, ")"),
            #        paste0("Percentile (Type ", q.type, ")"))

          ci.percentile[k, 3:4] <- quantile(boot.f2[, i],
                                            probs = c(alpha, 1 - alpha),
                                            na.rm = TRUE, names = FALSE,
                                            type = q.type[[j]])
        } # end quantile type loop j

        # Davison, Ch5, same as boot package percentile ci
        k <- k + 1
        ci.percentile[k, 1]   <- dimnames(boot.f2)[[2]][[i]]
        ci.percentile[k, 2]   <- "Percentile (boot)"
        ci.percentile[k, 3:4] <- normal.inter(boot.f2 = boot.f2[, i],
                                              alpha = c(alpha, 1 - alpha))
      } # end f2 type loop i
    } else {# type = 1, 2, ..., 9
      for (i in seq_len(NCOL(boot.f2) - 2)) {
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
  if (ci.type %in% c("all", "bca.jackknife")) {
    # function to obtain acceleration a by jackkinfe
    jackf2 <- function(data.t, data.r) {
      # remove time 0 point if any
      data.t <- data.t[data.t[, 1] != 0, ]
      data.r <- data.r[data.r[, 1] != 0, ]

      if (jackknife.type %in% c("all", "nt+nr")) {
        jk1.f2 <- matrix(NA, nrow = nt + nr, ncol = length(f2o) - 2,
                         dimnames = list(rep("", nt + nr),
                                         names(f2o)[1:(length(f2o) - 2)]))
        # jackknife with test data
        for (i in 1:nt) {
          jk1.f2[i, ] <- calcf2(test = data.t[, -(i + 1)], ref = data.r,
                                regulation = regulation, digits = digits,
                                cv.rule = FALSE, min.points = min.points,
                                both.TR.85 = both.TR.85, message = FALSE,
                                f2.type = f2.type, plot = FALSE,
                                time.unit = time.unit)$f2.value
        }
        # jackknife with reference data
        for (j in 1:nr) {
          jk1.f2[i + j, ] <- calcf2(test = data.t, ref = data.r[, -(j + 1)],
                                   regulation = regulation, digits = digits,
                                   cv.rule = FALSE, min.points = min.points,
                                   both.TR.85 = both.TR.85, message = FALSE,
                                   f2.type = f2.type, plot = FALSE,
                                   time.unit = time.unit)$f2.value
        }
        jk1.f2.mean <- as.data.frame(t(colMeans(jk1.f2)),
                                     stringsAsFactors = FALSE)
        jk1.f2.mean$type <- "nt+nr"
      }

      if (jackknife.type %in% c("all", "nt*nr")) {
        jk2.f2 <- matrix(NA, nrow = nt*nr, ncol = length(f2o) - 2,
                         dimnames = list(rep("", nt*nr),
                                         names(f2o)[1:(length(f2o) - 2)]))
        k <- 0
        for (i in 1:nt) {
          for (j in 1:nr) {
            k <- k + 1
            jk2.f2[k, ] <- calcf2(test = data.t[, -(i + 1)],
                                  ref = data.r[, -(j + 1)],
                                  regulation = regulation, digits = digits,
                                  cv.rule = FALSE, min.points = min.points,
                                  both.TR.85 = both.TR.85, message = FALSE,
                                  f2.type = f2.type, plot = FALSE,
                                  time.unit = time.unit)$f2.value
          }
        }
        jk2.f2.mean <- as.data.frame(t(colMeans(jk2.f2)),
                                     stringsAsFactors = FALSE)
        jk2.f2.mean$type <- "nt*nr"
      }

      if (jackknife.type %in% c("all", "nt=nr")) {# need nt = nr
        if (!all.equal(nt, nr)) {# usu. nt = nr = 12, not a problem
          stop("To use this type of jackknife, the number of test and ",
               "reference\ndata should be equal.")
        } else {
          jk3.f2 <- matrix(NA, nrow = nt, ncol = length(f2o) - 2,
                           dimnames = list(rep("", nt),
                                           names(f2o)[1:(length(f2o) - 2)]))
          for (i in 1:nt) {
            jk3.f2[i, ] <- calcf2(test = data.t[, -(i + 1)],
                                  ref = data.r[, -(i + 1)],
                                  regulation = regulation, digits = digits,
                                  cv.rule = FALSE, min.points = min.points,
                                  both.TR.85 = both.TR.85, message = FALSE,
                                  f2.type = f2.type, plot = FALSE,
                                  time.unit = time.unit)$f2.value
          }
        }
        jk3.f2.mean <- as.data.frame(t(colMeans(jk3.f2)),
                                     stringsAsFactors = FALSE)
        jk3.f2.mean$type <- "nt=nr"
      }

      # accelerated alpha and jackknife mean
      if (jackknife.type == "all") {
        jk.f2.mean <- rbind(jk1.f2.mean, jk2.f2.mean, jk3.f2.mean)
        jk.f2 <- list(jk1.f2, jk2.f2, jk3.f2)
        a <- matrix(NA, nrow = NROW(jk.f2.mean), ncol = NCOL(jk.f2.mean) - 1,
                    dimnames = list(
                      rep("", NROW(jk.f2.mean)),
                      paste0(names(jk.f2.mean)[1:(NCOL(jk.f2.mean) - 1)], ".a"))
        )

        for (i in 1:NROW(a)) {
          for (j in 1:NCOL(a)) {
            a[i, j] <- sum((jk.f2.mean[i, j] - jk.f2[[i]][, j])^3, na.rm=TRUE)/
              (6*(sum((jk.f2.mean[i, j] - jk.f2[[i]][, j])^2, na.rm=TRUE))^1.5)
          }
        }
      } else {
        if (jackknife.type == "nt+nr") {
          jk.f2 <- jk1.f2
          jk.f2.mean <- jk1.f2.mean
        } else if (jackknife.type == "nt*nr") {
          jk.f2 <- jk2.f2
          jk.f2.mean <- jk2.f2.mean
        } else {
          jk.f2 <- jk3.f2
          jk.f2.mean <- jk3.f2.mean
        }
        a <- matrix(NA, nrow = NROW(jk.f2.mean), ncol = NCOL(jk.f2.mean) - 1,
                    dimnames = list(
                      rep("", NROW(jk.f2.mean)),
                      paste0(names(jk.f2.mean)[1:(NCOL(jk.f2.mean) - 1)], ".a"))
        )

        for (i in 1:NCOL(a)) {
          a[, i] <- sum((jk.f2.mean[, i] - jk.f2[, i])^3, na.rm = TRUE)/
            (6*(sum((jk.f2.mean[, i] - jk.f2[, i])^2, na.rm = TRUE))^1.5)
        }
      }

      a <- as.data.frame(a)
      # jk.f2.mean <- jk.f2.mean[is.finite(a)]
      # a <- a[is.finite(a)]
      return(cbind(a, jk.f2.mean))
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
        if (all(is.finite(a1), is.finite(a2))) {
          return(normal.inter(boot.f2 = btf2, alpha = c(a1, a2)))
        } else return(c(NA, NA))
      } else return(c(NA, NA))
    }

    ci.bca.jackknife <- data.frame(f2.type = NA, ci.type = NA,
                                   ci.lower = NA, ci.upper = NA)

    k <- 0
    for (i in 1:NROW(a.jack)) {
      for (j in seq_len(NCOL(boot.f2) - 2)) {
        k <- k + 1
        ci.bca.jackknife[k, 1]   <- dimnames(boot.f2)[[2]][j]
        ci.bca.jackknife[k, 2]   <- paste0("BCa (jackknife, ",
                                           a.jack[i, "type"], ")")
        ci.bca.jackknife[k, 3:4] <- bca.ci.jack(boot.f2 = boot.f2[, j],
                                                f2o = f2o[[j]],
                                                a = a.jack[i, j],
                                                alpha = alpha)
      }
    }
  } # end BCa CI by jackknife

  # BCa, boot --------------------------------------------------------
  # BCa, by empirical regression, Davison, Sec 2.7.4
  # also, ref boot package, empinf.reg function. same as Mendyk's bootf2BCA
  if (ci.type %in% c("all", "bca.boot")) {
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
        if (all(is.finite(a1), is.finite(a2))) {
          return(normal.inter(boot.f2 = btf2, c(a1, a2)))
        } else return(c(NA, NA))
      } else return(c(NA, NA))
    } # end function bca.ci.boot

    ci.bca.boot <- data.frame(f2.type = NA, ci.type = NA,
                              ci.lower = NA, ci.upper = NA)

    for (i in seq_len(NCOL(boot.f2) - 2)) {
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
    } else {# if file.out provided with non-txt wrong extension
      file.out1 <-
        ifelse(regmatches(file.out, regexpr(".{3}$", file.out)) == "txt",
                          file.out, paste0(file.out, ".txt"))
    }
    file.out <- paste0(path.out, file.out1)
  } else {# no report
    path.out <- file.out <- file.out1 <- NA
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
    both.TR.85          = both.TR.85,
    f2.type             = f2.type,
    ci.type             = ci.type,
    quantile.type       = ifelse(ci.type %in% c("all", "percentile"),
                                 quantile.type,
                                 paste0(quantile.type, " (not used)")),
    jackknife.type      = ifelse(ci.type %in% c("all", "bca.jackknife"),
                                 jackknife.type,
                                 paste0(jackknife.type, " (not used)")),
    time.unit           = time.unit,
    digits              = digits,
    print.report        = print.report,
    report.style        = ifelse(isTRUE(print.report), report.style,
                                 paste0(report.style, " (not used)")),
    output.to.screen    = output.to.screen,
    sim.data.out        = sim.data.out,
    path.in             = ifelse(is.na(path.in), "Not provided", path.in),
    file.in             = ifelse(is.na(file.in1), "Not provided", file.in1),
    path.out            = ifelse(is.na(path.out), "Not provided", path.out),
    file.out            = ifelse(is.na(file.out1), "Not provided", file.out1),
    stringsAsFactors = FALSE
  )

  btf2.type <- dimnames(boot.f2)[[2]][1:(NCOL(boot.f2) - 2)]
  bt.mean <- apply(boot.f2[, 1:(NCOL(boot.f2) - 2), drop = FALSE], 2, mean,
                   na.rm = TRUE)
  bt.median <- apply(boot.f2[, 1:(NCOL(boot.f2) - 2), drop = FALSE], 2, median,
                     na.rm = TRUE)
  bt.na <- apply(boot.f2[, 1:(NCOL(boot.f2) - 2), drop = FALSE], 2,
                 function(x) sum(is.na(x)))
  f2.tp1 <- sum(boot.f2[, "f2.tp"] == 1)
  f2.tp2 <- sum(boot.f2[, "f2.tp"] == 2)
  d85at15 <- sum(boot.f2[, "d85at15"])

  btsum <- data.frame(f2.type = btf2.type, boot.mean = bt.mean,
                      boot.median = bt.median, boot.na = bt.na,
                      f2.tp1 = f2.tp1, f2.tp2 = f2.tp2, d85at15 = d85at15,
                      f2o = f2o[1:(length(f2o) - 2)],
                      row.names = 1:length(btf2.type),
                      stringsAsFactors = FALSE)

  # write report ---------------------------------------------------------------
  if (isTRUE(output.to.screen))
    rpt.screen(boot.f2.ci, boot.info, f2o, a.jack, btsum)

  if (isTRUE(print.report)) {
    sink(file.out, split = FALSE)

    if (report.style %in% c("concise", "intermediate", "detailed")) {
      rpt.concise(boot.f2.ci, boot.info, f2o, a.jack, btsum)
    }

    if (report.style %in%  c("intermediate", "detailed")) {
      rpt.intermediate(boot.info, boot.f2)
    } # end intermediate

    # this part only for detailed report. long process time. not recommended
    if (report.style == "detailed") {
      rpt.detailed(data.t, data.r, boot.t, boot.r, boot.f2, boot.info, f2o)
    } # end detailed report
    sink()
  } # end isTRUE(print.report)

  if (isTRUE(sim.data.out)) {
    invisible(list(boot.ci = boot.f2.ci[order(boot.f2.ci[, 1]), ],
                   boot.f2 = as.data.frame(boot.f2, row.names = ""),
                   boot.info = boot.info, boot.summary = btsum,
                   boot.t = boot.t, boot.r = boot.r))
  } else {# save some memory
    invisible(list(boot.ci = boot.f2.ci[order(boot.f2.ci[, 1]), ],
                   boot.f2 = as.data.frame(boot.f2, row.names = ""),
                   boot.info = boot.info,  boot.summary = btsum))
  }
}
