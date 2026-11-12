## ---------------------------------------------------------------------------
## The rank-ordered (exploded) logit and the test of the stability of
## preferences across ranking stages. Kernel: src/rologit_ll.cpp; design:
## the long layout of .pc_design() with the ranks in D$str$rank.
## ---------------------------------------------------------------------------

#' Rank-ordered (exploded) logit regression
#'
#' Fits the rank-ordered logit (Beggs, Cardell and Hausman 1981, who call it
#' the ordered logit), also called the exploded logit after Chapman and
#' Staelin's (1982) "explosion" of a ranking into independent choices (Allison
#' and Christakis 1994): the model for data in which each case ranks some or all of its available
#' alternatives, such as survey rankings of parties or candidates, or the
#' ballots of a ranked-choice election. With the utilities of [mnl()], a
#' ranking is a sequence of multinomial-logit choices: the first-ranked
#' alternative is chosen from the whole choice set, the second from those left,
#' and so on. Rankings carry more information than first choices, so the same
#' coefficients are estimated more precisely than from the first choices alone.
#'
#' @details
#' Case \eqn{i} ranks the alternatives \eqn{r_{i1}, r_{i2}, \dots} of its choice
#' set \eqn{A_i}; with \eqn{R_{il}} the alternatives not ranked before position
#' \eqn{l},
#' \deqn{\log L_i = \sum_l \Big[ V_{i r_{il}} - \log \sum_{m \in R_{il}} e^{V_{im}} \Big],}
#' with \eqn{V_{ij}} as in [mnl()]. The log-likelihood is globally concave.
#' Partial rankings (the first few positions ranked, the rest left unranked)
#' are handled exactly: the unranked alternatives are tied below every ranked
#' one and contribute only through the denominators. Ties among ranked
#' alternatives are handled by the exact likelihood of Allison and Christakis
#' (1994, eq. 8): the tied alternatives are taken to have an unknown strict
#' order, and the likelihood sums over every order (at most 7 tied
#' alternatives above others; Breslow's approximation, the default of
#' proportional-hazards software, attenuates the estimates and is not used).
#'
#' The model assumes that the relative preference for two alternatives does
#' not depend on which alternatives have already been ranked: independence of
#' irrelevant alternatives at every stage. Respondents may rank the first
#' positions more carefully than the later ones: Chapman and Staelin (1982)
#' and Hausman and Ruud (1987) argue, and find in their applications, that
#' later ranks can be noisier, Hausman and Ruud adding that the pattern
#' depends on the survey. [rank_test()] tests whether the coefficients differ
#' between early and late stages, and `depth` limits the fit to the first
#' positions (Chapman and Staelin's explosion depth, Hausman and Ruud's top-P
#' estimator; with only the first position, the model is the multinomial logit
#' of the first choices).
#'
#' The data are in the long layout: one row per case and available alternative,
#' with `id` and `alt` naming the case and alternative columns, as in
#' [mnl()]. The response is the rank: smaller numbers are preferred, equal
#' numbers are ties (the numbers need only be in the right order), and 0 marks
#' an alternative left unranked.
#' * An unranked alternative is in the choice set, below every ranked one: it
#'   enters the denominators of every stage (the respondent saw it and did not
#'   rank it).
#' * An alternative without a row for the case is unavailable and enters none.
#'   Drop a row only where the alternative was not on offer.
#' * Without a penalty, ties are exact up to 7 alternatives tied above others.
#'   A larger tie group above other alternatives is refused, and a tie at the
#'   bottom, below no other alternative, is unlimited.
#'
#' The formula is `rank ~ w | x`:
#' alternative-specific covariates with generic coefficients, then
#' case-specific covariates with alternative-specific coefficients. Cases in
#' which no alternative is ranked above another carry no information and are
#' dropped with a message. Weights, offsets, clusters, and standard errors are
#' as in [mnl()]; a case (one ranking) is the observation. Unit fixed effects
#' are not available.
#'
#' Fitted values, [predict()], [first_difference()], and [rrr()] describe
#' first choices: the multinomial-logit probabilities over the choice set.
#' `predict(type = "ranking")` gives the probability of each observed ranking.
#' `y` holds each case's first choice (with a tie at the top, the first of the
#' tied alternatives in level order), which is what [classification()], the
#' count measures of [fit_stats()], and [residuals()] use; the log-likelihood
#' measures of [fit_stats()] are those of the rankings, against the
#' rank-ordered logit with alternative-specific constants only. [simulate()]
#' draws rankings. The tests of the independence of irrelevant alternatives
#' and the goodness-of-fit test apply to single-choice fits and refuse a
#' rank-ordered one; [lr_test()] and [vuong()] compare two rank-ordered fits of
#' the same rankings.
#'
#' Separation occurs in rankings as in single choices: when a combination of
#' the covariates ranks every observed choice first among the alternatives
#' still unranked at its position, the maximum-likelihood estimates run off
#' (the monotone likelihood of the Cox model, of which this is the version
#' stratified by case; Bryson and Johnson 1981). Ranking fewer positions
#' (`depth`, or partial rankings) leaves fewer such constraints and makes it
#' more likely. The maximum-likelihood fit decides it exactly (the linear
#' program of [separation()], ties as groups) and names the diverging
#' coefficients in a warning; an alternative never ranked above another
#' available one is the simplest case (its constant runs to minus infinity).
#' `firth = TRUE` penalizes the log-likelihood by half the
#' log-determinant of its information, the sum over the stages of each
#' ranking of the information of the stage's choice given the alternatives
#' left: Heinze and Schemper's (2001) penalty for the Cox model, with one
#' stratum per case. A tie group of \eqn{d} alternatives enters with the
#' information of the exact tie probability (the probability that the tied
#' set is the top \eqn{d} of the alternatives left), the expected information
#' of that outcome over the possible sets, which reduces to the
#' conditional-logit information of an untied position. The estimates are
#' finite under separation; the rest is as for [mnl()] (section "Penalized
#' (Firth) estimation"): the log-likelihood at the penalized estimates,
#' standard errors from the inverse information, the separated coefficients
#' in `separation`, the run from a second start, profile intervals from [confint()], and
#' penalized likelihood-ratio tests from [term_test()], [lr_test()] and
#' [rank_test()]. Under the penalty or a prior a tie group is refused when
#' more than 62 alternatives remain at its positions, when it has more than
#' \eqn{10^5} subsets of them to enumerate, or (Firth's penalty) when the
#' second derivatives of its subsets would need more than 200 MB; without a
#' penalty only ties of more than 7 alternatives ranked above others are
#' refused.
#'
#' @inheritParams mnl
#' @param firth If `TRUE`, Firth's penalized likelihood (see Details); default
#'   `FALSE`, maximum likelihood. Weights are read as frequency weights and
#'   must total more than 2.
#' @param prior Optional: priors for the coefficients, made by [coef_prior()];
#'   the fit is the posterior mode, a penalized fit like `firth = TRUE`. Not
#'   with `firth = TRUE`.
#' @param formula `rank ~ w | x` (see Details); `rank ~ w` without
#'   case-specific covariates.
#' @param id,alt Names of the case and alternative columns of `data`.
#' @param ref The base alternative: a level of the alternative column `alt`
#'   (default: its first level), whose constant and case-specific
#'   coefficients are zero.
#' @param offset Optional offset on the utility scale: a column name of `data`
#'   or a vector with one entry per row.
#' @param se Standard errors: `"analytic"` (inverse observed information,
#'   default), `"opg"`, `"robust"`, `"cluster"` (needs `cluster`),
#'   `"bootstrap"` (nonparametric over cases, or over clusters when `cluster`
#'   is given, with `nboot` refits), or `"none"`.
#' @param start Optional starting values: the full parameter vector in the order
#'   of `coef()` (`object$theta` of an earlier fit); rarely needed, since the
#'   log-likelihood is globally concave.
#' @param depth Optional: use only the first `depth` positions of every
#'   ranking; alternatives ranked below become unranked (a tie group that
#'   starts within the first `depth` positions is kept whole). `depth = 1` is
#'   the multinomial logit of the first choices.
#' @return An object of class `c("rologit", "choicemodel")`, with the
#'   components of an [mnl()] fit; `design$str$rank` holds the rankings used
#'   (dense ranks by case, 0 for unranked or unavailable alternatives).
#' @references
#' Allison, P.D. and Christakis, N.A. (1994). Logit models for sets of ranked
#' items. \emph{Sociological Methodology}, 24, 199-228. \doi{10.2307/270983}
#'
#' Beggs, S., Cardell, S. and Hausman, J. (1981). Assessing the potential
#' demand for electric cars. \emph{Journal of Econometrics}, 17, 1-19.
#' \doi{10.1016/0304-4076(81)90056-7}
#'
#' Bryson, M.C. and Johnson, M.E. (1981). The incidence of monotone likelihood
#' in the Cox model. \emph{Technometrics}, 23, 381-383.
#' \doi{10.1080/00401706.1981.10487683}
#'
#' Chapman, R.G. and Staelin, R. (1982). Exploiting rank ordered choice set
#' data within the stochastic utility model. \emph{Journal of Marketing
#' Research}, 19, 288-301. \doi{10.1177/002224378201900302}
#'
#' Hausman, J.A. and Ruud, P.A. (1987). Specifying and testing econometric
#' models for rank-ordered data. \emph{Journal of Econometrics}, 34, 83-104.
#' \doi{10.1016/0304-4076(87)90068-6}
#'
#' Heinze, G. and Schemper, M. (2001). A solution to the problem of monotone
#' likelihood in Cox regression. \emph{Biometrics}, 57, 114-119.
#' \doi{10.1111/j.0006-341X.2001.00114.x}
#' @examples
#' set.seed(4)
#' d <- rrologit(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
#'               gamma = 0.7, depth = 3)
#' head(d)
#' m <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
#' summary(m)
#' rank_test(m, after = 1)
#' @seealso [rank_test()], [rrologit()], [mnl()]
#' @family estimators
#' @export
rologit <- function(formula, data, id, alt, ref = NULL, weights = NULL, offset = NULL, depth = NULL,
                    se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                    cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10, firth = FALSE, prior = NULL) {
  if (missing(id) || missing(alt)) stop("rologit() takes the long layout: give 'id' and 'alt' (one row per case and alternative; the response is the rank).", call. = FALSE)
  se_missing <- missing(se); se <- match.arg(se)
  if (!is.null(depth)) depth <- .pc_check_whole(depth, "depth", 1)
  .pc_fit(formula, data, ref = ref, id = id, alt = alt, weights = weights, offset = offset, cluster = cluster,
          se = se, nboot = nboot, cores = cores, start = start, maxit = maxit, reltol = reltol, call = match.call(),
          se_missing = se_missing, family = "rologit", str = list(depth = depth), firth = firth, prior = prior)
}

#' Stability of preferences across ranking stages
#'
#' Tests the assumption of the rank-ordered logit that the relative preference
#' for two alternatives is the same at every stage of the ranking: the model is
#' refit with a second set of coefficients for the stages after position
#' `after` (every coefficient, including the constants), and the hypothesis
#' that the two sets are equal is tested by the likelihood ratio, as Hausman
#' and Ruud (1987, section 3) do, and by the Wald statistic, as Allison and
#' Christakis (1994, section 5.4) do; both are chi-square with as many degrees
#' of freedom as the model has parameters. A coefficient that one block of
#' stages does not identify (an alternative that no ranking orders there, such
#' as one always ranked within the first `after` positions) is left out of
#' both statistics and their degrees of freedom, with a warning; a penalized
#' fit (`firth = TRUE` or a prior) is then refused, since its fit with
#' separate early and late coefficients has no finite penalty. Chapman and Staelin's (1982)
#' pooling test is the sequential version that compares the first \eqn{E}
#' stages with stage \eqn{E + 1}. Within a tie group the stages take the
#' positions of each order in the sum over orders.
#'
#' A rejection says that later rankings follow different coefficients. That
#' can mean careless ranking further down (the coefficients then shrink
#' towards zero in the late stages, the reading of Chapman and Staelin and of
#' Hausman and Ruud) or a change in valuation; Allison and Christakis also note
#' that apparent differences between stages can arise from heterogeneity
#' across respondents rather than from a change within respondents. Fitting
#' with `depth = after` uses only the early stages.
#'
#' The likelihood-ratio statistic assumes independent rankings (a fit's
#' `cluster` does not enter it); the Wald statistic uses the fit's covariance,
#' cluster-robust when the fit has one. Their chi-square reference is asymptotic; with few rankings
#' per coefficient it can be far off, and no parametric bootstrap of this test
#' is offered. The split design (separate coefficients after the first stage)
#' is outside [separation()]'s scope: when a split fit's data are separated,
#' its maximum-likelihood statistic uses the value where its optimizer
#' stopped, not the supremum (the penalized form is unaffected).
#'
#' @param object A [rologit()] fit.
#' @param after The last position of the early stages (default 1: the first
#'   choice against the rest of the ranking).
#' @return A data frame of class `"iia_test"` with the LR and Wald rows
#'   (`statistic`, `df`, `p.value`); the attribute `"coefficients"` holds the
#'   early and late estimates and their difference with its standard error.
#'   The Wald statistic uses the fit's covariance type (analytic for a
#'   bootstrap fit).
#' @references
#' Allison, P.D. and Christakis, N.A. (1994). Logit models for sets of ranked
#' items. \emph{Sociological Methodology}, 24, 199-228. \doi{10.2307/270983}
#'
#' Chapman, R.G. and Staelin, R. (1982). Exploiting rank ordered choice set
#' data within the stochastic utility model. \emph{Journal of Marketing
#' Research}, 19, 288-301. \doi{10.1177/002224378201900302}
#'
#' Hausman, J.A. and Ruud, P.A. (1987). Specifying and testing econometric
#' models for rank-ordered data. \emph{Journal of Econometrics}, 34, 83-104.
#' \doi{10.1016/0304-4076(87)90068-6}
#' @examples
#' set.seed(5)
#' d <- rrologit(800, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)))
#' rank_test(rologit(rank ~ 1 | x1, data = d, id = "id", alt = "alt"))
#' @seealso [rologit()]
#' @family specification tests
#' @export
rank_test <- function(object, after = 1L) {
  if (!inherits(object, "rologit")) stop("rank_test() takes a rologit() fit.", call. = FALSE)
  after <- .pc_check_whole(after, "after", 1)
  D <- object$design
  ## the position of the last stage of every ranking: the end of the last tie
  ## group that has an alternative (ranked or unranked) below it
  last_stage <- vapply(seq_len(D$n), function(i) {
    r <- D$str$rank[i, D$avail[i, ] == 1L]; ranked <- sort(r[r > 0L])
    if (!length(ranked)) return(0L)
    tb <- tabulate(ranked); start <- cumsum(c(0L, tb))[seq_along(tb)]
    has_below <- seq_along(tb) < length(tb) | any(r == 0L)
    if (!any(has_below)) 0L else as.integer(max((start + tb)[has_below]))
  }, integer(1))
  n_late <- sum(D$w[last_stage > after] > 0)
  if (n_late == 0L) stop("No ranking has a stage after position ", after, "; nothing to compare.", call. = FALSE)
  K0 <- D$K
  D2 <- D; D2$str$split <- as.integer(after)
  fn2 <- .pc_fns(D2)
  sc2 <- rep(object$scale, 2L)
  ## coefficients a block of stages does not identify (an alternative that no
  ## ranking orders within it): no information in that block's column. They
  ## leave the test (local review 2026-11-06, finding 5): its degrees of
  ## freedom count the identified differences only. A penalized fit with
  ## separate early and late coefficients has no finite penalty there, so a
  ## penalized fit is refused.
  I2 <- tryCatch(.pc_firth_fns(.pc_fns_family(D2), D2)$info(c(object$theta, object$theta)), error = function(e) { .pc_pass_limit(e); NULL })
  dI <- if (is.null(I2)) rep(1, 2L * K0) else diag(I2)
  unid_block <- !(dI > 1e-10 * max(dI))
  unid <- unid_block[seq_len(K0)] | unid_block[K0 + seq_len(K0)]
  if (any(unid)) {
    blocks <- c(if (any(unid_block[seq_len(K0)])) "early" , if (any(unid_block[K0 + seq_len(K0)])) "late")
    what <- paste0("the ", paste(blocks, collapse = " and "), " stages (", if (after == 1L) "position 1" else paste0("positions 1-", after), " against later) do not identify ",
                   paste(D$par_names[unid], collapse = ", "), " (no ranking orders those alternatives there)")
    if (.pc_penalized(object))
      stop("rank_test() of a penalized fit needs every coefficient identified in both blocks of stages, since the penalized fit with separate early and late ",
           "coefficients has no finite penalty otherwise; ", what, ". Use another 'after', or the maximum-likelihood fit.", call. = FALSE)
    if (all(unid)) stop("No coefficient is identified in both blocks of stages: ", what, ".", call. = FALSE)
    warning(toupper(substring(what, 1, 1)), substring(what, 2), "; they are left out of the test, which has ", sum(!unid), " degree(s) of freedom.", call. = FALSE)
  }
  op <- .pc_optimize(fn2, c(object$theta, object$theta), sc2, newton_first = .pc_concave(D2))
  ## penalized: a second run from the split design's maximum-likelihood
  ## estimate, the higher maximum kept (.pc_firth_second)
  if (.pc_penalized(object)) {
    f2 <- .pc_firth_second(D2, fn2, sc2, 1000L, 1e-10, op); op <- f2$op
    if (!is.null(f2$modes))
      warning("The penalized likelihood of the fit with separate early and late coefficients has two local maxima (penalized log-likelihoods ",
              formatC(f2$modes[1L], digits = 6, format = "f"), " from the fit's estimates and ", formatC(f2$modes[2L], digits = 6, format = "f"),
              " from the maximum-likelihood estimate); the higher is used.", call. = FALSE)
  }
  if (!isTRUE(op$converged)) warning("The fit with separate early and late coefficients did not converge (largest scaled gradient ", format(op$max_grad, digits = 3),
                                     "); there may be too few late stages for every coefficient.", call. = FALSE)
  lr <- if (.pc_penalized(object)) {
    ## penalized fit: the penalized likelihood ratio, both maxima under the split
    ## model's penalty (early equal to late imposed on it; Heinze and Schemper 2001)
    rr <- .pc_firth_restricted(fn2, sc2, cbind(-diag(K0), diag(K0)), c(object$theta, object$theta))
    if (!isTRUE(rr$converged)) warning("The restricted penalized fit (early equal to late) did not converge; the penalized likelihood-ratio statistic may be inaccurate.", call. = FALSE)
    2 * (-op$value + rr$value)
  } else 2 * (-op$value - object$loglik)
  type <- if (object$se %in% c("analytic", "opg", "robust", "cluster")) object$se else "analytic"
  obs <- .pc_obs(D2)
  V2 <- .pc_vcov(op$par, fn2, type, obs$cluster, sc2, obs$w)
  if (is.null(V2)) stop("The covariance of the fit with separate early and late coefficients could not be computed.", call. = FALSE)
  if (isTRUE(attr(V2, "ginv")) && !any(unid))
    warning("The information matrix of the fit with separate early and late coefficients is singular (a generalized inverse was used): ",
            "the late stages do not identify every coefficient.", call. = FALSE)
  V2 <- .pc_strip_vcov(V2)
  ## the identified differences only (the information is block-diagonal in an
  ## unidentified coefficient's zero row, so the rest of V2 is exact)
  Rm <- cbind(-diag(K0), diag(K0))[!unid, , drop = FALSE]
  df <- sum(!unid)
  delta <- as.vector(Rm %*% op$par)
  Vd <- Rm %*% V2 %*% t(Rm)
  wald <- .pc_wald(delta, Vd)
  if (is.na(wald)) warning("The covariance of the early-minus-late differences is singular (a coefficient at a boundary in one of the fits), so the Wald ",
                           "statistic is NA; the likelihood-ratio row does not need it.", call. = FALSE)
  out <- data.frame(test = c("Likelihood ratio", "Wald"), statistic = c(max(lr, 0), wald), df = df,
                    p.value = stats::pchisq(c(max(lr, 0), wald), df, lower.tail = FALSE), stringsAsFactors = FALSE)
  attr(out, "method") <- paste0("Stability of preferences across ranking stages (", if (after == 1L) "position 1" else paste0("positions 1-", after), " against later; ", n_late, " ranking(s) with later stages",
                                if (.pc_penalized(object)) paste0("; penalized (", .pc_penalty_label(object), ") fits, penalized likelihood ratio") else "", ")")
  dfull <- sefull <- rep(NA_real_, K0); dfull[!unid] <- delta; sefull[!unid] <- sqrt(pmax(diag(Vd), 0))
  attr(out, "coefficients") <- data.frame(early = op$par[seq_len(K0)], late = ifelse(unid_block[K0 + seq_len(K0)], NA_real_, op$par[K0 + seq_len(K0)]),
                                          difference = dfull, se_difference = sefull, row.names = D$par_names)
  class(out) <- c("iia_test", "data.frame")
  out
}
