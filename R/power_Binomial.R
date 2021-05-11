#' @title Power Calculations for Two-Sample Test for Proportions
#' @description  Compute power of test, or determine parameters to obtain target power for equal and unequal sample sizes.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, \code{p1}, \code{p2}, \code{power}, and \code{sig.level} must be passed as NULL, and that parameter is determined from the others.
#' Notice that \code{p1}, \code{p2}, \code{sig.level} have non-NULL defaults, so NULL must be explicitly expressed if you want to compute them.\cr\cr
#' If \code{equal.sample = TRUE} is used, N in output will denote the number in each group.\cr\cr
#' @usage power_Binomial(n1 = NULL, n2 = NULL, p1 = 0.5, p2 = 0.5, 
#' sig.level = 0.05, power = NULL, equal.sample = TRUE,
#' alternative = c("two.sided", "one.sided"))
#' @param n1 sample size in group 1, or sample size in each group if equal.sample = TRUE
#' @param n2 sample size in group 2
#' @param p1 probability in group 1
#' @param p2 probability in group 2
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param equal.sample equal sample sizes for two groups, see details
#' @param alternative one- or two-sided test
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with note and method elements.
#' @examples
#' # calculate power, equal sizes
#' power_Binomial(n1 = 100, p1 = 0.5, p2 = 0.7)
#' # calculate power, unequal sizes
#' power_Binomial(n1 = 150, n2 = 100, p1 = 0.5, p2 = 0.7)
#' # calculate n2
#' power_Binomial(n1 = 100, p1 = 0.5, p2 = 0.7, power = 0.9, equal.sample = FALSE)
#' @export
power_Binomial <- function(n1 = NULL, n2 = NULL, p1 = 0.5, p2 = 0.5, 
                           sig.level = 0.05, power = NULL, equal.sample = TRUE,
                           alternative = c("two.sided", "one.sided")) {
  # define the inputs
  alternative <- match.arg(alternative)
  ratio <- NULL
  if (!is.null(n1) & !is.null(n2)) {
    if (n1 != n2) {
      equal.sample <- FALSE
      ratio <- n2 / n1
    }
    else {
      equal.sample <- TRUE
    }
  }
  ## exactly one option must be empty
  if (equal.sample) {
    ratio <- 1
    if (sum(sapply(list(n1, p1, p2, power, sig.level), is.null)) != 1) {
      stop("exactly one of n1, p1, p2, power, sig.level must be NULL")
    }
  }
  else {
    if (sum(sapply(list(n1, n2, p1, p2, power, sig.level), is.null)) != 1) {
      stop("exactly one of n1, n2, p1, p2, power, sig.level must be NULL")
    }
  }
  ## significant level limit
  if (!is.null(sig.level) && (!is.numeric(sig.level) || (sig.level < 0 | sig.level > 1))) {
    stop("sig.level must be numeric in [0,1]")
  }
  ## n1 limit
  if (!is.null(n1) && (!is.numeric(n1) || n1 <= 0)) {
    stop("n1 must be positive number")
  }
  ## n2 limit
  if (!equal.sample && (!is.null(n2) && (!is.numeric(n2) || n2 <= 0))) {
    stop("n2 must be positive number")
  }
  ## p1 limit
  if (!is.null(p1) && (!is.numeric(p1) || (p1 < 0 | p1 > 1))) {
    stop("p1 must be numeric in [0,1]")
  }
  ## p2 limit
  if (!is.null(p2) && (!is.numeric(p2) || (p2 < 0 | p2 > 1))) {
    stop("p2 must be numeric in [0,1]")
  }
  ## power limit
  if (!is.null(power) && (!is.numeric(power) || (power < 0 | power > 1))) {
    stop("power must be numeric in [0,1]")
  }
  ## define test sided
  tside <- as.numeric(alternative == "two.sided") + 1

  # Calculate the outputs
  ## Equal sample
  if (equal.sample) {
    ### algorithm
    p.body <- quote({
      z_alpha <- -qnorm(sig.level / tside)
      d <- abs(p1 - p2)
      q1 <- 1 - p1
      q2 <- 1 - p2
      pbar <- (p1 + p2) / 2
      qbar <- 1 - pbar
      pnorm( (d * sqrt(n1) - z_alpha * sqrt(2 * pbar * qbar)) / sqrt(p1 * q1 + p2 * q2) )
    })
    ### calculation
    if (is.null(power)) {
      power <- eval(p.body)
    } 
    else if (is.null(p1)) {
      p1 <- uniroot(function(p1) eval(p.body) - power, c(0, p2), extendInt = "yes")$root
    } 
    else if (is.null(p2)) {
      p2 <- uniroot(function(p2) eval(p.body) - power, c(p1, 1), extendInt = "yes")$root
    } 
    else if (is.null(n1)) {
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2, 1e7))$root
    } 
    else if (is.null(sig.level)) {
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10), extendInt = "upX")$root
    }
    METHOD <- "Two-sample comparison of proportions power calculation"
    NOTE <- "N is number in *each* group"
    Output <- structure(list(
      N = n1, p1 = p1, p2 = p2,
      sig.level = sig.level, power = power, alternative = alternative,
      method = METHOD, note = NOTE
    ), class = "power.htest")
  }
  ## Unequal sample sizes
  else {
    ### algorithm
    p.body <- quote({
      z_alpha <- -qnorm(sig.level / tside)
      d <- abs(p1 - p2)
      q1 <- 1 - p1
      q2 <- 1 - p2
      pbar <- (p1 + ratio * p2) / (1 + ratio)
      qbar <- 1 - pbar
      pnorm( (d * sqrt(ratio * n1) - z_alpha * sqrt((1 + ratio) * pbar * qbar)) / sqrt(ratio * p1 * q1 + p2 * q2) )
    })
    ### the algorithm above may have issue when calculting n2, p1 & p2 using uniroot() function
    p.body.backup <- quote({
      z_alpha <- -qnorm(sig.level / tside)
      d <- abs(p1 - p2)
      q1 <- 1 - p1
      q2 <- 1 - p2
      pbar <- (p1 + ratio * p2) / (1 + ratio)
      qbar <- 1 - pbar
      ( - z_alpha * sqrt((1 + ratio) * pbar * qbar) + qnorm(1 - power) * sqrt(ratio * p1 * q1 + p2 * q2) ) ^ 2 / (ratio * d ^ 2)
    })
    ### calculation
    if (is.null(n2)) {
      ratio <- tryCatch(uniroot(function(ratio) eval(p.body) - power, c(2 / n1, 1e7))$root, 
                        error=function(err) NA)
      if(is.na(ratio))
        ratio <- uniroot(function(ratio) eval(p.body.backup) - power, c(2 / n1, 1e7))$root
      n2 <- n1 * ratio
    }
    else if (is.null(n1)) {
      ratio <- uniroot(function(ratio) eval(p.body) - power, c(2 / n2, 1e7))$root
      n1 <- n2 / ratio
    }
    else if (is.null(p1)) {
      p1 <- tryCatch(uniroot(function(p1) eval(p.body) - power, c(0, p2), extendInt = "yes")$root,
                     error=function(err) NA)
      if(is.na(p1))
        p1 <- uniroot(function(p1) eval(p.body.backup) - n1, c(0, p2), extendInt = "yes")$root
    } 
    else if (is.null(p2)) {
      p2 <- tryCatch(uniroot(function(p2) eval(p.body) - power, c(p1, 1), extendInt = "yes")$root,
                     error=function(err) NA)
      if(is.na(p2))
        p2 <- uniroot(function(p2) eval(p.body.backup) - n1, c(p1, 1), extendInt = "yes")$root
    } 
    else if (is.null(power)) {
      power <- eval(p.body)
    } 
    else if (is.null(sig.level)) {
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10), extendInt = "upX")$root
    }
    METHOD <- "Two-sample comparison of proportions power calculation with unequal sample sizes"
    Output <- structure(list(
      n1 = n1, n2 = n2, p1 = p1, p2 = p2,
      sig.level = sig.level, power = power, alternative = alternative,
      method = METHOD
    ), class = "power.htest")
  }
  return(Output)
}
