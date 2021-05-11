#' @title Power Calculations for One and Two Sample T-tests
#' @description  Compute power of t test, or determine parameters to obtain target power.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, \code{delta}, \code{sd1}, \code{sd2}, \code{power}, and \code{sig.level} must be passed as NULL, and that parameter is determined from the others.
#' Notice that \code{sd1}, \code{sd2}, \code{sig.level} have non-NULL defaults, so NULL must be explicitly expressed if you want to compute them.\cr\cr
#' If \code{equal.sample = TRUE} is used, N in output will denote the number in each group.\cr\cr
#' @usage power_Normal(n1 = NULL, n2 = NULL, delta = NULL, sd1 = 1, sd2 = 1, 
#' sig.level = 0.05, power = NULL, equal.sample = TRUE,
#' alternative = c("two.sided", "one.sided"),
#' type = c("two.sample", "one.sample", "paired"),
#' df.method = c("welch", "classical"), strict = FALSE)
#' @param n1 sample size in group 1, or sample size in each group if equal.sample = TRUE
#' @param n2 sample size in group 2
#' @param delta true difference in means
#' @param sd1 standard deviation for group 1
#' @param sd2 standard deviation for group 2
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param equal.sample equal sample sizes for two groups, see details
#' @param alternative one- or two-sided test
#' @param type Type of t test
#' @param df.method Method for calculating the degrees of default. Possibilities are welch (the default) or classical.
#' @param strict Use strict interpretation in two-sided case
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with note and method elements.
#' @examples
#' # Calculate power, equal sizes
#' power_Normal(n1 = 150, delta = 5, sd1 = 20, sd2 = 10)
#' # Calculate power, unequal sizes
#' power_Normal(n1 = 150, delta = 5, n2 = 120, sd1 = 10)
#' # Calculate n1, equal sizes
#' power_Normal(delta = 5, power = 0.9, sd1 = 10, sd2 = 12)
#' @note 'uniroot' is used to solve power equation for unknowns,
#' so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
#' @export
power_Normal <- function(n1 = NULL, n2 = NULL, delta = NULL, sd1 = 1, sd2 = 1, 
                         sig.level = 0.05, power = NULL, equal.sample = TRUE,
                         alternative = c("two.sided", "one.sided"),
                         type = c("two.sample", "one.sample", "paired"),
                         df.method = c("welch", "classical"), strict = FALSE) {
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  df.method <- match.arg(df.method)
  # apply the option "equal.sample"
  if (!is.null(n1) & !is.null(n2)) {
    if (n1 != n2) {
      equal.sample <- FALSE
    }
    else {
      equal.sample <- TRUE
    }
  }
  # define the inputs
  ## exactly one option must be empty
  if (equal.sample) {
    ratio <- 1
    if (sum(sapply(list(n1, delta, sd1, sd2, power, sig.level), is.null)) != 1) {
      stop("exactly one of n1, delta, sd1, sd2, power, sig.level must be NULL")
    }
  }
  else {
    if (type == "one.sample") {
      warning("equal.sample option is not available for one sample t test")
      type <- "two.sample"
    }
    if (type == "paired") {
      warning("equal.sample option is not available for paired t test")
      type <- "two.sample"
    }
    if (sum(sapply(list(n1, n2, delta, sd1, sd2, power, sig.level), is.null)) != 1) {
      stop("exactly one of n1, n2, delta, sd1, sd2, power, sig.level must be NULL")
    }
    if (!is.null(n1) & !is.null(n2)) {
      ratio <- n2 / n1
    }
    else {
      ratio <- NULL
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
  ## sd1 limit
  if (!is.null(sd1) && (!is.numeric(sd1) || sd1 <= 0)) {
    stop("sd1 must be positive number")
  }
  ## sd2 limit
  if (!is.null(sd2) && (!is.numeric(sd2) || sd2 <= 0)) {
    stop("sd2 must be positive number")
  }
  ## power limit
  if (!is.null(power) && (!is.numeric(power) || (power < 0 | power > 1))) {
    stop("power must be numeric in [0,1]")
  }

  tside <- as.numeric(alternative == "two.sided") + 1
  ttype <- switch(type,
    one.sample = 1,
    two.sample = 2,
    paired = 1
  )
  if (tside == 2 && !is.null(delta)) {
    delta <- abs(delta)
  }

  # algorithm
  p.body <- quote({
    nu <- switch(ttype,
      n1 - 1,
      switch(df.method,
        welch = (sd1^2 / n1 + sd2^2 / (n1 * ratio))^2 / ((sd1^2 / n1)^2 / (n1 - 1) + (sd2^2 / (n1 * ratio))^2 / (n1 * ratio - 1)),
        classical = n1 * (1 + ratio) - 2
      )
    )
    qu <- -qt(sig.level / tside, nu)
    ncp <- delta / sd1 * switch(ttype,
      sqrt(n1 / ttype),
      sqrt(n1 / (1 + (sd2 / sd1)^2 / ratio))
    )
    pt(qu, nu, ncp = ncp, lower.tail = FALSE)
  })
  ## if strict = TRUE
  if (strict == TRUE & tside == 2) {
    p.body <- quote({
      nu <- switch(ttype,
        n1 - 1,
        switch(df.method,
          welch = (sd1^2 / n1 + sd2^2 / (n1 * ratio))^2 / ((sd1^2 / n1)^2 / (n1 - 1) + (sd2^2 / (n1 * ratio))^2 / (n1 * ratio - 1)),
          classical = n1 * (1 + ratio) - 2
        )
      )
      qu <- -qt(sig.level / tside, nu)
      ncp <- delta / sd1 * switch(ttype,
        sqrt(n1 / ttype),
        sqrt(n1 / (1 + (sd2 / sd1)^2 / ratio))
      )
      pt(qu, nu, ncp = ncp, lower.tail = FALSE) + pt(-qu, nu, ncp = ncp, lower.tail = TRUE)
    })
  }

  # Calculate the outputs
  if (equal.sample) {
    if (is.null(power)) {
      power <- eval(p.body)
    } else if (is.null(n1)) {
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2, 1e7))$root
    } else if (is.null(delta)) {
      delta <- uniroot(function(delta) eval(p.body) - power, sd1 * c(1e-7, 1e7))$root
    } else if (is.null(sd1)) {
      sd1 <- uniroot(function(sd1) eval(p.body) - power, delta * c(1e-7, 1e7))$root
    } else if (is.null(sd2)) {
      sd2 <- uniroot(function(sd2) eval(p.body) - power, delta * c(1e-7, 1e7))$root
    } else if (is.null(sig.level)) {
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10), extendInt = "upX")$root
    }
    METHOD <- switch(type,
      one.sample = "One-sample t test power calculation",
      two.sample = "Two-sample t test power calculation",
      paired = "Paired t test power calculation"
    )
    NOTE <- switch(type,
      one.sample = NULL,
      two.sample = "N is number in *each* group",
      paired = "N is number of *pairs*, sd is std dev of *differences* within pairs"
    )
    Output <- structure(list(
      N = n1, delta = delta, sd1 = sd1, sd2 = sd2,
      sig.level = sig.level, power = power, alternative = alternative,
      method = METHOD, note = NOTE
    ), class = "power.htest")
  }
  else {
    if (is.null(power)) {
      power <- eval(p.body)
    } else if (is.null(n1)) {
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2, 1e7))$root
    } else if (is.null(n2)) {
      ratio <- uniroot(function(ratio) eval(p.body) - power, c(2 / n1, 1e7))$root
      n2 <- n1 * ratio
    }
    else if (is.null(delta)) {
      delta <- uniroot(function(delta) eval(p.body) - power, sd1 * c(1e-7, 1e7))$root
    } else if (is.null(sd1)) {
      sd1 <- uniroot(function(sd1) eval(p.body) - power, delta * c(1e-7, 1e7))$root
    } else if (is.null(sd2)) {
      sd2 <- uniroot(function(sd2) eval(p.body) - power, delta * c(1e-7, 1e7))$root
    } else if (is.null(sig.level)) {
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10), extendInt = "upX")$root
    }
    METHOD <- "Two-sample t test power calculation with unequal sample sizes"
    Output <- structure(list(
      n1 = n1, n2 = n2, delta = delta, sd1 = sd1, sd2 = sd2,
      sig.level = sig.level, power = power, alternative = alternative,
      method = METHOD
    ), class = "power.htest")
  }
  return(Output)
}
