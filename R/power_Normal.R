#' @title Power Calculations for One and Two Sample T-tests
#' @description  Compute power of t test, or determine parameters to obtain target power.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, \code{delta}, \code{sd1}, \code{sd2}, \code{power}, and \code{sig.level} must be passed as NULL, and that parameter is determined from the others.
#' Notice that \code{sd1}, \code{sd2}, \code{sig.level} have non-NULL defaults, so NULL must be explicitly passed if you want to compute them.\cr\cr
#' If \code{equal.sample = TRUE} is used, N in output will denote the number in each group.\cr\cr
#' See \link[MESS]{power_t_test} for more details.
#' @usage power_Normal(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
#' delta = NULL, sd1 = 1, sd2 = 1, equal.sample = TRUE,
#' alternative = c("two.sided", "one.sided"),
#' type = c("two.sample", "one.sample", "paired"),
#' df.method = c("welch", "classical"), strict = FALSE)
#' @param n1 sample size in group 1, or sample size in each group if equal.sample = TRUE
#' @param n2 sample size in group 2 
#' @param power power of test (1 minus Type II error probability)
#' @param sig.level significance level (Type I error probability)
#' @param delta true difference in means
#' @param sd1 standard deviation for group 1
#' @param sd2 standard deviation for group 2
#' @param equal.sample equal sample sizes for two groups, see details
#' @param alternative one- or two-sided test
#' @param type Type of t test
#' @param df.method Method for calculating the degrees of default. Possibilities are welch (the default) or classical.
#' @param strict Use strict interpretation in two-sided case
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with note and method elements.
#' @examples 
#' # Calculate power, equal sizes
#' power_Normal(n1 = 150, delta = 5, sd1 = 20, sd2 = 10)
#' # calculate power, unequal sizes
#' power_Normal(n1 = 150, delta = 5, n2 = 120, sd1 = 10)
#' # calculate n1, equal sizes 
#' power_Normal(delta = 5,  power = 0.9, sd1 = 10, sd2 = 12)
#' @note 'uniroot' is used to solve power equation for unknowns, 
#' so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
#' @importFrom MESS power_t_test
#' @export
power_Normal <- function(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
                         delta = NULL, sd1 = 1, sd2 = 1, equal.sample = TRUE,
                         alternative = c("two.sided", "one.sided"),
                         type = c("two.sample", "one.sample", "paired"),
                         df.method = c("welch", "classical"), strict = FALSE){
  alternative <- match.arg(alternative)
  type <- match.arg(type)
  df.method <- match.arg(df.method)
  if(!is.null(n1)&!is.null(n2)){
    if(n1 != n2){
      equal.sample <- FALSE
    }
    else{
      equal.sample <- TRUE
    }
  }
  if(equal.sample){
    ratio <- 1
    if (sum(sapply(list(n1, delta, sd1, sd2, power, sig.level), is.null)) != 1) 
      stop("exactly one of n1, delta, sd1, sd2, power, sig.level must be NULL")
  }
  else{
    if (sum(sapply(list(n1, n2, delta, sd1, sd2, power, sig.level), is.null)) != 1) 
      stop("exactly one of n1, n2, delta, sd1, sd2, power, sig.level must be NULL")
    if(!is.null(n1)&!is.null(n2)){
      ratio <- n2/n1
    }
    else{
      ratio <- NULL
    }
  }
  if(!is.null(sd1)&!is.null(sd2)){
    sd.ratio <- sd2/sd1
  }
  else{
    sd.ratio <- NULL
  }
  Output <- NULL
  # Case 1: when n1 is missing, ratio would be NULL; delta would be -delta; and group numbers would be exchanged during calculation
  if(is.null(n1)|is.null(sd1)){
    if(is.null(sd.ratio)){
      new.sd.ratio <- NULL
    }
    else{
      new.sd.ratio <- 1/sd.ratio
    }
    if(is.null(ratio)){
      new.ratio <- NULL
    }
    else{
      new.ratio <- 1/ratio
    }
    # ratio = NULL
    Output <- MESS::power_t_test(n = n2, delta = -delta, sd = sd2, sig.level = sig.level,
                                 power = power, ratio = new.ratio, sd.ratio = new.sd.ratio, type = type,
                                 alternative = alternative, df.method = df.method, strict = strict)
    Output$n <- rev(Output$n)
    Output$sd <- rev(Output$sd)
  }
  else{
    # Case 2: when ratio < 1, group numbers would be exchanged during calculation
    if(!is.null(ratio)){
      if(ratio < 1){
        if(is.null(sd.ratio)){
          new.sd.ratio <- NULL
        }
        else{
          new.sd.ratio <- 1/sd.ratio
        }
        Output <- MESS::power_t_test(n = n2, delta = -delta, sd = sd2, sig.level = sig.level,
                                     power = power, ratio = 1/ratio, sd.ratio = new.sd.ratio, type = type,
                                     alternative = alternative, df.method = df.method, strict = strict)
        Output$n <- rev(Output$n)
        Output$sd <- rev(Output$sd)
      }
    }
    if(is.null(Output)){
      Output <- MESS::power_t_test(n = n1, delta = delta, sd = sd1, sig.level = sig.level,
                                   power = power, ratio = ratio, sd.ratio = sd.ratio, type = type,
                                   alternative = alternative, df.method = df.method, strict = strict)
    }
  }
  return(Output)
}
