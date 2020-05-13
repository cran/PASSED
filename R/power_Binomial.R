#' @title Power Calculations for Two-Sample Test for Proportions
#' @description  Compute power of test, or determine parameters to obtain target power for equal and unequal sample sizes.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, \code{p1}, \code{p2}, \code{power}, and \code{sig.level} must be passed as NULL, and that parameter is determined from the others.
#' Notice that \code{p1}, \code{p2}, \code{sig.level} have non-NULL defaults, so NULL must be explicitly passed if you want to compute them.\cr\cr
#' If \code{equal.sample = TRUE} is used, N in output will denote the number in each group.\cr\cr
#' See \link[MESS]{power_prop_test} for more details.
#' @usage power_Binomial(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
#' p1 = 0.5, p2 = 0.5, equal.sample = TRUE, alternative = c("two.sided", "one.sided"))
#' @param n1 sample size in group 1, or sample size in each group if equal.sample = TRUE
#' @param n2 sample size in group 2 
#' @param power power of test (1 minus Type II error probability)
#' @param sig.level significance level (Type I error probability)
#' @param p1 probability in group 1
#' @param p2 probability in group 2
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
#' @importFrom MESS power_prop_test
#' @export
power_Binomial <- function(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
                           p1 = 0.5, p2 = 0.5, equal.sample = TRUE, 
                           alternative = c("two.sided", "one.sided")){
  alternative <- match.arg(alternative)
  ratio <- NULL
  if(!is.null(n1)&!is.null(n2)){
    if(n1 != n2){
      equal.sample <- FALSE
      ratio <- n2/n1
    }
    else{
      equal.sample <- TRUE
    }
  }
  if(equal.sample){
    ratio <- 1
    if (sum(sapply(list(n1, p1, p2, power, sig.level), is.null)) != 1) 
      stop("exactly one of n1, p1, p2, power, sig.level must be NULL")
  }
  else{
    if (sum(sapply(list(n1, n2, p1, p2, power, sig.level), is.null)) != 1) 
      stop("exactly one of n1, n2, p1, p2, power, sig.level must be NULL")
  }
  Output <- NULL
  if(is.null(n1)){
    if(is.null(ratio)){
      new.ratio <- NULL
    }
    else{
      new.ratio <- 1/ratio
    }
    Output <- do.call(MESS::power_prop_test,list(n = n2, p1 = p2, p2 = p1, sig.level = sig.level,
                                                 power = power, ratio = new.ratio, alternative = alternative))
    Output$n <- rev(Output$n)
    Tmp_p1 <- Output$p1
    Output$p1 <- Output$p2
    Output$p2 <- Tmp_p1
  }
  else{
    if(!is.null(ratio)){
      if(ratio<1){
        Output <- do.call(MESS::power_prop_test,list(n = n2, p1 = p2, p2 = p1, sig.level = sig.level,
                                                     power = power, ratio = 1/ratio, alternative = alternative))
        Output$n <- rev(Output$n)
        Tmp_p1 <- Output$p1
        Output$p1 <- Output$p2
        Output$p2 <- Tmp_p1
      }
    }
    if(is.null(Output)){
      Output <- do.call(MESS::power_prop_test,list(n = n1, p1 = p1, p2 = p2, sig.level = sig.level,
                                                   power = power, ratio = ratio, alternative = alternative))
    }
  }
  return(Output)
}