#' @title Power Calculations for Test of Two Poisson Ratios
#' @description  Compute the power for a test of two sample means with Poisson distributions, or determine parameters to obtain a target power.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, \code{lambda1}, \code{lambda2}, \code{power}, and \code{sig.level} must be passed as NULL, and that parameter is determined from the others.
#' Notice that \code{sig.level} has non-NULL defaults, so NULL must be explicitly passed if you want to compute them.\cr\cr
#' If \code{equal.sample = TRUE} is used, \code{n2} would be ignored and N in output denotes the number in each group.
#' @usage power_Poisson(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
#' lambda1 = NULL, lambda2 = NULL, t1 = 1, t2 = 1, RR0 = 1,
#' equal.sample = TRUE, alternative = c("two.sided", "one.sided"))
#' @param n1 sample size in group 1, or sample size in each group if \code{equal.sample = TRUE}
#' @param n2 sample size in group 2 
#' @param power power of test (1 minus Type II error probability)
#' @param sig.level significance level (Type I error probability)
#' @param lambda1 Poisson rate for group 1
#' @param lambda2 Poisson rate for group 2
#' @param t1 observed time period for group 1
#' @param t2 observed time period for group 2
#' @param RR0 the ratio of lambda2 and lambda1 under null hypothesis
#' @param equal.sample equal sample sizes for two groups, see details
#' @param alternative one- or two-sided test
#' @references Gu et al. (2008). Testing the ratio of two poisson rates. \emph{Biometrical Journal: Journal of Mathematical Methods in Biosciences}. \bold{50}:283-298.
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with method element.
#' @examples 
#' power_Poisson(lambda1 = 0.0005, lambda2 = 0.003, n1 = 2000, t1 = 2, t2 = 2)
#' power_Poisson(lambda1 = 0.0005, lambda2 = 0.003, power = 0.8, t1 = 2, t2 = 2)
#' power_Poisson(n1 = 2000, lambda1 = 0.0005, lambda2 = 0.003, power = 0.8, 
#' t1 = 2, t2 = 2, equal.sample = FALSE)
#' @note 'uniroot' is used to solve power equation for unknowns, 
#' so you may see errors from it, notably about inability to bracket the root when invalid arguments are given.
#' @importFrom stats uniroot pnorm
#' @export
power_Poisson <- function(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
                          lambda1 = NULL, lambda2 = NULL, t1 = 1, t2 = 1, RR0 = 1,
                          equal.sample = TRUE, alternative = c("two.sided", "one.sided")){
  if(equal.sample==TRUE){
    if (sum(sapply(list(n1, lambda1, lambda2, power, sig.level), is.null)) != 1) 
      stop("exactly one of n1, lambda1, lambda2, power, and sig.level must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
      stop(sQuote("sig.level"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1)) 
      stop(sQuote("power"), " must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    test.side <- switch(alternative,one.sided = 1, two.sided = 2)
    p.body <-quote({
      RR_a <- lambda2/lambda1
      para_d <- (t1)/(t2)
      para_A <- 2*(1-sqrt(RR0/RR_a))
      para_B <- lambda1*t1*n1 + 3/8
      para_C <- sqrt((RR0+para_d)/RR_a)
      para_D <- sqrt((RR_a+para_d)/RR_a)
      pnorm((abs(para_A)*sqrt(para_B)-qnorm(sig.level/test.side,lower.tail = FALSE)*para_C)/para_D)
    })
    if (is.null(power)) 
      power <- eval(p.body)
    else if(is.null(n1))
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    else if(is.null(lambda1))
      lambda1 <- uniroot(function(lambda1) eval(p.body) - power, c(1e-07, 1e7))$root
    else if(is.null(lambda2))
      lambda2 <- uniroot(function(lambda2) eval(p.body) - power, c(1e-07, 1e7))$root
    else if(is.null(sig.level))
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1-1e-10))$root
    
    METHOD <- "Two-sample Poisson Ratio Tests (Equal Sizes)"
    NOTE <- "N is number in *each* group"
    return(structure(list(N = n1, lambda1 = lambda1, lambda2 = lambda2, sig.level = sig.level, power = power, 
                          alternative = alternative,  method = METHOD, note = NOTE), 
                     class = "power.htest"))
  }
  else{
    if (sum(sapply(list(n1, n2, lambda1, lambda2, power, sig.level), is.null)) != 1) 
      stop("exactly one of n1, n2, lambda1, lambda2, power, and sig.level must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)) 
      stop(sQuote("sig.level"), " must be numeric in [0, 1]")
    if (!is.null(power) && !is.numeric(power) || any(0 > power | power > 1)) 
      stop(sQuote("power"), " must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    test.side <- switch(alternative,one.sided = 1, two.sided = 2)
    p.body <-quote({
      RR_a <- lambda2/lambda1
      para_d <- (t1*n1)/(t2*n2)
      para_A <- 2*(1-sqrt(RR0/RR_a))
      para_B <- lambda1*t1*n1 + 3/8
      para_C <- sqrt((RR0+para_d)/RR_a)
      para_D <- sqrt((RR_a+para_d)/RR_a)
      pnorm((abs(para_A)*sqrt(para_B)-qnorm(sig.level/test.side,lower.tail = FALSE)*para_C)/para_D)
    })
    if (is.null(power)) 
      power <- eval(p.body)
    else if(is.null(n1))
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    else if(is.null(n2))
      n2 <- uniroot(function(n2) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
    else if(is.null(lambda1))
      lambda1 <- uniroot(function(lambda1) eval(p.body) - power, c(1e-7, 1e7))$root
    else if(is.null(lambda2))
      lambda2 <- uniroot(function(lambda2) eval(p.body) - power, c(1e-7, 1e7))$root
    else if(is.null(sig.level))
      sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1-1e-10))$root
    
    METHOD <- "Two-sample Poisson Ratio Tests (Different Sizes)"
    return(structure(list(n1 = n1, n2 = n2, lambda1 = lambda1, lambda2 = lambda2, sig.level = sig.level, power = power, 
                          alternative = alternative,  method = METHOD), 
                     class = "power.htest"))
  }
}