#' @title Power Calculation for Comparing Two Negative Binomial Rates
#' @description  Compute sample size or power for comparing two negative binomial rates.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, and \code{power} must be passed as NULL, and that parameter is determined from the others.\cr\cr
#' If \code{equal.sample = TRUE} is used, N in output will denote the number in each group.\cr\cr
#' The computations are based on the formulas given in Zhu and Lakkis (2014). See \link[MKmisc]{power.nb.test} for more details.
#' @usage power_NegativeBinomial(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
#' mu1 = NULL, mu2 = NULL, duration = 1, theta = NULL, equal.sample = TRUE, 
#' alternative = c("two.sided", "one.sided"), approach = 3)
#' @param n1 sample size in group 1, or sample size in each group if \code{equal.sample = TRUE}
#' @param n2 sample size in group 2 
#' @param power power of test (1 minus Type II error probability)
#' @param sig.level significance level (Type I error probability)
#' @param mu1 expected rate of events per time unit for group 1
#' @param mu2 expected rate of events per time unit for group 2
#' @param duration (average) treatment duration
#' @param theta theta parameter of negative binomial distribution; see \link[MASS]{rnegbin}
#' @param equal.sample equal sample sizes for two groups, see details
#' @param alternative one- or two-sided test
#' @param approach 1, 2, or 3; see Zhu and Lakkis (2014).
#' @references H. Zhu and H. Lakkis (2014). Sample size calculation for comparing two negative binomial rates. \emph{Statistics in Medicine}, \bold{33}:376-387.
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with note and method elements.
#' @examples 
#' # calculate power, equal sizes
#' power_NegativeBinomial(n1 = 20, mu1 = 1, mu2 = 2, theta = 0.8)
#' # calculate power, unequal sizes
#' power_NegativeBinomial(n1 = 80, n2 = 40, mu1 = 1, mu2 = 2, theta = 0.8)
#' # calculate n
#' power_NegativeBinomial( mu1 = 1, mu2 = 2, theta = 0.8, power = 0.8)
#' @importFrom stats qnorm pnorm uniroot 
#' @export
power_NegativeBinomial <- function(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05,
                                   mu1 = NULL, mu2 = NULL, duration = 1, theta = NULL,
                                   equal.sample = TRUE, alternative = c("two.sided", "one.sided"),
                                   approach = 3){
  if(!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in (0, 1)")
  if(!is.null(power) && (!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1)))
    stop("'power' must be numeric in (0, 1)")
  if(missing(mu1))
    stop("'mu1' must be given")
  if(missing(mu2))
    stop("'mu2' must be given")
  if(missing(theta))
    stop("'theta' must be given")
  if(!is.numeric(mu1) || mu1 <= 0)
    stop("'mu1' must be numeric and > 0")
  if(!is.numeric(mu2) || mu2 <= 0)
    stop("'mu2' must be numeric and > 0")
  if(!is.numeric(duration) || duration <= 0)
    stop("'duration' must be numeric and > 0")
  if(!is.numeric(theta) || theta <= 0)
    stop("'theta' must be numeric and > 0")
  if(!(approach %in% 1:3))
    stop("'approach' must be equal to 1, 2, or 3")
  alternative <- match.arg(alternative)
  sided <- switch(alternative, one.sided = 1, two.sided = 2)
  new.sig.level <- sig.level/sided
  if(!is.null(n1)&!is.null(n2)){
    if(n1 != n2){
      equal.sample <- FALSE
    }
    else{
      equal.sample <- TRUE
    }
  }
  if(equal.sample){
    METHOD <- "Two-sample Negative Binomial rates Tests (Equal Sizes)"
    n2 <- n1
    if (sum(sapply(list(n1, power), is.null)) != 1) 
      stop("exactly one of n1 and power must be NULL")
    if(approach == 1){
      p.body <- quote({
        V0 <- 2/(duration*mu1) + 2/theta
        V1 <- 1/duration*(1/mu1 + 1/mu2) + 2/theta
        pnorm((sqrt(n1)*abs(log(mu2/mu1)) - qnorm(1-new.sig.level)*sqrt(V0))/sqrt(V1))
      })
    }
    if(approach == 2){
      p.body <- quote({
        V0 <- 1/duration*(1/mu1 + 1/mu2) + 2/theta
        pnorm((sqrt(n1)*abs(log(mu2/mu1)) - qnorm(1-new.sig.level)*sqrt(V0))/sqrt(V0))
      })
    }
    if(approach == 3){
      p.body <- quote({
        V0 <- 2^2/(duration*(mu1 + mu2)) + 2/theta
        V1 <- 1/duration*(1/mu1 + 1/mu2) + 2/theta
        pnorm((sqrt(n1)*abs(log(mu2/mu1)) - qnorm(1-new.sig.level)*sqrt(V0))/sqrt(V1))
      })
    }
    if(is.null(power)){
      power <- eval(p.body)
    }
    if(is.null(n1)){
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2+1e-7, 1e7))$root
    }
    NOTE <- "N is number in *each* group"
    Output <- structure(list(N = n1, mu1 =  mu1, mu2 =  mu2,
                             theta =  theta, duration =  duration, sig.level =  sig.level, 
                             power =  power, alternative =  alternative,
                             method = METHOD, note = NOTE), class = "power.htest")
  }
  else{
    METHOD <- "Two-sample Negative Binomial rates Tests (Different Sizes)"
    if (sum(sapply(list(n1, n2, power), is.null)) != 1) 
      stop("exactly one of n1, n2, and power must be NULL")
    if(approach == 1){
      p.body <- quote({
        ratio <- n2/n1
        V0 <- (1 + ratio)/(duration*ratio*mu1) + (1 + ratio)/theta/ratio
        V1 <- 1/duration*(1/mu1 + 1/(ratio*mu2)) + (1 + ratio)/theta/ratio
        pnorm((sqrt(n1)*abs(log(mu2/mu1)) - qnorm(1-new.sig.level)*sqrt(V0))/sqrt(V1))
      })
    }
    if(approach == 2){
      p.body <- quote({
        ratio <- n2/n1
        V0 <- 1/duration*(1/mu1 + 1/(ratio*mu2)) + (1 + ratio)/theta/ratio
        pnorm((sqrt(n1)*abs(log(mu2/mu1)) - qnorm(1-new.sig.level)*sqrt(V0))/sqrt(V0))
      })
    }
    if(approach == 3){
      p.body <- quote({
        ratio <- n2/n1
        V0 <- (1+ratio)^2/(duration*ratio*(mu1 + ratio*mu2)) + (1 + ratio)/theta/ratio
        V1 <- 1/duration*(1/mu1 + 1/(ratio*mu2)) + (1 + ratio)/theta/ratio
        pnorm((sqrt(n1)*abs(log(mu2/mu1)) - qnorm(1-new.sig.level)*sqrt(V0))/sqrt(V1))
      })
    }
    if(is.null(n2)){
      n2 <- uniroot(function(n2) eval(p.body) - power, c(2+1e-7, 1e7))$root
    }
    if(is.null(power)){
      power <- eval(p.body)
    }
    if(is.null(n1)){
      n1 <- uniroot(function(n1) eval(p.body) - power, c(2+1e-7, 1e7))$root
    }
    Output <- structure(list(n1 = n1, n2 =  n2, mu1 =  mu1, mu2 =  mu2,
                             theta =  theta, duration =  duration, sig.level =  sig.level, 
                             power =  power, alternative =  alternative, 
                             method = METHOD), class = "power.htest")
  }
  
  return(Output)
}