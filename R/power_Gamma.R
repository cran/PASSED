#' @title Power Calculations for Test of Two Gamma Means
#' @description  Compute the power for a test of two sample means with Gamma distributions, or determine parameters to obtain a target power.
#' @details Exactly one of the parameters \code{n1}, \code{n2}, and \code{power} must be passed as NULL, and that parameter is determined from the others.
#' Notice that \code{sig.level} has non-NULL defaults, so NULL must be explicitly passed if you want to compute it.\cr\cr
#' If \code{equal.sample = TRUE} is used, N in output will denote the number in each group.\cr\cr
#' The equal shape parameter assumption will be tested automatically; otherwise it could be set manually with \code{equal.shape}.
#' @usage power_Gamma(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05, 
#' mu1 = NULL, mu2 = NULL, gmu1 = NULL, gmu2 = NULL, trials = 100, 
#' M = 10000, equal.sample = TRUE, equal.shape = NULL, trace = FALSE)
#' @param n1 sample size in group 1, or sample size in each group if equal.sample = TRUE
#' @param n2 sample size in group 2
#' @param power power of test (1 minus Type II error probability)
#' @param sig.level significance level (Type I error probability)
#' @param mu1 arithmetic mean of group 1
#' @param mu2 arithmetic mean of group 2
#' @param gmu1 geometric mean of group 1
#' @param gmu2 geometric mean of group 2
#' @param trials number of trials in simulation
#' @param M number of simulations used in CAT method, see Chang (2011)
#' @param equal.sample equal sample sizes for two groups, see details
#' @param equal.shape assume the shape parameters are equal for two groups, see details
#' @param trace if positive, sample size and power are printed during the running of each simulation
#' @references Chang et al. (2011). Testing the equality of several gamma means: a parametric bootstrap method with applications. \emph{Computational Statistics}, \bold{26}:55-76.
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with method element.
#' @examples
#' # Calculate power, equal sizes
#' power_Gamma(n1 = 50, mu1 = 1, mu2 = 1.5, gmu1 = 0.6, gmu2 = 0.6, M = 100)
#' @importFrom stats uniroot rgamma quantile
#' @importFrom rootSolve multiroot
#' @export
power_Gamma <- function(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05, 
                        mu1 = NULL, mu2 = NULL, gmu1 = NULL, gmu2 = NULL,
                        trials = 100, M = 10000, equal.sample = TRUE, equal.shape = NULL,
                        trace = FALSE){
  if(!is.null(n1)&!is.null(n2)){
    if(n1==n2)
      equal.sample <- TRUE
    else
      equal.sample <- FALSE
  }
  if(equal.sample){
    n2 <- n1
  }
  # Test if shape parametes are equal
  equal_shape_test <- function(n_1, n_2){
    ## Step 1
    R1 <- mu1/gmu1
    R2 <- mu2/gmu2
    delta1hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(R1), c(1e-7,1e7))$root
    delta2hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(R2), c(1e-7,1e7))$root
    gamma1hat <- log(delta1hat)
    gamma2hat <- log(delta2hat)
    gammabarhat <- mean(c(gamma1hat,gamma2hat))
    etahat <- (gamma1hat-gammabarhat)^2 + (gamma2hat-gammabarhat)^2

    ## Step 2
    ### i
    w1 <- n_1/(n_1+n_2)
    w2 <- n_2/(n_1+n_2)
    Step2deltahat <- uniroot(function(delta) w1*log(delta/mu1)+w2*log(delta/mu2)-digamma(delta)+w1*log(gmu1)+w2*log(gmu2), c(1e-7,1e7))$root
    Step2lambda1hat <- Step2deltahat/mu1
    Step2lambda2hat <- Step2deltahat/mu2
    ### ii
    Step2etahat <- sapply(1:M, function(s){
      Step2X1 <- rgamma(n_1,shape = Step2deltahat,rate = Step2lambda1hat)
      Step2X2 <- rgamma(n_2,shape = Step2deltahat,rate = Step2lambda2hat)

      Step2X1AM <- mean(Step2X1)
      Step2X2AM <- mean(Step2X2)
      Step2X1GM <- exp(mean(log(Step2X1)))
      Step2X2GM <- exp(mean(log(Step2X2)))
      Step2R1 <- Step2X1AM/Step2X1GM
      Step2R2 <- Step2X2AM/Step2X2GM

      Step2iidelta1hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(Step2R1), c(1e-7,1e7))$root
      Step2iidelta2hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(Step2R2), c(1e-7,1e7))$root
      Step2gamma1hat <- log(Step2iidelta1hat)
      Step2gamma2hat <- log(Step2iidelta2hat)
      Step2gammabarhat <- mean(c(Step2gamma1hat,Step2gamma2hat))
      Step2etahat <- (Step2gamma1hat-Step2gammabarhat)^2 + (Step2gamma2hat-Step2gammabarhat)^2
      return(Step2etahat)
    })

    ## Step 3
    eta_pvalue <- mean(Step2etahat>etahat)
    return(eta_pvalue>=sig.level)
  }

  # Case 1: equal shape parameter delta1 = delta2
  power_equal_shape <- function(n_1, n_2){
    RejectTest <- sapply(1:trials, function(i){
      ## Step 1
      w1 <- n_1/(n_1+n_2)
      w2 <- n_2/(n_1+n_2)
      R1 <- mu1/gmu1
      R2 <- mu2/gmu2
      Step1C <- w1*log(R1)+w2*log(R2)
      deltahat <- uniroot(function(delta) log(delta)-digamma(delta)-Step1C, c(1e-7,1e7))$root
      lambda1hat <- deltahat/mu1
      lambda2hat <- deltahat/mu2
      beta1hat <- log(lambda1hat)
      beta2hat <- log(lambda2hat)
      betabarhat <- mean(c(beta1hat,beta2hat))
      etahat <- (beta1hat-betabarhat)^2 + (beta2hat-betabarhat)^2
      
      ## Step 2
      ### i
      XAM <- w1*mu1+w2*mu2
      XGM <- gmu1^w1*gmu2^w2
      R0 <- XAM/XGM
      Step2C <- log(R0)
      deltahat <- uniroot(function(delta) log(delta)-digamma(delta)-Step2C, c(1e-7,1e7))$root
      lambdahat <- deltahat/XAM
      ### ii
      Step2etahat <- sapply(1:M, function(s){
        Step2X <- rgamma(n_1+n_2,shape = deltahat,rate = lambdahat)
        Step2X1 <- Step2X[1:n_1]
        Step2X2 <- Step2X[-c(1:n_1)]
        
        Step2X1AM <- mean(Step2X1)
        Step2X2AM <- mean(Step2X2)
        Step2X1GM <- exp(mean(log(Step2X1)))
        Step2X2GM <- exp(mean(log(Step2X2)))
        Step2R1 <- Step2X1AM/Step2X1GM
        Step2R2 <- Step2X2AM/Step2X2GM
        Step2Cii <- w1*log(Step2R1)+w2*log(Step2R2)
        Step2deltahat <- uniroot(function(delta) log(delta)-digamma(delta)-Step2Cii, c(1e-7,1e7))$root
        Step2lambda1hat <- Step2deltahat/Step2X1AM
        Step2lambda2hat <- Step2deltahat/Step2X2AM
        Step2beta1hat <- log(Step2lambda1hat)
        Step2beta2hat <- log(Step2lambda2hat)
        Step2betabarhat <- mean(c(Step2beta1hat,Step2beta2hat))
        Step2etahat <- (Step2beta1hat-Step2betabarhat)^2 + (Step2beta2hat-Step2betabarhat)^2
        return(Step2etahat)
      })
      
      ## Step 3
      etaU <- quantile(Step2etahat,(1-sig.level))
      return(as.numeric(etahat>etaU))
    })
    
    return(mean(RejectTest))
  }
  
  # Case 2: unequal shape parameter delta1 != delta2
  power_unequal_shape <- function(n_1, n_2){
    RejectTest <- sapply(1:trials, function(i){
      ## Step 1
      R1 <- mu1/gmu1
      R2 <- mu2/gmu2
      delta1hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(R1), c(1e-7,1e7))$root
      delta2hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(R2), c(1e-7,1e7))$root
      lambda1hat <- delta1hat/mu1
      lambda2hat <- delta2hat/mu2
      mu1hat <- delta1hat/lambda1hat
      mu2hat <- delta2hat/lambda2hat
      beta1hat <- log(mu1hat)
      beta2hat <- log(mu2hat)
      betabarhat <- mean(c(beta1hat,beta2hat))
      etahat <- (beta1hat-betabarhat)^2 + (beta2hat-betabarhat)^2
      
      ## Step 2
      ### i
      lambda_model <- function(lambda) c(F1 = 1+log(lambda[1]) - digamma(lambda[3]*lambda[1]) + log(gmu1) - mu1/lambda[3],
                                         F2 = 1+log(lambda[2]) - digamma(lambda[3]*lambda[2]) + log(gmu2) - mu2/lambda[3],
                                         F3 = lambda[1] * n_1 * (mu1/lambda[3] - 1) + lambda[2] * n_2 * (mu2/lambda[3] - 1))
      Step2lambdas <- rootSolve::multiroot(f = lambda_model, start = c(lambda1hat, lambda2hat, mu1hat))$root
      Step2delta1hat <- Step2lambdas[1]*Step2lambdas[3]
      Step2delta2hat <- Step2lambdas[2]*Step2lambdas[3]
      ### ii
      Step2etahat <- sapply(1:M, function(s){
        Step2X1 <- rgamma(n_1,shape = Step2delta1hat,rate = Step2lambdas[1])
        Step2X2 <- rgamma(n_2,shape = Step2delta2hat,rate = Step2lambdas[2])
        
        Step2X1AM <- mean(Step2X1)
        Step2X2AM <- mean(Step2X2)
        Step2X1GM <- exp(mean(log(Step2X1)))
        Step2X2GM <- exp(mean(log(Step2X2)))
        Step2R1 <- Step2X1AM/Step2X1GM
        Step2R2 <- Step2X2AM/Step2X2GM
        
        Step2iidelta1hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(Step2R1), c(1e-7,1e7))$root
        Step2iidelta2hat <- uniroot(function(delta) log(delta)-digamma(delta)-log(Step2R2), c(1e-7,1e7))$root
        Step2iilambda1hat <- Step2iidelta1hat/Step2X1AM
        Step2iilambda2hat <- Step2iidelta2hat/Step2X2AM
        Step2mu1hat <- Step2iidelta1hat/Step2iilambda1hat
        Step2mu2hat <- Step2iidelta2hat/Step2iilambda2hat
        Step2beta1hat <- log(Step2mu1hat)
        Step2beta2hat <- log(Step2mu2hat)
        Step2betabarhat <- mean(c(Step2beta1hat,Step2beta2hat))
        Step2etahat <- (Step2beta1hat-Step2betabarhat)^2 + (Step2beta2hat-Step2betabarhat)^2
        return(Step2etahat)
      })
      
      ## Step 3
      etaU <- quantile(Step2etahat,(1-sig.level))
      return(as.numeric(etahat>etaU))
    })
    
    return(mean(RejectTest))
  }
  
  # calculate power with known sample size
  if(is.null(power)){
    if(is.null(equal.shape)){
      equal.shape <- equal_shape_test(n_1 = n1, n_2 = n2)
    }
    if(equal.sample){
      if(equal.shape){
        METHOD = "Two-sample Gamma Means Tests (Equal Sizes, Equal Shape Parameters)"
      }
      else{
        METHOD = "Two-sample Gamma Means Tests (Equal Sizes, Unequal Shape Parameters)"
      }
    }
    else{
      if(equal.shape){
        METHOD = "Two-sample Gamma Means Tests (Different Sizes, Equal Shape Parameters)"
      }
      else{
        METHOD = "Two-sample Gamma Means Tests (Different Sizes, Unequal Shape Parameters)"
      }
    }
    if(equal.shape){
      power <- power_equal_shape(n_1 = n1, n_2 = n2)
      return(structure(list(n1 = n1, n2 = n2, mu1 = mu1, mu2 = mu2,
                            sig.level = sig.level, power = power, method = METHOD),
                       class = "power.htest"))
    }
    else{
      power <- power_unequal_shape(n_1 = n1, n_2 = n2)
      return(structure(list(n1 = n1, n2 = n2, mu1 = mu1, mu2 = mu2,
                            sig.level = sig.level, power = power, method = METHOD),
                       class = "power.htest"))
    }
  }
  
  # calculate sample size with known power
  if(is.null(n1)|is.null(n2)){
    if(equal.sample){
      ## step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
      reach.flag <- 0
      ss.lower <- ss.upper <- sample.size.starting <- 100
      if(is.null(equal.shape)){
        equal.shape <- equal_shape_test(n_1 = ss.lower, n_2 = ss.upper)
      }
      if(equal.shape){
        METHOD = "Two-sample Gamma Means Tests (Equal Sizes, Equal Shape Parameters)"
        power.lower <- power.upper <- do.call("power_equal_shape",list(n_1 = sample.size.starting, n_2 = sample.size.starting))
      }
      else{
        METHOD = "Two-sample Gamma Means Tests (Equal Sizes, Unequal Shape Parameters)"
        power.lower <- power.upper <- do.call("power_unequal_shape",list(n_1 = sample.size.starting, n_2 = sample.size.starting))
      }
      while(reach.flag == 0){
        if(power.lower > power){
          ss.upper <- ss.lower
          power.upper <- power.lower
          # sample size should be greater than or equal to 3
          ss.lower <- max(floor(ss.lower/2),3)
          if(equal.shape){
            power.lower <-  do.call("power_equal_shape",list(n_1 = ss.lower, n_2 = ss.lower))
          }
          else{
            power.lower <-  do.call("power_unequal_shape",list(n_1 = ss.lower, n_2 = ss.lower))
          }
        }
        if(power.upper < power){
          ss.lower <- ss.upper
          power.lower <- power.upper
          ss.upper <- ss.upper*2
          if(equal.shape){
            power.upper <-  do.call("power_equal_shape",list(n_1 = ss.upper, n_2 = ss.upper))
          }
          else{
            power.upper <-  do.call("power_unequal_shape",list(n_1 = ss.upper, n_2 = ss.upper))
          }
        }
        if((power.lower <= power & power <= power.upper)|(ss.upper <= 3)){
          reach.flag <- 1
        }
        if(ss.upper>1e5){
          stop("The required sample size is greater than 100,000.")
        }
        if(ss.upper>1e3){
          warning("The target power is not reached with a large sample size (> 1000).\n A test with normal approximation might be considered")
        }
        if(trace){
          cat(paste0("Current range of sample size (power) is [",
                     ss.lower, ", ", ss.upper, "] ([", 
                     power.lower, ", ", power.upper, "]).\n"))
        }
      }
      ## step 2: find exact sample size
      reach.flag <- 0
      while(reach.flag == 0){
        if((ss.upper-ss.lower)<=1){
          reach.flag <- 1
          sample.size.min <- ss.upper
          power.output <- power.upper
        }
        else{
          sample.size.midpt <- (ss.upper+ss.lower)%/%2
          if(equal.shape){
            power.midpt <- do.call("power_equal_shape",list(n_1 = sample.size.midpt, n_2 = sample.size.midpt))
          }
          else{
            power.midpt <- do.call("power_unequal_shape",list(n_1 = sample.size.midpt, n_2 = sample.size.midpt))
          }
          
          if(power.midpt < power){
            ss.lower <- sample.size.midpt
            power.lower <- power.midpt
          }
          else{
            ss.upper <- sample.size.midpt
            power.upper <- power.midpt
          }
        }
        if(trace){
          cat(paste0("Current range of sample size (power) is [",
                     ss.lower, ", ", ss.upper, "] ([", 
                     power.lower, ", ", power.upper, "]).\n"))
        }
      }
      return(structure(list(n1 = ss.upper, n2 = ss.upper, mu1 = mu1, mu2 = mu2,
                            sig.level = sig.level, power = power, method = METHOD),
                       class = "power.htest"))
    }
    
    else{
      if(is.null(n1)){
        ## step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
        reach.flag <- 0
        ss.lower <- ss.upper <- sample.size.starting <- 100
        if(is.null(equal.shape)){
          equal.shape <- equal_shape_test(n_1 = ss.lower, n_2 = n2)
        }
        if(equal.shape){
          METHOD = "Two-sample Gamma Means Tests (Different Sizes, Equal Shape Parameters)"
          power.lower <- power.upper <- do.call("power_equal_shape",list(n_1 = sample.size.starting, n_2 = n2))
        }
        else{
          METHOD = "Two-sample Gamma Means Tests (Different Sizes, Unequal Shape Parameters)"
          power.lower <- power.upper <- do.call("power_unequal_shape",list(n_1 = sample.size.starting, n_2 = n2))
        }
        while(reach.flag == 0){
          if(power.lower > power){
            ss.upper <- ss.lower
            power.upper <- power.lower
            # sample size should be greater than or equal to 3
            ss.lower <- max(floor(ss.lower/2),3)
            if(equal.shape){
              power.lower <-  do.call("power_equal_shape",list(n_1 = ss.lower, n_2 = n2))
            }
            else{
              power.lower <-  do.call("power_unequal_shape",list(n_1 = ss.lower, n_2 = n2))
            }
          }
          if(power.upper < power){
            ss.lower <- ss.upper
            power.lower <- power.upper
            ss.upper <- ss.upper*2
            if(equal.shape){
              power.upper <-  do.call("power_equal_shape",list(n_1 = ss.upper, n_2 = n2))
            }
            else{
              power.upper <-  do.call("power_unequal_shape",list(n_1 = ss.upper, n_2 = n2))
            }
          }
          if((power.lower <= power & power <= power.upper)|(ss.upper <= 3)){
            reach.flag <- 1
          }
          if(ss.upper>1e5){
            stop("The required sample size is greater than 100,000")
          }
          if(ss.upper>1e3){
            warning("The target power is not reached with a large sample size (> 1000).\n A test with normal approximation might be considered")
          }
          if(trace){
            cat(paste0("Current range of sample size (power) is [",
                       ss.lower, ", ", ss.upper, "] ([", 
                       power.lower, ", ", power.upper, "]).\n"))
          }
        }
        ## step 2: find exact sample size
        reach.flag <- 0
        while(reach.flag == 0){
          if((ss.upper-ss.lower)<=1){
            reach.flag <- 1
            sample.size.min <- ss.upper
            power.output <- power.upper
          }
          else{
            sample.size.midpt <- (ss.upper+ss.lower)%/%2
            if(equal.shape){
              power.midpt <- do.call("power_equal_shape",list(n_1 = sample.size.midpt, n_2 = n2))
            }
            else{
              power.midpt <- do.call("power_unequal_shape",list(n_1 = sample.size.midpt, n_2 = n2))
            }
            
            if(power.midpt < power){
              ss.lower <- sample.size.midpt
              power.lower <- power.midpt
            }
            else{
              ss.upper <- sample.size.midpt
              power.upper <- power.midpt
            }
          }
          if(trace){
            cat(paste0("Current range of sample size (power) is [",
                       ss.lower, ", ", ss.upper, "] ([", 
                       power.lower, ", ", power.upper, "]).\n"))
          }
        }
        return(structure(list(n1 = ss.upper, n2 = n2, mu1 = mu1, mu2 = mu2,
                              sig.level = sig.level, power = power, method = METHOD),
                         class = "power.htest"))
      }
      
      if(is.null(n2)){
        ## step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
        reach.flag <- 0
        ss.lower <- ss.upper <- sample.size.starting <- 100
        if(is.null(equal.shape)){
          equal.shape <- equal_shape_test(n_1 = n1, n_2 = ss.upper)
        }
        if(equal.shape){
          METHOD = "Two-sample Gamma Means Tests (Different Sizes, Equal Shape Parameters)"
          power.lower <- power.upper <- do.call("power_equal_shape",list(n_1 = n1, n_2 = sample.size.starting))
        }
        else{
          METHOD = "Two-sample Gamma Means Tests (Different Sizes, Unequal Shape Parameters)"
          power.lower <- power.upper <- do.call("power_unequal_shape",list(n_1 = n1, n_2 = sample.size.starting))
        }
        while(reach.flag == 0){
          if(power.lower > power){
            ss.upper <- ss.lower
            power.upper <- power.lower
            # sample size should be greater than or equal to 3
            ss.lower <- max(floor(ss.lower/2),3)
            if(equal.shape){
              power.lower <-  do.call("power_equal_shape",list(n_1 = n1, n_2 = ss.lower))
            }
            else{
              power.lower <-  do.call("power_unequal_shape",list(n_1 = n1, n_2 = ss.lower))
            }
          }
          if(power.upper < power){
            ss.lower <- ss.upper
            power.lower <- power.upper
            ss.upper <- ss.upper*2
            if(equal.shape){
              power.upper <-  do.call("power_equal_shape",list(n_1 = n1, n_2 = ss.upper))
            }
            else{
              power.upper <-  do.call("power_unequal_shape",list(n_1 = n1, n_2 = ss.upper))
            }
          }
          if((power.lower <= power & power <= power.upper)|(ss.upper <= 3)){
            reach.flag <- 1
          }
          if(ss.upper>1e5){
            stop("The required sample size is greater than 100,000")
          }
          if(ss.upper>1e3){
            warning("The target power is not reached with a large sample size (> 1000).\n A test with normal approximation might be considered")
          }
          if(trace){
            cat(paste0("Current range of sample size (power) is [",
                       ss.lower, ", ", ss.upper, "] ([", 
                       power.lower, ", ", power.upper, "]).\n"))
          }
        }
        ## step 2: find exact sample size
        reach.flag <- 0
        while(reach.flag == 0){
          if((ss.upper-ss.lower)<=1){
            reach.flag <- 1
            sample.size.min <- ss.upper
            power.output <- power.upper
          }
          else{
            sample.size.midpt <- (ss.upper+ss.lower)%/%2
            if(equal.shape){
              power.midpt <- do.call("power_equal_shape",list(n_1 = n1, n_2 = sample.size.midpt))
            }
            else{
              power.midpt <- do.call("power_unequal_shape",list(n_1 = n1, n_2 = sample.size.midpt))
            }
            
            if(power.midpt < power){
              ss.lower <- sample.size.midpt
              power.lower <- power.midpt
            }
            else{
              ss.upper <- sample.size.midpt
              power.upper <- power.midpt
            }
          }
          if(trace){
            cat(paste0("Current range of sample size (power) is [",
                       ss.lower, ", ", ss.upper, "] ([", 
                       power.lower, ", ", power.upper, "]).\n"))
          }
        }
        return(structure(list(n1 = n1, n2 = ss.upper, mu1 = mu1, mu2 = mu2,
                              sig.level = sig.level, power = power, method = METHOD),
                         class = "power.htest"))
        
      }
    }
    
  }

}
