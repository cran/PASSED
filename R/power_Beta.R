#' @title Power Calculations for Test of Two Beta Means
#' @description  Compute the power for a test of two sample means with beta distributions, or determine the minimum sample size to obtain a target power.
#' @details Exactly one of the parameters \code{n1}, \code{n2} and \code{power} must be passed as NULL, and that parameter is determined from the others.\cr\cr
#' This function allows you to set the number of trials in the simulation to control the result accuracy, 
#' and type of link used in the beta regression. You can choose one of the following: "logit", "probit", "cloglog", "cauchit", "log", "loglog".
#' @usage power_Beta(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05, 
#' mu1 = NULL, sd1 = NULL, mu2 = NULL, equal.sample = TRUE,
#' trials = 100, equal.precision = TRUE, sd2 = NULL, 
#' link.type = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"))
#' @param n1 sample size in group 1, or sample size in each group if \code{equal.sample = TRUE}
#' @param n2 sample size in group 2 
#' @param power power of test (1 minus Type II error probability)
#' @param sig.level significance level (Type I error probability)
#' @param mu1 sample mean of group 1
#' @param sd1 standard deviation for group 1
#' @param mu2 sample mean of group 2
#' @param equal.sample equal sample sizes for two groups, see details
#' @param trials number of trials in simulation
#' @param equal.precision equal dispersion parameter assumption in simulation
#' @param sd2 standard deviation for group 2. Only applicable when \code{equal.precision = FALSE}
#' @param link.type type of link used in the beta regression, see details
#' @return Object of class "power.htest", a list of the arguments (including the computed one) augmented with method and note elements.
#' @examples 
#' # calculate power
#' power_Beta(mu1 = 0.5, mu2 = 0.80, sd1 = 0.25, n1 = 60)
#' # calculate sample size for both groups
#' power_Beta(mu1 = 0.5, mu2 = 0.80, sd1 = 0.25, power=0.8)
#' @importFrom stats rbeta wilcox.test pnorm
#' @export

power_Beta <- function(n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05, 
                       mu1 = NULL, sd1 = NULL, mu2 = NULL,  equal.sample = TRUE,
                       trials = 100, equal.precision = TRUE, sd2 = NULL, 
                       link.type = c("logit", "probit", "cloglog", "cauchit", "log", "loglog")){
  if(is.null(mu1))
    stop("'mu1' must be given")
  if(is.null(mu2))
    stop("'mu2' must be given")
  if(is.null(sd1))
    stop("'sd1' must be given")
  if(!is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
    stop("'sig.level' must be numeric in (0, 1)")
  link.type <- match.arg(link.type)
  # Set parameters
  phi<- ((mu1*(1-mu1))/(sd1*sd1))-1
  if(phi < 0){
    stop("phi must be greater than 0")
  }
  a0 <- mu1*phi
  b0 <- (1-mu1)*phi
  if(equal.precision == TRUE){
    a1 <- mu2*phi
    b1 <- (1-mu2)*phi
  }
  else{
    if(is.null(sd2)==TRUE){
      stop("'sd2' must be given with equal precision assumption")
    }
    else{
      phi1 <- ((mu2*(1-mu2))/(sd2*sd2))-1
      if(phi1 < 0){
        stop("phi1 must be greater than 0")
      }
      a1 <- mu2*phi1
      b1 <- (1-mu2)*phi1
    }
  }
  
  if(!is.null(n1)&!is.null(n2)){
    if(n1 != n2){
      equal.sample <- FALSE
    }
    else{
      equal.sample <- TRUE
    }
  }
  if(equal.sample){
    if (sum(sapply(list(n1, power), is.null)) != 1) 
      stop("exactly one of n1 and power must be NULL")
    NOTE <- "N is number in *each* group"
    # define power calculation function
    betapwr_one_samp <- function(sampsize){
      betapwr.base <- function(){
        Y.H0 <- cbind(rep(1:sampsize,trials),rep(1:trials,rep(sampsize,trials)),rbeta(sampsize*trials,a0,b0))
        Y.Ha <- cbind(rep(1:sampsize,trials),rep(1:trials,rep(sampsize,trials)),rbeta(sampsize*trials,a1,b1))
        
        ## Combine Y.H0 and Y.H1
        Y.mat <- rbind(Y.H0,Y.Ha)
        colnames(Y.mat) <- c( "sample","trials","y")
        tmt <-c(rep(0,(trials*sampsize)),rep(1,(trials*sampsize)))
        ## Combine "sample trial y" with "tmt"(0,1)
        ## Set simulation matrix as sim, ordered by trials
        sim <- data.frame(Y.mat,tmt)  
        
        if(max(sim[,3]) > (1-1e-16) | min(sim[,3]) < 1e-16){
          sim[,3] <- (sim[,3] * (sampsize * 2 - 1) + 0.5) / (sampsize * 2)
        }
        
        outtest <- sapply(1:trials, function(i){
          sub.sim <-  subset(sim, trials == i)
          X <- cbind(rep(1,nrow(sub.sim)),sub.sim$tmt)
          colnames(X) <- c("(Intercept)","tmt")
          
          fit1 <- suppressWarnings(do.call(betareg::betareg.fit,list(x=X, y=as.numeric(sub.sim$y), link = link.type,type ="ML")))
          cf <- as.vector(do.call("c",fit1$coefficients))
          se <- sqrt(diag(fit1$vcov))
          wald.pvalue <- 2*pnorm(-abs(cf/se))[2]
          
          return(wald.pvalue)
        })
        Power = mean(as.numeric(outtest<sig.level))
        return(Power)
      }
      Power <- tryCatch(betapwr.base(),error=function(e){return(NA)})
      while(is.na(Power[1])){
        warning("Simulation failed with current seed because of 0 or 1 aggregations, re-running with new seed")
        Power <- tryCatch(betapwr.base(),error=function(e){return(NA)})
      }
      return(Power)
    }
    
    # Calculate power
    if(is.null(power)){
      Power <- betapwr_one_samp(n1)
      if(equal.precision==TRUE){
        METHOD = paste0("Two-sample Beta Means Tests (Equal Sizes) (",link.type," link, equal precision)")
        return(structure(list(N = n1, mu1 = mu1, mu2 = mu2, sd1 = sd1, sig.level = sig.level, power = Power, 
                              method = METHOD, note = NOTE), 
                         class = "power.htest"))
      }
      else{
        METHOD = paste0("Two-sample Beta Means Tests (Equal Sizes) (",link.type," link, unequal precision)")
        return(structure(list(N = n1, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, sig.level = sig.level, power = Power, 
                              method = METHOD, note = NOTE), 
                         class = "power.htest"))
      }
    }
    
    #Calculate sample size
    else{
      # use two-sample t test to get a starting value estimation
      sample.size.starting <- round(do.call("power.t.test",list(delta = (mu2-mu1), sd = sd1, sig.level = sig.level,power = power))$n,0)
      
      # step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
      reach.flag <- 0
      ss.lower <- ss.upper <- sample.size.starting
      power.lower <- power.upper <- do.call("betapwr_one_samp",list(sampsize = sample.size.starting))
      while(reach.flag == 0){
        if(power.lower > power){
          ss.upper <- ss.lower
          power.upper <- power.lower
          # sample size should be greater than or equal to 3
          ss.lower <- max(floor(ss.lower/2),3)
          power.lower <-  do.call("betapwr_one_samp",list(sampsize = ss.lower))
        }
        if(power.upper < power){
          ss.lower <- ss.upper
          power.lower <- power.upper
          ss.upper <- ss.upper*2
          power.upper <-  do.call("betapwr_one_samp",list(sampsize = ss.upper))
        }
        if((power.lower <= power & power <= power.upper)|(ss.upper <= 3)){
          reach.flag <- 1
        }
      }
      
      # step 2: find exact sample size
      reach.flag <- 0
      while(reach.flag == 0){
        if((ss.upper-ss.lower)<=1){
          reach.flag <- 1
          sample.size.min <- ss.upper
          power.output <- power.upper
        }
        else{
          sample.size.midpt <- (ss.upper+ss.lower)%/%2
          power.midpt <- do.call("betapwr_one_samp",list(sampsize = sample.size.midpt))
          if(power.midpt < power){
            ss.lower <- sample.size.midpt
            power.lower <- power.midpt
          }
          else{
            ss.upper <- sample.size.midpt
            power.upper <- power.midpt
          }
        }
      }
      
      
      if(equal.precision==TRUE){
        METHOD = paste0("Two-sample Beta Means Tests (Equal Sizes) (",link.type," link, equal precision)")
        return(structure(list(N = sample.size.min, mu1 = mu1, mu2 = mu2, sd1 = sd1, sig.level = sig.level, power = power.output, 
                              method = METHOD, note = NOTE), 
                         class = "power.htest"))
      }
      else{
        METHOD = paste0("Two-sample Beta Means Tests (Equal Sizes) (",link.type," link, unequal precision)")
        return(structure(list(N = sample.size.min, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, sig.level = sig.level, power = power.output, 
                              method = METHOD, note = NOTE), 
                         class = "power.htest"))
      }
    }
  }
  else{
    if (sum(sapply(list(n1, n2, power), is.null)) != 1) 
      stop("exactly one of n1, n2 and power must be NULL")
    # define power calculation function
    betapwr <- function(sampsize1, sampsize2){
      betapwr.base <- function(){
        Y.H0 <- cbind(rep(1:sampsize1,trials),rep(1:trials,rep(sampsize1,trials)),rbeta(sampsize1*trials,a0,b0))
        Y.Ha <- cbind(rep(1:sampsize2,trials),rep(1:trials,rep(sampsize2,trials)),rbeta(sampsize2*trials,a1,b1))
        
        ## Combine Y.H0 and Y.H1
        Y.mat <- rbind(Y.H0,Y.Ha)
        colnames(Y.mat) <- c( "sample","trials","y")
        tmt <-c(rep(0,(trials*sampsize1)),rep(1,(trials*sampsize2)))
        ## Combine "sample trial y" with "tmt"(0,1)
        ## Set simulation matrix as sim, ordered by trials
        sim <- data.frame(Y.mat,tmt)  
        
        if(max(sim[,3]) > (1-1e-16) | min(sim[,3]) < 1e-16){
          sim[,3] <- (sim[,3] * (sampsize1 + sampsize2 - 1) + 0.5) / (sampsize1 + sampsize2)
        }
        
        outtest <- sapply(1:trials, function(i){
          sub.sim <-  subset(sim, trials == i)
          X <- cbind(rep(1,nrow(sub.sim)),sub.sim$tmt)
          colnames(X) <- c("(Intercept)","tmt")
          
          fit1 <- suppressWarnings(do.call(betareg::betareg.fit,list(x=X, y=as.numeric(sub.sim$y), link = link.type,type ="ML")))
          cf <- as.vector(do.call("c",fit1$coefficients))
          se <- sqrt(diag(fit1$vcov))
          wald.pvalue <- 2*pnorm(-abs(cf/se))[2]
          
          return(wald.pvalue)
        })
        Power = mean(as.numeric(outtest<sig.level))
        return(Power)
      }
      Power <- tryCatch(betapwr.base(),error=function(e){return(NA)})
      while(is.na(Power[1])){
        warning("Simulation failed with current seed because of 0 or 1 aggregations, re-running with new seed")
        Power <- tryCatch(betapwr.base(),error=function(e){return(NA)})
      }
      return(Power)
    }
    
    # Calculate power
    if(is.null(power)){
      Power <- betapwr(n1, n2)
      if(equal.precision==TRUE){
        METHOD = paste0("Two-sample Beta Means Tests (Different Sizes) (",link.type," link, equal precision)")
        return(structure(list(n1 = n1, n2 = n2, mu1 = mu1, mu2 = mu2, sd1 = sd1, sig.level = sig.level, power = Power, 
                              method = METHOD), 
                         class = "power.htest"))
      }
      else{
        METHOD = paste0("Two-sample Beta Means Tests (Different Sizes) (",link.type," link, unequal precision)")
        return(structure(list(n1 = n1, n2 = n2, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, sig.level = sig.level, power = Power, 
                              method = METHOD), 
                         class = "power.htest"))
      }
    }
    
    #Calculate sample size
    else{
      if(is.null(n1)){
        # use two-sample t test to get a starting value estimation
        sample.size.starting <- round(do.call("power.t.test",list(delta = (mu2-mu1), sd = sd1, sig.level = sig.level,power = power))$n,0)
        
        # step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
        reach.flag <- 0
        ss.lower <- ss.upper <- sample.size.starting
        power.lower <- power.upper <- do.call("betapwr",list(sampsize1 = sample.size.starting, sampsize2 = n2))
        while(reach.flag == 0){
          if(power.lower > power){
            ss.upper <- ss.lower
            power.upper <- power.lower
            # sample size should be greater than or equal to 3
            ss.lower <- max(floor(ss.lower/2),3)
            power.lower <-  do.call("betapwr",list(sampsize1 = ss.lower, sampsize2 = n2))
          }
          if(power.upper < power){
            ss.lower <- ss.upper
            power.lower <- power.upper
            ss.upper <- ss.upper*2
            power.upper <-  do.call("betapwr",list(sampsize1 = ss.upper, sampsize2 = n2))
          }
          if((power.lower <= power & power <= power.upper)|(ss.upper <= 3)){
            reach.flag <- 1
          }
        }
        
        # step 2: find exact sample size
        reach.flag <- 0
        while(reach.flag == 0){
          if((ss.upper-ss.lower)<=1){
            reach.flag <- 1
            sample.size.min <- ss.upper
            power.output <- power.upper
          }
          else{
            sample.size.midpt <- (ss.upper+ss.lower)%/%2
            power.midpt <- do.call("betapwr",list(sampsize1 = sample.size.midpt, sampsize2 = n2))
            if(power.midpt < power){
              ss.lower <- sample.size.midpt
              power.lower <- power.midpt
            }
            else{
              ss.upper <- sample.size.midpt
              power.upper <- power.midpt
            }
          }
        }
        if(equal.precision==TRUE){
          METHOD = paste0("Two-sample Beta Means Tests (Different Sizes) (",link.type," link, equal precision)")
          return(structure(list(n1 = sample.size.min, n2 = n2, mu1 = mu1, mu2 = mu2, sd1 = sd1, sig.level = sig.level, power = power.output, 
                                method = METHOD), 
                           class = "power.htest"))
        }
        else{
          METHOD = paste0("Two-sample Beta Means Tests (Different Sizes) (",link.type," link, unequal precision)")
          return(structure(list(n1 = sample.size.min, n2 = n2, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, sig.level = sig.level, power = power.output, 
                                method = METHOD), 
                           class = "power.htest"))
        }
      }
      if(is.null(n2)){
        # use two-sample t test to get a starting value estimation
        sample.size.starting <- round(do.call("power.t.test",list(delta = (mu2-mu1), sd = sd1, sig.level = sig.level,power = power))$n,0)
        
        # step 1: get an interval of sample size [ss.lower, ss.upper], which satisfies power.ss.lower < target power < power.ss.upper
        reach.flag <- 0
        ss.lower <- ss.upper <- sample.size.starting
        power.lower <- power.upper <- do.call("betapwr",list(sampsize1 = n1, sampsize2 = sample.size.starting))
        while(reach.flag == 0){
          if(power.lower > power){
            ss.upper <- ss.lower
            power.upper <- power.lower
            # sample size should be greater than or equal to 3
            ss.lower <- max(floor(ss.lower/2),3)
            power.lower <-  do.call("betapwr",list(sampsize1 = n1, sampsize2 = ss.lower))
          }
          if(power.upper < power){
            ss.lower <- ss.upper
            power.lower <- power.upper
            ss.upper <- ss.upper*2
            power.upper <-  do.call("betapwr",list(sampsize1 = n1, sampsize2 = ss.upper))
          }
          if((power.lower <= power & power <= power.upper)|(ss.upper <= 3)){
            reach.flag <- 1
          }
        }
        
        # step 2: find exact sample size
        reach.flag <- 0
        while(reach.flag == 0){
          if((ss.upper-ss.lower)<=1){
            reach.flag <- 1
            sample.size.min <- ss.upper
            power.output <- power.upper
          }
          else{
            sample.size.midpt <- (ss.upper+ss.lower)%/%2
            power.midpt <- do.call("betapwr",list(sampsize1 = n1, sampsize2 = sample.size.midpt))
            if(power.midpt < power){
              ss.lower <- sample.size.midpt
              power.lower <- power.midpt
            }
            else{
              ss.upper <- sample.size.midpt
              power.upper <- power.midpt
            }
          }
        }
        if(equal.precision==TRUE){
          METHOD = paste0("Two-sample Beta Means Tests (Different Sizes) (",link.type," link, equal precision)")
          return(structure(list(n1 = n1, n2 = sample.size.min, mu1 = mu1, mu2 = mu2, sd1 = sd1, sig.level = sig.level, power = power.output, 
                                method = METHOD), 
                           class = "power.htest"))
        }
        else{
          METHOD = paste0("Two-sample Beta Means Tests (Different Sizes) (",link.type," link, unequal precision)")
          return(structure(list(n1 = n1, n2 = sample.size.min, mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, sig.level = sig.level, power = power.output, 
                                method = METHOD), 
                           class = "power.htest"))
        }
      }
    }
  }
}