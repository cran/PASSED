## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  knitr::opts_chunk$set(dev = 'pdf')
)

## -----------------------------------------------------------------------------
# Load the PASSED package
library(PASSED)

# Set a random seed for reproducibility
set.seed(1)

# Input data for control group
mu1 = 0.0174  # Mean
sd1 = 0.0211  # Standard Deviation

# Target alternative mean (or the mean for treatment group) (25% reduction)
mu2 = 0.0131

# Calculate sample size
result <- power_Beta(mu1 = mu1, sd1 = sd1, mu2 = mu2, power = 0.80, link.type = "logit", trials = 100, equal.precision = TRUE)

# Print the result
print(result)

## -----------------------------------------------------------------------------
# Set seed for the simulation below
set.seed(1)

Ex1 <- mapply(
  function(mu2, sample_size){
    Betapower <- power_Beta(mu1 = 0.0174, sd1 = 0.0211,
                            mu2 = mu2, n1 = sample_size,
                            link.type = "logit", trials = 100,
                            equal.precision = TRUE)
    Normalpower <- power_Normal(delta = (0.0174 - mu2), n1 = sample_size,
                                sd1 = 0.0211, sd2 = 0.0211)
    return(c(Betapower$power,
             round(Normalpower$power,3),
             sample_size,
             mu2,
             0.0174))
  },
  # Range of mu2 was set as [0.0120, 0.0140] by 0.0010
  rep(seq(0.0120, 0.0140, 0.0010), 5),
  # Range of sample size was set as [100, 200] by 25
  rep(seq(100, 200, 25), rep(3, 5))
)

# Reform the output
Ex1 <- as.data.frame(t(Ex1))
# Set column names
colnames(Ex1) <- c("Power (Beta)",
                   "Power (Normal)",
                   "Sample Size",
                   "mu2",
                   "mu1")
# Display the results
print(Ex1)


## -----------------------------------------------------------------------------
# Set seed for the simulation below
set.seed(1)
# Set the simulations
Ex2 <- mapply(
  function(mu2, sample_size){
    Betapower <- power_Beta(mu1 = 0.0174, sd1 = 0.0211, sd2 = 0.030,
                            mu2 = mu2, n1 = sample_size,
                            link.type = "logit", trials = 100,
                            equal.precision = FALSE)
    Normalpower <- power_Normal(delta = (0.0174 - mu2), n1 = sample_size,
                                sd1 = 0.0211, sd2 = 0.030)
    return(c(Betapower$power,
             round(Normalpower$power,3),
             sample_size,
             mu2,
             0.0174))
  },
  # Range of mu2 was set as [0.0120, 0.0140] by 0.0010
  rep(seq(0.0120, 0.0140, 0.0010), 5),
  # Range of sample size was set as [100, 200] by 25
  rep(seq(100, 200, 25), rep(3, 5))
)
# Reform the output
Ex2 <- as.data.frame(t(Ex2))
# Set column names
colnames(Ex2) <- c("Power (Beta)",
                   "Power (Normal)",
                   "Sample Size",
                   "mu2",
                   "mu1")
# Display the results
print(Ex2)


## ---- out.width="100%", fig.width=7, fig.height=6-----------------------------

## Parameter settings
mu1 <- 0.0174
mu2 <- 0.0131
sd1 <- 0.0211
sd2 <- 0.0300
phi<- ((mu1*(1-mu1))/(sd1*sd1))-1
a0 <- mu1*phi
b0 <- (1-mu1)*phi
a1 <- mu2*phi
b1 <- (1-mu2)*phi

## Draw lines
Ex1_eq_p_1 <- sapply(seq(0.001,0.100,0.001), function(i) return(dbeta(i,a0, b0)))
plot(seq(0.001,0.100,0.001), Ex1_eq_p_1, xlab="percentage of SNF residents with pressure ulcers", ylab="density", type ="l", xlim = c(-0.001, 0.10), main = "equal precision/standard deviation", col="darkgreen", lwd=3,lty = 1)
Ex1_eq_p_2 <- sapply(seq(0.001,0.100,0.001), function(i) return(dbeta(i,a1, b1)))
lines(seq(0.001,0.100,0.001), Ex1_eq_p_2, col="darkgreen", lwd=3,lty=2)
Ex1_nm_p_1 <- sapply(seq(-0.00,0.100,0.001), function(i) return(dnorm(i,mu1, sd1)))
lines(seq(-0.00,0.100,0.001), Ex1_nm_p_1, col="darkred", lwd=3)
Ex1_nm_p_2 <- sapply(seq(-0.00,0.100,0.001), function(i) return(dnorm(i,mu2, sd1)))
lines(seq(-0.00,0.100,0.001), Ex1_nm_p_2, col="darkred", lwd=3,lty=2)
# Shading area
polygon(c(seq(0.001,0.100,0.001),rev(seq(0.001,0.100,0.001))),c(Ex1_eq_p_2,rev(Ex1_eq_p_1)),col=rgb(0, 1, 0,0.5), border = NA)
polygon(c(seq(-0.00,0.100,0.001),rev(seq(-0.00,0.100,0.001))),c(Ex1_nm_p_1,rev(Ex1_nm_p_2)),col=rgb(1, 0, 0,0.5), border = NA)
# Add legend
legend(0.050,80, c("Beta (control)","Beta (alternative)","Normal (control)", "Normal (alternative)"),lty=c(1,2,1,3),col=c("darkgreen","darkgreen","darkred","darkred"),lwd=c(3,3,3,3))

## Parameter settings
phi<- ((mu1*(1-mu1))/(sd1*sd1))-1
phi1 <- ((mu2*(1-mu2))/(sd2*sd2))-1
a0 <- mu1*phi
b0 <- (1-mu1)*phi
a1 <- mu2*phi1
b1 <- (1-mu2)*phi1

## Draw lines
Ex1_ne_p_1 <- sapply(seq(0.001,0.100,0.001), function(i) return(dbeta(i,a0, b0)))
plot(seq(0.001,0.100,0.001), Ex1_ne_p_1, xlab="percentage of SNF residents with pressure ulcers", ylab="density", type ="l",  xlim = c(-0.001, 0.10), main = "unequal precision/standard deviation", col="darkgreen", lwd=3,lty = 1)
Ex1_ne_p_2 <- sapply(seq(0.001,0.100,0.001), function(i) return(dbeta(i,a1, b1)))
lines(seq(0.001,0.100,0.001), Ex1_ne_p_2, col="darkgreen", lwd=3,lty=2)
Ex1_nm_p_1 <- sapply(seq(-0.00,0.100,0.001), function(i) return(dnorm(i,mu1, sd1)))
lines(seq(-0.00,0.100,0.001), Ex1_nm_p_1, col="darkred", lwd=3)
Ex1_nm_p_2 <- sapply(seq(-0.00,0.100,0.001), function(i) return(dnorm(i,mu2, sd2)))
lines(seq(-0.00,0.100,0.001), Ex1_nm_p_2, col="darkred", lwd=3,lty=2)
# Shading area
polygon(c(seq(0.001,0.100,0.001),rev(seq(0.001,0.100,0.001))),c(Ex1_ne_p_2,rev(Ex1_ne_p_1)),col=rgb(0, 1, 0,0.5), border = NA)
polygon(c(seq(-0.00,0.100,0.001),rev(seq(-0.00,0.100,0.001))),c(Ex1_nm_p_1,rev(Ex1_nm_p_2)),col=rgb(1, 0, 0,0.5), border = NA)
# Add legend
legend(0.050,80, c("Beta (control)","Beta (alternative)","Normal (control)", "Normal (alternative)"),lty=c(1,2,1,3),col=c("darkgreen","darkgreen","darkred","darkred"),lwd=c(3,3,3,3))




