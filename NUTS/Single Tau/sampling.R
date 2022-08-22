#############################################################
# R script for generating time series from artifical/inferred
# parameter values for single-tau model
#############################################################

#### Necessary libraries ####
library(reshape2)

result_fit <- function(n, time, tau, sigma_x, x_hat, t_switch) {
  # dt <- 1 time difference is one second
  # kbT <- 1 set Boltzmann-temperature as 1
  n <- n # determines the amount of intervals after t_sw (approx. 10)


  #### Dummy parameters ####
  sigma_x_d <- sigma_x
  t_d <- time
  tau_d <- vector(mode = "numeric", length = n)
  t_switch_d <- vector(mode = "numeric", length = n)
  x_hat_d <- vector(mode = "numeric", length = n)
  for (i in 1:n) {
    tau_d[i] <- tau[i]
    t_switch_d[i] <- t_switch[i]
    x_hat_d[i] <- x_hat[i]
  }
  t_switch_d[n] <- t_d - 1

  #### Print overview ####
  cat("Tau: ", tau_d, "\n")
  cat("t_sw: ", t_switch_d, "\n")
  cat("x_hat: ", x_hat_d, "\n")
  cat("Sigma: ", sigma_x[1], "\n")

  #### Generating fake data ####
  x_d <- vector(mode = "numeric", length = t_d)
  x_d[1] <- rnorm(1, x_hat_d[1], sigma_x_d)
  m <- 1
  for (i in 2:t_d) {
    if (i - 1 <= t_switch_d[m]) {
      x_d[i] <- rnorm(
        1,
        (x_d[i - 1] - x_hat_d[m]) * exp(-1 / tau_d[m]) + x_hat_d[m],
        sigma_x_d * sqrt(1 - exp(-2 / tau_d[m]))
      )
    } else {
      t1_d <- t_switch_d[m] - (i - 2)
      t2_d <- (i - 1) - t_switch_d[m]
      x_d[i] <- rnorm(
        1,
        ((((x_d[i - 1] - x_hat_d[m]) * exp(-t1_d / tau_d[m]) + x_hat_d[m]) -
          x_hat_d[m + 1]) * exp(-t2_d / tau_d[m + 1])) + x_hat_d[m + 1],
        sigma_x_d * sqrt(1 - exp(-4 / (tau_d[m] + tau_d[m + 1])))
      )
      m <- m + 1
    }
  }
  return(x_d)
}