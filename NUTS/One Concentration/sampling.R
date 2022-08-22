#############################################################
# R script for generating time series from artifical/inferred
# parameter values for one-concentration model
#############################################################

#### Necessary libraries ####
library(reshape2)

result_fit <- function(n, time, tau, sigma_x2, x1, x2, t_switch) {
  # dt <- 1 time difference is one second
  # kbT <- 1 set Boltzmann-temperature as 1
  n <- n

  #### Dummy parameters ####
  tau_d <- rnorm(n, tau[1], tau[2])
  sigma_x_d <- rnorm(n, sigma_x2[1], sigma_x2[2])
  t_d <- time
  t_switch_d <- runif(
    n,
    t_switch[1] - 2 * t_switch[2],
    t_switch[1] + 2 * t_switch[2]
  )
  x1_d <- rnorm(n, x1[1], x1[2])
  x2_d <- rnorm(n, x2[1], x2[2])

  params_d <- data.frame(tau_d, sigma_x_d, t_switch_d, x1_d, x2_d)
  print(params_d)

  #### Generating fake data ####
  x_res <- data.frame(matrix(NA, nrow = t_d, ncol = n))
  for (j in 1:n) {
    x_d <- vector(mode = "numeric", length = t_d)

    x_d[1] <- rnorm(1, x1_d[j], sigma_x_d[j])
    for (i in 2:t_d) {
      t1_d <- t_switch_d[j] - i - 1
      t2_d <- i - t_switch_d[j]
      if (i - 1 <= t_switch_d[j]) {
        x_d[i] <- rnorm(
          1,
          (x_d[i - 1] - x1_d[j]) * exp(-1 / tau_d[j]) + x1_d[j],
          sigma_x_d[j] * sqrt(1 - exp(-2 / tau_d[j]))
        )
      } else if (i - 2 >= t_switch_d[j]) {
        x_d[i] <- rnorm(
          1,
          (x_d[i - 1] - x2_d[j]) * exp(-1 / tau_d[j]) + x2_d[j],
          sigma_x_d[j] * sqrt(1 - exp(-2 / tau_d[j]))
        )
      } else {
        t1_d <- t_switch_d[j] - (i - 2)
        t2_d <- (i - 1) - t_switch_d[j]
        x_d[i] <- rnorm(
          1,
          ((((x_d[i - 1] - x1_d[j]) * exp(-t1_d / tau_d[j]) + x1_d[j]) - x2_d[j]) * exp(-t2_d / tau_d[j])) + x2_d[j],
          sigma_x_d[j] * sqrt((1 - exp(-2 / tau_d[j])))
        )
      }
    }
    x_res[j] <- x_d
  }
  return(x_res)
}