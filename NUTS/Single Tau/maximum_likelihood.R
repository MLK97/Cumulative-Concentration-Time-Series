#############################################################
# R script for calculating the maximum likelihood that
# can be found and was actually found.
#############################################################

maximum_log_likelihood <- function(time, x_hat, t_switch, tau, sigma_x, x) {
    # time:<numeric>: Number of timesteps x_i
    # x_hat:<vector:numeric>: Values of equilibrium for time interval 1, 2 ...
    # t_switch:<vector:numeric>: Timevalues at which it will switch to the next x_hat equilibrium value
    # tau:<numeric>: Time required to reach x_hat of current time interval
    # sigma_x:<numeric>: Fluctuation
    # x:<vector:numeric>: Data generated from Stan or DummyData program and position values of tension graph

    log_likelihood <- 0 # Value of the log likelihood
    m <- 1 # Keeps track of the corresponding x_hat (and tau) values of the current time interval
    res <- vector(mode = "numeric", length = time)

    log_likelihood <- log_likelihood + dnorm(
        x[1],
        x_hat[m],
        sigma_x
    )
    res[1] <- log_likelihood

    for (i in 2:time) {
        if (i - 1 <= t_switch[m]) {
            log_likelihood <- log_likelihood + dnorm(
                x[i],
                (x[i - 1] - x_hat[m]) * exp(-1 / tau) + x_hat[m],
                sigma_x * sqrt(1 - exp(-2 / tau))
            )
        } else {
            t1 <- t_switch[m] - (i - 2)
            t2 <- (i - 1) - t_switch[m]
            log_likelihood <- log_likelihood + dnorm(
                x[i],
                ((((x[i - 1] - x_hat[m]) * exp(-t1 / tau) + x_hat[m]) -
                    x_hat[m + 1]) * exp(-t2 / tau)) + x_hat[m + 1],
                sigma_x * sqrt(1 - exp(-2 / tau))
            )
            m <- m + 1
        }
        res[i] <- log_likelihood
    }
    return(res)
}