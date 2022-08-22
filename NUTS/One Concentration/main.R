#############################################################
# R interface script for calling and testing STAN model for
# one-concentration model
#############################################################

#### Necessary libraries ####
library(rstan)
library(GGally)
library(readxl)
library(latex2exp)

source("sampling.R")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

t <- 1000
tau <- c(600, 1)
sigma_x <- c(0.1, 0.001)
x1 <- c(0.2, 0.01)
x2 <- c(8, 2)
t_switch <- c(250, 1)

result <- result_fit(1, t, tau, sigma_x, x1, x2, t_switch)

data_dummy <- data.frame(seq(t), result[1])
data_d <- melt(data_dummy,
    id.vars = "seq.t.",
    variable.name = "Patient",
    value.name = "Concentration"
)


p <- ggplot(
    data = data_d,
    aes(x = seq.t., y = Concentration, group = Patient, colour = Patient)
)
plot <- p + geom_line() +
    xlab("Time (s)") +
    ylab(TeX("Tension $\\left(N/m\\right)$")) +
    theme(legend.position = "none")


pdf("testing.pdf")
print(plot)
dev.off()

#### Stan Model & Fit ####
model <- stan_model("General_adv.stan")
for (fake_data in result) {
    fit <- sampling(
        model,
        list(t = t, x = fake_data),
        iter = 1000, chains = 4, algorithm = "NUTS"
    )
}
print(fit)