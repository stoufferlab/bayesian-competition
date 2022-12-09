library(greta)
x <- iris$Petal.Length
y <- as_data(iris$Sepal.Length)
int <- normal(0, 10)
coef <- normal(0, 10)
sd <- student(3, 0, 10)
mean <- int + coef * x
distribution(y) <- normal(mean, sd)
m <- model(int, coef, sd)
draws <- mcmc(m, n_samples = 1000, chains = 4)
# summary(draws)
# generate predictions from prior draws
y_hat_prior <- calculate(y, nsim = 100)$y[, , 1]
# generate predictions from posterior draws
y_hat_posterior <- calculate(y, values = draws, nsim = 100)$y[, , 1]
# visualise y_hat generated from prior vs. posterior
# using row median for now
y_hat_prior_median     <- apply(y_hat_prior, 2, median)
y_hat_posterior_median <- apply(y_hat_posterior, 2, median)
# plot(y_hat_prior_median, y_hat_posterior_median)
plot(density(y_hat_prior_median), ylim = c(0, 0.7), col = "red", 
     main = "Comparing prior and posterior prediction")
lines(density(y_hat_posterior_median), add = TRUE, col = "blue")