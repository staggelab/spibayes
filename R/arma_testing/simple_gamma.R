

require(rstan)
require(tidyverse)

n <- 500
rand_y <- rgamma(n, shape = 0.5, scale = 5)

fit_data <- list(N = n, y = rand_y)


fit <- stan(file = 'my_attempt_solver2.stan', data = fit_data,  iter = 2000, init =list(list(shape = 0.5, scale = 10)), chains = 1)

print(fit)
plot(fit)


est_df <- data.frame(shape =  rstan::extract(fit, "shape")$shape, scale = rstan::extract(fit, "scale")$scale, qgamma_draw = rstan::extract(fit, "qgamma_result")$qgamma_result)
est_df <- est_df %>% mutate(qgamma_correct = qgamma(0.4, shape = shape, scale = scale))
head(est_df)

ggplot(est_df, aes(x=shape)) + geom_density() + geom_vline(xintercept = 0.5, colour = "red") + theme_classic()
ggplot(est_df, aes(x=scale)) + geom_density() + geom_vline(xintercept = 5, colour = "red") + theme_classic()
ggplot(est_df, aes(x=qgamma_result)) + geom_density() + geom_vline(xintercept = qgamma(0.4, shape = 0.5, scale = 5), colour = "red") + theme_classic()


