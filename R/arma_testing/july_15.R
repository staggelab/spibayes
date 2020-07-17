

###########################################################################
###  Load functions
###########################################################################
require(tidyverse)
require(here)

### 
require(rstan)

### To save in SVG
require(svglite)
require(viridis)

### Packages for spi
require(fitdistrplus)
require(lubridate)

select <- dplyr::select

theme_set(theme_classic(8))


###########################################################################
## Set initial values
###########################################################################
### Set seed so everyone gets the same random numbers (reproducible example)
set.seed(7890)

### SPI-3 3 months equals 92 days
n_roll <- 92



###########################################################################
###  Create a true SPI time series
###########################################################################
### Assume 100 years from 1920 to 2019, will then take subsets of this
spi_true <- data.frame(date = seq(as.Date("1920-01-01"), as.Date("2019-12-31"), by = "1 day")) %>%
	mutate(jdate = yday(date)) %>%
	mutate(month_day = paste0(month(date),"-",day(date))) %>%
	mutate(year = year(date))

n_days <- dim(spi_true)[1]

### Use a 92 day moving average MA(91) with innovations of sqrt(92)/n_roll, which produces N(0,1)
### Technically, the moving average would be on precip, not SPI, but this is a very close approximation
innov_c <- rnorm(n_days, 0, sqrt(n_roll))

### Generate SPI series using an MA(91) model with coef of 1. Need to remove the first 91 values
sim_spi <- arima.sim(list(order = c(0,0,(n_roll - 1)), ma = rep(1,(n_roll -1))), n = n_days, innov=innov_c/n_roll)
sim_spi[seq(1,(n_roll-1))] <- NA

spi_true <- spi_true %>%
	mutate(innov = innov_c/n_roll) %>%
	mutate(spi = c(sim_spi))

ggplot(spi_true, aes(x=date, y=spi)) + geom_line() + theme_classic()

mean(spi_true$spi, na.rm=TRUE)
sd(spi_true$spi, na.rm=TRUE)

###########################################################################
###  Make up synthetic gamma distribution (stationary)
###########################################################################
### Make some daily synthetic distributions that vary smoothly through time
true_param_stat <- tibble(jdate = seq(1,365), shape = rep(1.7, 365)) %>%
	mutate(mean = rep(4, 365)) %>%
	mutate(scale = mean / shape) %>%
	mutate(rate = 1/scale) %>%
	mutate(disp = 1/shape) %>%
	select(jdate, shape, scale, scale, rate, mean, disp)


###########################################################################
###  Sample from this distribution
###########################################################################
### Apply day 365 to day 366, join with true parameters, and then return to the original format
synth_stat_df <- spi_true %>%
	mutate(jdate = case_when(
			jdate == 366 ~ 365,
			TRUE ~ jdate)
	) %>%
	left_join(true_param_stat) %>%
	mutate(jdate = yday(date)) 

### First transform SPI into percentile. Then set zeros when percentile is below theta threshold
### Stretch the remaining probability space so that it fills the space (0 to 1)
### Then sample from gamma distribution and combine with zeros
synth_stat_df <- synth_stat_df %>% 
	mutate(p_overall = pnorm(spi)) %>%
	mutate(precip = qgamma(p_overall, shape = shape, scale = scale))



###########################################################################
###  Prepare data for fitting
###########################################################################
fitting_df <- synth_stat_df %>%
	drop_na(precip) %>%
	filter(year > 2009)



fit_data <- list(N = nrow(fitting_df), y = fitting_df$precip)
str(fit_data)


fit_simple <- stan(file = 'gamma2.stan', data = fit_data,  init = init_vals, iter = 200, chains = 1)

est_simple <- data.frame(shape =  rstan::extract(fit_simple, "shape")$shape, mean = rstan::extract(fit_simple, "mean_param")$mean_param) %>%
	mutate(scale = mean/shape)


est_simple <- data.frame(shape =  rstan::extract(fit_simple, "shape")$shape, rate = rstan::extract(fit_simple, "rate")$rate) %>%
	mutate(scale = 1/rate)


init_vals <- list(list(shape = 1.5, rate = 0.5, sigma_error = 0.1))
##fit <- stan(file = 'gamma_ar.stan', data = fit_data,  iter = 100, chains = 1)

fit <- stan(file = 'gamma_ar3.stan', data = fit_data,  init = init_vals, iter = 2000, chains = 1, control = list(metric = "dense_e"), algorithm = "HMC")




print(fit)
plot(fit)


est_df <- data.frame(shape =  rstan::extract(fit, "shape")$shape, rate = rstan::extract(fit, "rate")$rate) %>%
	mutate(scale = 1/rate)


ggplot(est_df, aes(x=shape)) + geom_density() + geom_vline(xintercept = synth_stat_df$shape[1], colour="red") + geom_density(data = est_simple, fill = "blue", alpha =0.5)

ggplot(est_df, aes(x=scale)) + geom_density() + geom_vline(xintercept = synth_stat_df$scale[1], colour="red")


yoop <- rstan::extract(fit, "mu")

plot(fitting_df$spi, type="l")
lines(yoop[[1]][4,], col="red")





, qgamma_draw = rstan::extract(fit, "qgamma_result")$qgamma_result)
est_df <- est_df %>% mutate(qgamma_correct = qgamma(0.4, shape = shape, scale = scale))
head(est_df)

ggplot(est_df, aes(x=shape)) + geom_density() + geom_vline(xintercept = 0.5, colour = "red") + theme_classic()
ggplot(est_df, aes(x=scale)) + geom_density() + geom_vline(xintercept = 5, colour = "red") + theme_classic()
ggplot(est_df, aes(x=qgamma_result)) + geom_density() + geom_vline(xintercept = qgamma(0.4, shape = 0.5, scale = 5), colour = "red") + theme_classic()












###############################################
###  Try with innovations
###############################################

fit_data <- list(N = nrow(fitting_df), y = fitting_df$precip, n_roll = 92)
str(fit_data)





init_vals <- list(list(shape = 1.5, rate = 0.5, sigma_error = 0.1))
##fit <- stan(file = 'gamma_ar.stan', data = fit_data,  iter = 100, chains = 1)

fit <- stan(file = 'gamma_ar2c.stan', data = fit_data,  init = init_vals, iter = 600, chains = 1, control = list(metric = "dense_e"), algorithm = "HMC")

print(fit)
plot(fit)

fit <- stan(file = 'gamma_ar2d.stan', data = fit_data,  init = init_vals, iter = 200, chains = 1, control = list(metric = "dense_e"), algorithm = "HMC")



###############################################
###  Try with covariance matrix
###############################################





fit_data <- list(N = nrow(fitting_df), y = fitting_df$precip)
str(fit_data)

require(hierband)
fit_data$Sigma <- ma(fit_data$N, n_roll - 1)$Sig 



init_vals <- list(list(shape = 1.5, rate = 0.5, mu = fitting_df$spi, sigma_error = 0.1))
##fit <- stan(file = 'gamma_ar.stan', data = fit_data,  iter = 100, chains = 1)

fit <- stan(file = 'gamma_ar4.stan', data = fit_data,  init = init_vals, iter = 100, chains = 1)





require(hierband)
sigma_mat <- ma(1000, 92)$Sig 


yup <- mvrnorm(n = 100, mu = rep(0, 1000) , Sigma = sigma_mat)





