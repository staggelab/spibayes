# spibayes package


### To install
Package is not on cran, but you can install directly from github.

```{r}
library(devtools)
install_github("staggelab/spibayes")
```

### Example showing how to use the package to fit the tensor
In order to fit a new SPI distribution, first you must preprocess the data. This will give initial estimates based on an MLE fitting of a cyclic spline, i.e no change through the years.


```{r}
library(spibayes)
library(tidyverse)

### Set the output folder. You can modify this
write_output_path <- getwd()

### Precip data
precip_df <- expand.grid(jdate = seq(1,365), year = seq(1940, 2020,1 ))
precip_df <- precip_df %>%
	mutate(precip = rgamma(dim(precip_df)[1], shape = 3, scale = 5/3))

### Set the knots
knot_loc <- list(jdate = seq(1,365, 8), year = seq(1940, 2020, 15))

### Run preprocess command
tensor_init <- pre_proc(data = precip_df, type = "tensor", knot_loc = knot_loc)

### Run the non-Bayesian code
tensor_optimize <- model_fit(spi_input = tensor_init, iter = 2000, engine = "optimize", output_dir = write_output_path)

### Estimate the distribution parameters at new locations
newdata_df <- expand.grid(jdate = seq(1,365,1), year = seq(1940, 2020, 1))
param_est_optimize <- predict_vals(tensor_optimize, newdata = newdata_df)	

names(param_est_optimize)
names(param_est_optimize$estimate)

head(param_est_optimize$estimate$gamma)
head(param_est_optimize$estimate$theta)

names(param_est_optimize$marginal)

head(param_est_optimize$marginal$mean)
head(param_est_optimize$marginal$disp)
head(param_est_optimize$marginal$theta)
```




