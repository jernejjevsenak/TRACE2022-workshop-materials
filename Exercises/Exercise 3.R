library("dendroTools")
library("dplR")

data <- read.rwl("rwls/czec010.rwl") 
data_detrended <- detrend(data, method = "Spline", nyrs = 32)
data_chron <- chron(data_detrended, prewhiten = TRUE)
chron_res <- data_chron[,2, drop = FALSE]

# Precipitation data, source E-OBS
precipitation <- read.table("climate/Psum_czec010.txt", header = TRUE)
precipitation$lat <- NULL
precipitation$lon <- NULL
precipitation <- data_transform(precipitation, format = 'daily')

temperature <- read.table("climate/Tavg_czec010.txt", header = TRUE)
temperature$lat <- NULL
temperature$lon <- NULL
temperature <- data_transform(temperature, format = 'daily')

# Detrend climate
# dc_method : a character string to determine the method to detrend climate 
# (environmental) data. Possible values are c("Spline", "ModNegExp", "Mean", 
# "Friedman", "ModHugershoff"). Defaults to "none". See dplR R package for more
# details.

# dc_nyrs	a number giving the rigidity of the smoothing spline, defaults to 0.67
# of series length if nyrs is NULL.

# dc_f	a number between 0 and 1 giving the frequency response or wavelength 
# cutoff. Defaults to 0.5 (see dplR R package).

# dc_pos.slope a logical flag. Will allow for a positive slope to be used in 
# method "ModNegExp" and "ModHugershoff". If FALSE the line will be horizontal.

# dc_constrain.nls a character string which controls the constraints of the 
# "ModNegExp" model and the "ModHugershoff" (see dplR R package)

# dc_span	a numeric value controlling method "Friedman", or "cv" (default) for 
# automatic choice by cross-validation (see dplR R package).

# dc_bass	a numeric value controlling the smoothness of the fitted curve in 
# method "Friedman" (see dplR R package).

# dc_difference	a logical flag. Compute residuals by subtraction if TRUE, 
# otherwise use division (see dplR R package).

ds1 <- daily_response(response = chron_res, 
                          env_data = temperature,
                          method = "cor",
                          lower_limit = 21, upper_limit = 35,
                          row_names_subset = TRUE,
                          remove_insignificant = FALSE)

ds2 <- daily_response(response = chron_res, 
                      env_data = temperature,
                      method = "cor",
                      lower_limit = 21, upper_limit = 35,
                      dc_method = "Spline",
                      row_names_subset = TRUE,
                      remove_insignificant = FALSE)

library(ggpubr)
ggarrange(plot(ds1, type = 2) + ggtitle('Without detrending') + scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)), 
          plot(ds2, type = 2) + ggtitle('With detrending') + scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)))

summary(ds1)
summary(ds2)
