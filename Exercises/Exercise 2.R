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

# A) use day_interval: a vector of two values: lower and upper time interval of 
# days that will be used to calculate statistical metrics. Negative values 
# indicate previous growing season days. This argument overwrites the 
# calculation limits defined by lower_limit and upper_limit arguments.

ds1 <- daily_response(response = chron_res, 
                          env_data = temperature,
                          method = "cor",
                          row_names_subset = TRUE,
                          remove_insignificant = FALSE)

plot(ds1, type = 2)
summary(ds1)

ds2 <- daily_response(response = chron_res, 
                      env_data = temperature,
                      method = "cor",
                      day_interval = c(120, 300),
                      row_names_subset = TRUE,
                      remove_insignificant = FALSE)

library(ggpubr)
ggarrange(plot(ds1, type = 2), 
          plot(ds2, type = 2),
          ncol = 1)


# Hint: use scale_fill_gradient2
ggarrange(plot(ds1, type = 2) + scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)), 
          plot(ds2, type = 2) + scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)),
          ncol = 1)


# including previous year
ds3 <- daily_response(response = chron_res, 
                      env_data = temperature,
                      method = "cor",
                      previous_year = TRUE,
                      row_names_subset = TRUE,
                      remove_insignificant = FALSE)

plot(ds3, type = 2)

ds4 <- daily_response(response = chron_res, 
                      env_data = temperature,
                      method = "cor",
                      day_interval = c(-120, 300), # minus sign indicates previous year time
                      previous_year = TRUE,
                      row_names_subset = TRUE,
                      remove_insignificant = FALSE)

ggarrange(plot(ds3, type = 2) + scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)), 
          plot(ds4, type = 2) + scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)),
          ncol = 1)


