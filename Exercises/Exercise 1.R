# install.packages("dendroTools") - the most recent version?

library("dendroTools")
library("dplR")

getwd()

# source https://doi.org/10.1016/j.dendro.2021.125845
# Picea abies
data <- read.rwl("rwls/czec013.rwl")
data_detrended <- detrend(data, method = "Spline", nyrs = 32, f = 0.5)
data_chron <- chron(data_detrended, prewhiten = TRUE)
chron_res <- data_chron[,2, drop = FALSE]

# Precipitation data, source E-OBS
precipitation <- read.table("climate/Psum_czec013.txt", header = TRUE)
precipitation$lat <- NULL
precipitation$lon <- NULL
precipitation <- data_transform(precipitation, format = 'daily')

temperature <- read.table("climate/Tavg_czec013.txt", header = TRUE)
temperature$lat <- NULL
temperature$lon <- NULL
temperature <- data_transform(temperature, format = 'daily')

glimpse_daily_data(precipitation)
glimpse_daily_data(temperature)

# 1 Precipitation correlations
ds_prec <- daily_response(response = chron_res, 
               env_data = precipitation,
               method = "cor",
               remove_insignificant = TRUE,
               alpha = 0.01,
               row_names_subset = TRUE)

plot(ds_prec, type = 2)
summary(ds_temp)

# 2 Temperature correlations
ds_temp <- daily_response(response = chron_res, 
                          env_data = temperature,
                          method = "cor",
                          remove_insignificant = TRUE,
                          alpha = 0.01,
                          row_names_subset = TRUE)

plot(ds_temp, type = 2)

summary(ds_temp)


