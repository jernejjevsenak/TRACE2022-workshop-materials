library(dendroTools)
library(dplR)

# Create a list of all rwl chronologies
list_rwl <- list.files(path = "rwls", pattern = ".rwl")

# Create empty list where you later save correlations
list_calculations <- list()
b <- 1

for (i in list_rwl){
  
  print(paste0("I currently work with ", i))

  # Open temporal rwl
  data_rwl <- read.rwl(paste0("rwls/", i)) 
  data_detrended <- detrend(data_rwl, method = "Spline", nyrs = 32)
  data_chron <- chron(data_detrended, prewhiten = TRUE)
  chron_res <- data_chron[,2, drop = FALSE]
  
  chron_res
  
  # Extract site name
  site <- gsub(".rwl", "", i)

  # Open temporal precipitation data
  precipitation <- read.table(paste0("climate/Psum_", site, ".txt"), header = TRUE)
  precipitation$lat <- NULL
  precipitation$lon <- NULL
  precipitation <- data_transform(precipitation, format = 'daily')

  # Open temporal temperature data  
  temperature <- read.table(paste0("climate/Tavg_", site,".txt"), header = TRUE)
  temperature$lat <- NULL
  temperature$lon <- NULL
  temperature <- data_transform(temperature, format = 'daily')
  
  # Calculate daily response with precipitation data
  dsP <- daily_response(response = chron_res, 
                        env_data = precipitation,
                        method = "cor",
                        lower_limit = 21, upper_limit = 35,
                        row_names_subset = TRUE,
                        remove_insignificant = TRUE)
  
  # extract daily correlations (precipitation)
  temp_calculations <- data.frame(dsP$calculations)
  temp_calculations$climate <- "Precipitation"
  temp_calculations$site <- site
  temp_calculations$season <- seq(21, 35)
  
  list_calculations[[b]] <- temp_calculations
  b = b +1
  
  # Calculate daily response with temperature data
  dsT <- daily_response(response = chron_res, 
                        env_data = temperature,
                        method = "cor",
                        lower_limit = 21, upper_limit = 35,
                        row_names_subset = TRUE,
                        remove_insignificant = TRUE)
  
  # extract daily correlations (temperature)
  temp_calculations <- data.frame(dsT$calculations)
  temp_calculations$climate <- "Temperature"
  temp_calculations$site <- site
  temp_calculations$season <- seq(21, 35)
  
  list_calculations[[b]] <- temp_calculations
  b = b +1
  
}

# rbind all list elements into one data frame
binded <- do.call(rbind, list_calculations)
binded[,c(340:349)]



# transform data into long format using the melt function from reshape2 r package
library(reshape2)
binded <- melt(binded, id.vars = c("site", "season", "climate"))

library(ggplot2)
ggplot(binded,aes_(x = ~as.numeric(variable), y = ~season, fill = ~value)) +
  geom_tile() +
  facet_grid(site  ~  climate, scales = "free") +
  xlab("Day of Year") +
  ylab("Season Length") +
  scale_x_continuous(expand=c(0,0), breaks = c(75,150, 225, 300), labels = c("March", "June", "August", "November")) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.5, 0.5)) +
    theme_minimal() +  theme(
                           axis.text = element_text(size = 10),
                           axis.title = element_blank(),
                           strip.background = element_blank(),
                           plot.title = element_text(size = 16),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom", legend.key.width = unit(3, "line"),
                           panel.background = element_rect(fill = "gray97", colour = "gray80",
                                                           size = 0.5, linetype = "solid")) 

ggsave("Exercise4.png")
