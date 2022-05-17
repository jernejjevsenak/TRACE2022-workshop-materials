library(dendroTools)
library(dplR)
library(ggplot2)
library(reshape2)

# Create a list of all rwl chronologies
list_rwl <- list.files(path = "rwls", pattern = ".rwl")

n_years <- 30 # years
selected_fixed_window <- 45 # days

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

start_years <- seq(1950, max(as.numeric(row.names(chron_res))) - n_years + 1)

for (t_stab in start_years){
  
# temperature
temp_dsT <- daily_response(response = chron_res, env_data = temperature,
                           method = "cor",
                           row_names_subset = TRUE, 
                           remove_insignificant = FALSE,
                           
                           fixed_width = selected_fixed_window,
                           subset_years = c(t_stab, t_stab + n_years - 1)
                           )
  
  result_temp_T <- data.frame(temp_dsT$calculations)
  result_temp_T$climate <- "Temperature"
  result_temp_T$s_year <- t_stab
  result_temp_T$site <- site
  list_results[[b]] <- result_temp_T
  b = b + 1

# precipitation
temp_dsP <- daily_response(response = chron_res, env_data = precipitation,
                             method = "cor",
                             row_names_subset = TRUE, 
                             remove_insignificant = FALSE,
                             
                             fixed_width = selected_fixed_window,
                             subset_years = c(t_stab, t_stab + n_years - 1)
  )
  
  result_temp_P <- data.frame(temp_dsP$calculations)
  result_temp_P$climate <- "Precipitation"
  result_temp_P$s_year <- t_stab
  result_temp_P$site <- site
  list_results[[b]] <- result_temp_P
  b = b + 1
  
}

}

binded <- do.call(rbind, list_results)
binded <- melt(binded, id.vars = c("climate", "s_year", "site"))

# remove correlations with low significance
critical_r <- function(n, alpha = .05) {
  df <- n - 2
  critical_t <- qt(alpha / 2, df, lower.tail = FALSE)
  critical_r_value <- sqrt((critical_t ^ 2) / ((critical_t ^ 2) + df))
  return(critical_r_value)
}

binded$value <- ifelse(abs(binded$value) < critical_r(n_years, alpha = 0.05), NA, binded$value)
binded$period <- paste0(binded$s_year, " - ", binded$s_year + n_years - 1)

summary(binded$value)

ggplot(binded,aes_(x = ~as.numeric(variable), y = ~period, fill = ~value)) +
  geom_tile() +
  facet_grid(site  ~  climate) +
  xlab("Day of Year") +
  ylab("Period") +
  scale_x_continuous(expand=c(0,0), breaks = c(75,150, 225, 300), labels = c("March", "June", "August", "November")) +
  scale_y_discrete(expand=c(0,0),
                   breaks = unique(binded$period)[seq(1,40, by = 3)]) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", na.value = 'gray97', midpoint = 0, limits = c(-0.60, 0.60)) +
  
  theme_minimal() +  theme(axis.text = element_text(size = 10),
                           strip.background = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           legend.title = element_blank(),
                           legend.position = "bottom", legend.key.width = unit(3, "line"),
                           panel.background = element_rect(fill = "gray97",
                                                           colour = "gray80",
                                                           size = 0.5, linetype = "solid")) 

ggsave("temporal_stabilty_combined.png", height = 20, width = 15)




