library(rTG)
library(dplyr)
library(ggplot2)
library(MLmetrics)

# Open data from rTG
data(data_trees)
data(parameters)

# Select one tree for XYLEM
data_trees <- dplyr::filter(data_trees, Species == "QUPU", Tissue == "XYLEM", Tree == 2)

# Simulation with no correction
simulation_A <- XPSgrowth(data_trees, parameters,
                ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                fitting_method = c("brnn", "gam"),
                search_initial_gom = TRUE, add_zeros_before = 75,
                fitted_save = FALSE, unified_parameters = FALSE,
                add_zeros = FALSE, post_process = FALSE)

# Simulation with added zeros
simulation_B <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam"),
                          search_initial_gom = TRUE, add_zeros_before = 75,
                          fitted_save = FALSE, unified_parameters = FALSE,
                          add_zeros = TRUE, post_process = FALSE)

# Simulation with added zeros and post-correction
simulation_C <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam"),
                          search_initial_gom = TRUE, add_zeros_before = 75,
                          fitted_save = FALSE, unified_parameters = FALSE,
                          add_zeros = TRUE, post_process = TRUE)

# Extract the first elements from each simulation and give description
df_a <- simulation_A[[1]]; df_a$app <- "No correction" 
df_b <- simulation_B[[1]]; df_b$app <- "Added zeros"
df_c <- simulation_C[[1]]; df_c$app <- "Added zeros and\npost-correction"

# The next one is only to produce labels
df_abc <- simulation_B[[1]]; df_abc$method <- ifelse(df_abc$method == "gam", "GAM", "BRNN")

# Rbind a, b and c data frames and improve labels
temp_data <- rbind( df_a, df_b, df_c)
temp_data$app <- factor(temp_data$app, levels = c("No correction", "Added zeros", "Added zeros and\npost-correction"))
temp_data$method <- ifelse(temp_data$method == "gam", "GAM", "BRNN")

# create a data frame with labels 
label_df <- data.frame(ID = c("(a)","(b)","(c)","(d)","(e)","(f)"), 
                       method = c("BRNN", "BRNN","BRNN", "GAM","GAM", "GAM"), 
                       app = c("No correction", "Added zeros", "Added zeros and\npost-correction",
                               "No correction", "Added zeros", "Added zeros and\npost-correction"),
                       doy = 1, width_pred = 490)

label_df$app <- factor(label_df$app, levels = c("No correction", "Added zeros", "Added zeros and\npost-correction"))

# Create ggplot object
ggplot(temp_data, aes(x = doy, y = width_pred)) +
    geom_line() + facet_grid(app~method) + theme_bw() +
    geom_point(df_abc, mapping = aes(x = doy, y = width, alpha = note)) +
    scale_alpha_discrete(range=c(0.1, 1)) + guides(alpha = "none") +
    ylab("Xylem Width [µm]") + xlab("Day of Year") + theme_minimal() +
    theme(panel.background = element_rect(fill = "gray97",
    colour = "gray80", size = 0.5, linetype = "solid"),
    text = element_text(size = 14)) +
    geom_text(data = label_df, mapping = aes(label = ID))
