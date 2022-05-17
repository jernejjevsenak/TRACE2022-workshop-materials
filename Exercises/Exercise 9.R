library(dplyr)
library(ggplot2)
library(rTG)
library(MLmetrics)

# Open data from rTG
data(data_trees)
data(parameters)

# run simulation 1 
simulation_1 <- XPSgrowth(data_trees, parameters,
                          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
                          fitting_method = c("brnn", "gam", "gompertz"),
                          search_initial_gom = TRUE,
                          fitted_save = TRUE, add_zeros = TRUE,
                          add_zeros_before = 75, post_process = TRUE)

plot(simulation_1)
summary(simulation_1)

# save simulation_1 output in separate data frame and remove manually added zeros
df_A <- simulation_1[[1]]
df_A <- mutate(df_A, width = ifelse(note == "added zero", NA, width))


# Improve labels for plotting purposes
df_A <- mutate(df_A, method = ifelse(method == "brnn", "BRNN", 
                                     ifelse(method == "gam", "GAM", 
                                            ifelse(method == "gompertz", "Gompertz", NA))))

# Sort conifer and broadlevaes
df_A$Species <- factor(df_A$Species, levels = c("FASY", "QUPU", "PCAB"))

# Create and save ggplot object
ggplot(df_A, aes(x = doy, y = width_pred, col = factor(Tree), group = key)) + geom_line() +
  geom_point(df_A, mapping = aes(x = doy, y = width),alpha = 0.2) +
  facet_grid(Species + Tissue ~ method, scales = "free") +
  scale_colour_brewer(palette = "Dark2") + theme_minimal() +
  theme(legend.position = "none", text = element_text(size=18)) +
  ylab("Width [µm] / Cell Number") + xlab("Day of Year") +
  theme(panel.background = element_rect(fill = "gray97",
                                        colour = "gray80",
                                        size = 0.5, linetype = "solid")) 

# Calculate RMSE, %RMSE and R2
group_by(dplyr::filter(df_A, !is.na(width)), Tissue, method, Species) %>%
                 summarise(RMSE = RMSE(width_pred, width),
                           rRMSE = RMSE/mean(width)*100,
                           R2_Score = MLmetrics::R2_Score(width, width_pred)) %>%
                 arrange(-RMSE) %>% arrange(Tissue, Species, method)
