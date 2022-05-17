library(rTG)
library(dplyr)
library(ggplot2)

# open dendrometer data
data("data_dendrometers")

# apply the XPSgrowth() with the post-processing algorithm
sim1 <- XPSgrowth(data_dendrometers, unified_parameters = TRUE,
                  ID_vars = c("site", "species", "year", "tree"),
                  fitting_method = c("brnn", "gam", "gompertz"), 
                  brnn_neurons = 2, gam_k = 9, gam_sp = 0.5, search_initial_gom = TRUE,
                  add_zeros = FALSE, post_process = TRUE)

# extract the fitted values and add label "With post-process"
simulation_A <- sim1$fitted
simulation_A$label <- "With post-process"

# apply the XPSgrowth() without the post-processing algorithm
sim2 <- XPSgrowth(data_dendrometers, unified_parameters = TRUE,
                  ID_vars = c("site", "species", "year", "tree"),
                  fitting_method = c("brnn", "gam", "gompertz"), 
                  brnn_neurons = 2, gam_k = 9, gam_sp = 0.5, search_initial_gom = TRUE,
                  add_zeros = FALSE, post_process = FALSE)

# extract the fitted values and add label "Without post-process"
simulation_B <- sim2$fitted
simulation_B$label <- "Without post-process"

# rbind both simulations and improve method labels
combined <- rbind(simulation_A, simulation_B)
combined <- dplyr::mutate(combined, method = ifelse(method == "brnn", "BRNN", 
                                                    ifelse(method == "gam", "GAM", 
                                                           ifelse(method == "gompertz", "Gompertz", NA))))
# create ggplot
ggplot(combined, aes(x = doy, y = width)) + geom_point() +
  geom_line(aes(x = doy, y = width_pred), col = "red") +
  facet_grid(species + label ~ method) + scale_colour_brewer(palette = "Dark2") + theme_minimal() +
  theme(legend.position = "none", text = element_text(size=18)) +
  ylab("Width") + xlab("Day of Year") +
  theme(panel.background = element_rect(fill = "gray97",
                                        colour = "gray80",
                                        size = 0.5, linetype = "solid"))

ggsave("dendrometer_example.png", height = 9, width = 10)

