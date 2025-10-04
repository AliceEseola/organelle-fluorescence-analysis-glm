# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# Read data
df <- read_csv("~/Desktop/sp8/organelle7.csv")

# Summary of raw data
summary(df)

# Rename region_index to Compartment for clarity
df <- df %>% rename(Compartment = region_index)

# Group and summarize data
organelle7_df <- df %>%
  group_by(Reps, Compartment, Time, Marker) %>%
  summarize(mean_fluor = mean(Norm_fluoro), .groups = "drop")

# Save summarized data
write.csv(organelle7_df, file = "~/Desktop/sp8/organelles7_df.csv", row.names = FALSE)

# Convert columns to factors and create numeric time
organelle7_df <- organelle7_df %>%
  mutate(
    Compartment = as.factor(Compartment),
    Time = factor(Time, levels = c("0h", "4h", "8h", "24h"), ordered = TRUE),
    Marker = factor(Marker, levels = c("Scad2", "Pex6", "Sec7", "Grh1", "Vac8", "Rpl25", "HDEL"), ordered = TRUE),
    time_c = as.numeric(str_replace(Time, "h", ""))
  )

# Initial plot with jitter and GLM smoothing
ggplot(organelle7_df, aes(x = Time, y = mean_fluor, colour = Compartment)) +
  geom_jitter() +
  geom_smooth(method = "glm") +
  facet_grid(~Marker)

# Final publication-style plot
ggplot(organelle7_df, aes(x = time_c, y = mean_fluor, colour = Compartment)) +
  geom_smooth(method = "glm") +
  facet_wrap(~Marker, ncol = 2) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 16, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 16, face = "bold", colour = "black", vjust = 1.3),
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.title.x = element_text(size = 16, face = "bold", colour = "black", vjust = 0.8),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.position = c(0.8, 0.11),
    strip.background = element_rect(colour = "black", fill = "white")
  ) +
  scale_y_continuous(name = "Normalised GFP intensity") +
  scale_x_continuous(name = "Time (hour)", limits = c(0, 24), breaks = seq(0, 24, by = 4))

# Save plot
ggsave("organelle7_df.tiff", dpi = 900, height = 297, width = 210, units = "mm")