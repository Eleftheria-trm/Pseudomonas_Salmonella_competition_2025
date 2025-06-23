# Load required packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(emmeans)
library(lme4)
library(ggtext)

# Define function to read and process data from each sheet
read_and_process_data <- function(file_path, sheet_name) {
  data <- read_excel(file_path, sheet = sheet_name)
  data <- data %>%
    mutate(
      Strain = as.factor(Strain),
      `In competition` = as.factor(`In competition`),
      `Cell counted` = as.factor(`Cell counted`),
      `Time (hrs)` = as.numeric(`Time (hrs)`),
      tech.rep = as.factor(tech.rep),
      cfu.seedling = as.numeric(cfu.seedling)
    ) %>%
    mutate(Experiment = case_when(
      sheet_name == "Experiment1" ~ "Co-inoculation",
      sheet_name == "Experiment2" ~ "48hrs pre-incubation with BAHP",
      sheet_name == "Experiment3" ~ "48hrs pre-incubation with Salmonella",
      TRUE ~ sheet_name
    ))
  return(data)
}

# Read and process data from all sheets
file_path <- "All_data_Sal_treated.xlsx"
sheet_names <- c("Experiment1", "Experiment2", "Experiment3")

all_data <- lapply(sheet_names, function(sheet) {
  read_and_process_data(file_path, sheet)
}) %>%
  bind_rows()

# Convert Experiment to a factor with specified levels
all_data$Experiment <- factor(all_data$Experiment, 
                              levels = c("Co-inoculation", 
                                         "48hrs pre-incubation with BAHP", 
                                         "48hrs pre-incubation with Salmonella"))

# Replace "Pseudomonas" with "BAHP" in 'Cell counted'
all_data$`Cell counted` <- as.character(all_data$`Cell counted`)
all_data$`Cell counted`[all_data$`Cell counted` == "Pseudomonas"] <- "BAHP"
all_data$`Cell counted` <- as.factor(all_data$`Cell counted`)

# Summarize CFU/ml Data
summary_data <- all_data %>%
  group_by(`Time (hrs)`, `Cell counted`, Experiment) %>%  
  summarise(
    mean_cfu = mean(cfu.seedling, na.rm = TRUE),
    se_cfu = sd(cfu.seedling, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Visualization: CFU/ml Over Time
ggplot() +
  geom_point(data = all_data,
             aes(x = `Time (hrs)`, y = cfu.seedling, color = `Cell counted`),
             alpha = 0.3, position = position_jitter(width = 0.5), size = 3) +
  geom_line(data = summary_data,
            aes(x = `Time (hrs)`, y = mean_cfu, color = `Cell counted`),
            size = 1) +
  geom_ribbon(data = summary_data,
              aes(x = `Time (hrs)`, ymin = mean_cfu - se_cfu,
                  ymax = mean_cfu + se_cfu, fill = `Cell counted`),
              alpha = 0.2, colour = NA) +
  labs(title = "CFU/Seedling Over Time",
       x = "Time (hrs)",
       y = "CFU/seedling") +
  scale_y_log10() +
  facet_wrap(~ Experiment, nrow = 1) +
  theme_minimal() +
  theme(legend.text = element_markdown(size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.text = element_text(size = 20)) +
  scale_color_discrete(
    name = "Species",
    labels = c("BA<sub>HP</sub>", "<i>Salmonella</i>")
  ) +
  scale_fill_discrete(
    name = "Species",
    labels = c("BA<sub>HP</sub>", "<i>Salmonella</i>")
  )

# Save the plot
ggsave("cfu_ml_over_time.png", width = 22, height = 7, dpi = 300)

# This script was generated with the assistance of Perplexity AI.
# All code has been reviewed and tested for accuracy and functionality by the project team.