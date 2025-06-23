library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Define function to read and process data
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
    )
  data_coincubated <- data %>%
    filter(`In competition` == "STM50")
  data_coincubated <- data_coincubated %>%
    group_by(`Time (hrs)`, tech.rep) %>%
    summarise(
      salmonella_cfu_seedling = sum(cfu.seedling[`Cell counted` == "Salmonella"], na.rm = TRUE),
      pseudomonas_cfu_seedling = sum(cfu.seedling[`Cell counted` == "Pseudomonas"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      total_cfu_seedling = salmonella_cfu_seedling +  pseudomonas_cfu_seedling,
      percentage_salmonella = (salmonella_cfu_seedling / total_cfu_seedling) * 100,
      percentage_pseudomonas = (pseudomonas_cfu_seedling / total_cfu_seedling) * 100
    )
  data_long <- data_coincubated %>%
    select(`Time (hrs)`, percentage_salmonella, percentage_pseudomonas) %>%
    pivot_longer(cols = c(percentage_salmonella, percentage_pseudomonas),
                 names_to = "Species",
                 values_to = "Percentage") %>%
    mutate(Species = ifelse(Species == "percentage_salmonella", "Salmonella", "Pseudomonas"))
  data_long <- data_long %>%
    mutate(Experiment = case_when(
      sheet_name == "Experiment1" ~ "Co-inoculation",
      sheet_name == "Experiment2" ~ "BAHP[HP]~pre-incubation", # Use a valid plotmath expression
      sheet_name == "Experiment3" ~ "Salmonella~pre-incubation",
      TRUE ~ sheet_name
    ))
  return(data_long)
}

# 2. Read and process all sheets
file_path <- "All_data_Sal_treated.xlsx"
sheet_names <- c("Experiment1", "Experiment2", "Experiment3")
all_data <- lapply(sheet_names, function(sheet) {
  read_and_process_data(file_path, sheet)
}) %>%
  bind_rows()

# 3. Summarize data
summary_data <- all_data %>%
  group_by(`Time (hrs)`, Species, Experiment) %>%
  summarise(
    mean_percentage = mean(Percentage, na.rm = TRUE),
    se_percentage = sd(Percentage, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 4. Set factor levels with valid plotmath labels
all_data$Experiment <- factor(all_data$Experiment,
  levels = c("Co-inoculation", "BAHP[HP]~pre-incubation", "Salmonella~pre-incubation"),
  labels = c("Co-inoculation", "BAHP[HP]~pre-incubation", "Salmonella~pre-incubation")
)

# 5. Plot
ggplot() +
  geom_point(data = all_data,
             aes(x = `Time (hrs)`, y = Percentage, color = Species),
             alpha = 0.3, position = position_jitter(width = 0.5), size=3) +
  geom_line(data = summary_data,
            aes(x = `Time (hrs)`, y = mean_percentage, color = Species),
            size = 1) +
  geom_ribbon(data = summary_data,
              aes(x = `Time (hrs)`, ymin = mean_percentage - se_percentage,
                  ymax = mean_percentage + se_percentage, fill = Species),
              alpha = 0.2, colour = NA) +
  labs(title = "Percentage of Salmonella and Pseudomonas Over Time (Co-Incubation)",
       x = "Time (hrs)",
       y = "Percentage (%)",
       color = "Species",
       fill = "Species") +
  ylim(0, 100) +
  facet_wrap(~ Experiment, labeller = label_parsed) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "bottom",
    panel.grid = element_blank(),
    strip.text = element_text(size = 20)
  )

# 6. Save plot
ggsave("percentage_overtime.png", width = 20, height = 6, dpi = 300)

# This script was generated with the assistance of Perplexity AI.
# All code has been reviewed and tested for accuracy and functionality by the project team.
