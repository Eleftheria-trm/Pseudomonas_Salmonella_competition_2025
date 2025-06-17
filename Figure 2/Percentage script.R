# Load required packages
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(emmeans)
library(lme4)

# 1. Define Function to Read and Process Data from Each Sheet
read_and_process_data <- function(file_path, sheet_name) {
  # Read data from the sheet
  data <- read_excel(file_path, sheet = sheet_name)
  
  # Ensure correct data types and rename strains
  data <- data %>%
    mutate(
      Strain = as.factor(Strain),
      `In competition` = as.factor(`In competition`),
      `Cell counted` = as.factor(`Cell counted`),
      `Time (hrs)` = as.numeric(`Time (hrs)`),
      tech.rep = as.factor(tech.rep),
      cfu.ml = as.numeric(cfu.ml)
    )
  
  # Filter for Co-Incubation (STM50 present)
  data_coincubated <- data %>%
    filter(`In competition` == "STM50")
  
  # Calculate Total CFU/ml and Percentages
  data_coincubated <- data_coincubated %>%
    group_by(`Time (hrs)`, tech.rep) %>%
    summarise(
      salmonella_cfu_ml = sum(cfu.ml[`Cell counted` == "Salmonella"], na.rm = TRUE),
      pseudomonas_cfu_ml = sum(cfu.ml[`Cell counted` == "Pseudomonas"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      total_cfu_ml = salmonella_cfu_ml + pseudomonas_cfu_ml,
      percentage_salmonella = (salmonella_cfu_ml / total_cfu_ml) * 100,
      percentage_pseudomonas = (pseudomonas_cfu_ml / total_cfu_ml) * 100
    )
  
  
  # Reshape Data for Combined Plot
  data_long <- data_coincubated %>%
    select(`Time (hrs)`, percentage_salmonella, percentage_pseudomonas) %>%
    pivot_longer(cols = c(percentage_salmonella, percentage_pseudomonas),
                 names_to = "Species",
                 values_to = "Percentage") %>%
    mutate(Species = ifelse(Species == "percentage_salmonella", "Salmonella", "Pseudomonas"))
  
  # Rename Experiment Names
  data_long <- data_long %>%
    mutate(Experiment = case_when(
      sheet_name == "Experiment1" ~ "Co-inoculation",
      sheet_name == "Experiment2" ~ "48hrs pre-incubation with BAHP",
      sheet_name == "Experiment3" ~ "48hrs pre-incubation with Salmonella",
      TRUE ~ sheet_name  # Keep original name if not one of the specified experiments
    ))
  
  return(data_long)
}

# 2. Read and Process Data from All Sheets
  file_path <- "All_data_treated.xlsx"
sheet_names <- c("Experiment1", "Experiment2", "Experiment3")

all_data <- lapply(sheet_names, function(sheet) {
  read_and_process_data(file_path, sheet)
}) %>%
  bind_rows()  # Combine all data frames into one

# 3. Convert Experiment to a factor with specified levels
all_data$Experiment <- factor(all_data$Experiment, levels = c("Co-inoculation", "48hrs pre-incubation with BAHP", "48hrs pre-incubation with Salmonella"))

# 4. Summarize Percentage Data
summary_data <- all_data %>%
  group_by(`Time (hrs)`, Species, Experiment) %>%
  summarise(
    mean_percentage = mean(Percentage, na.rm = TRUE),
    se_percentage = sd(Percentage, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 5. Visualization: Combined Percentages Over Time
ggplot() +
  # Add individual data points
  geom_point(data = all_data,
             aes(x = `Time (hrs)`, y = Percentage, color = Species),
             alpha = 0.3, position = position_jitter(width = 0.5)) +
  # Add mean lines
  geom_line(data = summary_data,
            aes(x = `Time (hrs)`, y = mean_percentage, color = Species),
            size = 1) +
  # Add error shading
  geom_ribbon(data = summary_data,
              aes(x = `Time (hrs)`, ymin = mean_percentage - se_percentage,
                  ymax = mean_percentage + se_percentage, fill = Species),
              alpha = 0.2, colour = NA) +
  # Set labels
  labs(title = "Percentage of Salmonella and Pseudomonas Over Time (Co-Incubation)",
       x = "Time (hrs)",
       y = "Percentage (%)",
       color = "Species",
       fill = "Species") +
  # Set y-axis limits
  ylim(0, 100) +
  # Facet by Experiment
  facet_wrap(~ Experiment) +
  # Customize theme
  theme_minimal() +
  theme(legend.text = element_text(size = 20),  # Adjust the size as needed
        axis.title = element_text(size = 16),   # Size of axis titles
        axis.text = element_text(size = 16),    # Size of axis text (numbers)
        plot.title = element_text(size = 20, hjust = 0.5), # Increase title size and center it
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.text = element_text(size = 20)) # Size of facet titles


# Save the plot longer
ggsave("percentage_salmonella_pseudomonas_combined_all_experiments_longer.png", width = 20, height = 6, dpi = 300)
