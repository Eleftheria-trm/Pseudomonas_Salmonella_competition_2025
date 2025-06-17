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
  
  # Rename "None" to "Untreated" in 'In competition'
  data <- data %>%
    mutate(`In competition` = ifelse(`In competition` == "None", "Untreated", as.character(`In competition`)))
  
  # Add Experiment identifier and rename
  data <- data %>%
    mutate(Experiment = case_when(
      sheet_name == "Experiment1" ~ "Co-inoculation",
      sheet_name == "Experiment2" ~ "48hrs pre-incubation with BAHP",
      sheet_name == "Experiment3" ~ "48hrs pre-incubation with Salmonella",
      TRUE ~ sheet_name  # Keep original name if not one of the specified experiments
    ))
  return(data)
}

# 2. Read and Process Data from All Sheets
file_path <- "All_data_treated.xlsx"
sheet_names <- c("Experiment1", "Experiment2", "Experiment3")

all_data <- lapply(sheet_names, function(sheet) {
  read_and_process_data(file_path, sheet)
}) %>%
  bind_rows()  # Combine all data frames into one

# 3. Filter out "Untreated"
all_data <- all_data %>%
  filter(`In competition` == "STM50")

# 4. Convert Experiment to a factor with specified levels
all_data$Experiment <- factor(all_data$Experiment, levels = c("Co-inoculation", "48hrs pre-incubation with BAHP", "48hrs pre-incubation with Salmonella"))

# 5. Summarize CFU/ml Data
summary_data <- all_data %>%
  group_by(`Time (hrs)`, `Cell counted`, Experiment, `In competition`) %>%
  summarise(
    mean_cfu = mean(cfu.ml, na.rm = TRUE),
    se_cfu = sd(cfu.ml, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 6. Visualization: CFU/ml Over Time
ggplot() +
  # Add individual data points
  geom_point(data = all_data,
             aes(x = `Time (hrs)`, y = cfu.ml, color = `Cell counted`),
             alpha = 0.3, position = position_jitter(width = 0.5)) +
  # Add mean lines
  geom_line(data = summary_data,
            aes(x = `Time (hrs)`, y = mean_cfu, color = `Cell counted`),
            size = 1) +
  # Add error shading
  geom_ribbon(data = summary_data,
              aes(x = `Time (hrs)`, ymin = mean_cfu - se_cfu,
                  ymax = mean_cfu + se_cfu, fill = `Cell counted`),
              alpha = 0.2, colour = NA) +
  # Set labels
  labs(title = "CFU/Seedling Over Time",
       x = "Time (hrs)",
       y = "CFU/seedling",
       color = "Species",
       fill = "Species") +
  # Set y-axis scale to log10
  scale_y_log10() +
  # Facet by Experiment
  facet_wrap(~ Experiment,  nrow = 1) +
  # Customize theme
  theme_minimal() +
  theme(legend.text = element_text(size = 20),  # Adjust the size as needed
        axis.title = element_text(size = 16),   # Size of axis titles
        axis.text = element_text(size = 20),    # Size of axis text (numbers)
        plot.title = element_text(size = 20, hjust = 0.5), # Increase title size and center it
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.text = element_text(size = 20)) # Size of facet titles

# Save the plot
#ggsave("cfu_ml_over_time_stm50_experiments_side_by_side.png", width = 14, height = 8, dpi = 300)

# Save the plot with a shorter x-axis
#ggsave("cfu_ml_over_time_stm50_experiments_short_x_axis.png", width = 10, height = 8, dpi = 300)

# Save the plot with a longer x-axis
ggsave("cfu_ml_over_time_stm50_experiments_long_x_axis.png", width = 22, height = 7, dpi = 300)


