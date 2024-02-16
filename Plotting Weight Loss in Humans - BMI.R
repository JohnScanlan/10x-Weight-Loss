# Load necessary library
library(ggplot2)

# Define the data
bmi_data <- data.frame(
  Participant = factor(1:4), # Assign participant IDs
  Before = c(65, 62, 68, 80),
  After = c(61, 59, 65, 74)
)

# Reshape the data from wide to long format for ggplot2
bmi_long <- tidyr::pivot_longer(bmi_data, cols = c("Before", "After"), names_to = "Time", values_to = "BMI")

# Assuming bmi_long is already created
bmi_long$Time <- factor(bmi_long$Time, levels = c("Before", "After"))

# Create the adjusted line graph with specified customizations
ggplot(bmi_long, aes(x = Time, y = BMI, group = Participant)) +
  geom_line(aes(color = "black"), size = 1.5) + # Set lines to black and make them bolder
  geom_point(aes(color = "black"), size = 4) + # Increase point size and set to black
  scale_color_manual(values = c("black")) + # Ensure all colors are set to black
  theme_minimal(base_size = 14) + # Use a minimal theme
  labs(title = "Change in BMI Over 6 Weeks of Weight Loss",
       x = "Time",
       y = "BMI") + # Customize labels
  theme(legend.title = element_blank(), # Remove the legend title
        legend.position = "none", # Remove legend as color distinction is no longer needed
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        axis.line = element_line(colour = "black")) # Enhance axis lines for clarity
