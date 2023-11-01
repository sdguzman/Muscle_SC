
#The following code was used to generate the plots in Figure 2


library(ggplot2)
library(ggfortify)
library(dplyr)
library(plotly)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(factoextra)
library(cluster)
library(broom)
library(jtools)
library(umap)
library(forcats)
library(fmsb) # radar plots



setwd("/Users/steed/Dropbox (University of Michigan)/bi_working/guzman")
getwd()


colors= c('#59C7EB', '#2C83C7', '#FDF1AE', '#FEA071', '#18c0c0', '#077187')

#load data
data <- read.csv("nmj_data.csv", header = TRUE)
data$group <- factor(data$group, levels = c("1mo_wt", "1mo_ko", "2mo_wt", "2mo_ko", "3mo_wt", "3mo_ko"))

#Removes any rows that have a missing value. Ensures downstream analysis will continue.
data <- na.omit(data)


# Slice data frame from 3rd to last column, removing labels.
data_num <- data[, 4:ncol(data)]

# Scale the numeric data
data_scaled <- scale(data_num)


# Perform UMAP Analysis
set.seed(42)
umap_result <- umap(data_scaled)

# Convert to a data frame
umap_data <- as.data.frame(umap_result$layout)


# Check if the number of rows are the same
if (nrow(data) == nrow(umap_data)) {
  # Combine the data with the UMAP results
  data <- cbind(data, umap_data)
} else {
  stop("Row dimension in data and umap_data do not match.")
}



umap_plot <- ggplot(data, aes(x = V1 , y = V2, color = group)) +
  geom_point(size = 0.75) +
  scale_color_manual(values = colors) +
  theme_minimal() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  ggtitle("UMAP plot") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        legend.position = 'none')


# Perform k-means clustering
set.seed(123)  # Reset seed for k-means clustering
km_clusters <- kmeans(umap_data[, 1:2], centers = 4)

# Add the cluster assignments to your UMAP dataframe
umap_data$cluster <- as.factor(km_clusters$cluster)

# Compute cluster centers
cluster_centers <- aggregate(cbind(V1, V2) ~ cluster, umap_data, mean)


# Create UMAP plot with k-means cluster and cluster labels

colors = c("#003F5C", "#58508D", "#BC5090", "#FF6361")
# Create mapping of old to new labels
label_mapping <- c("1" = "4", "2" = "1", "3" = "3", "4" = "2")

# Relabel clusters
umap_data <- umap_data %>%
  mutate(cluster = factor(label_mapping[as.character(cluster)]))


pdf("umap_kmeans.pdf", width = 2, height = 2)

umap_kmeans <- ggplot(umap_data, aes(x = V1, y = V2, color = cluster)) +
  geom_point(size = 0.75) +
  theme_minimal() +
  ggtitle("UMAP plot with k-means clusters") +
  labs(x = "UMAP_1", y = "UMAP_2") +
  scale_color_manual(values = colors) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        legend.position = 'none')
umap_kmeans

dev.off()


################################


##Code below is Extended Data Fig 2c


# Extract the variable order from the original dataframe (excluding 'Animal.ID' and 'group' columns)
variable_order <- colnames(data)[!(colnames(data) %in% c("Animal.ID", "group"))]


data_long <- data %>%
  pivot_longer(cols = -c(Animal.ID, group, cluster),
               names_to = "variable",
               values_to = "value")

# Set the factor levels for 'variable' in the long dataframe
data_long$variable <- factor(data_long$variable, levels = variable_order)

# Calculate the mean for each Animal_ID for each variable
animal_means <- data_long %>%
  group_by(Animal.ID, group, variable) %>%
  summarise(animal_mean = mean(value, na.rm = TRUE), .groups = 'drop')

openxlsx::write.xlsx(animal_means, "animal_means.xlsx", rowNames=FALSE)


# Calculate the group mean and SEM using means of each Animal_ID
group_summary <- animal_means %>%
  group_by(group, variable) %>%
  summarise(
    group_mean = mean(animal_mean, na.rm = TRUE),
    sem = sd(animal_mean, na.rm = TRUE) / sqrt(n()), .groups = 'drop'
  )

colors= c('#59C7EB', '#2C83C7', '#FDF1AE', '#FEA071', '#18c0c0', '#077187')

group_levels <- levels(data_long$group)
vline_positions <- seq(2.5, (length(group_levels) - 1.5), by = 2)

#Extended data Figure 2c
pdf("all_means.pdf", width = 11, height = 6)
ggplot(data_long, aes(x = group, y = value)) +
  
  # Individual data points
  geom_jitter(aes(color = group), width = 0.2, size = 1, alpha = 0.6, shape = 16) +
   # Mean of each animal
  geom_point(data = animal_means, aes(y = animal_mean), position = position_jitter(width = 0.4), color = "black", shape = 21, fill = NA, size = 1.5) +
  # Horizontal line for group means
  geom_crossbar(data = group_summary, aes(y = group_mean, ymin = group_mean, ymax = group_mean, fill = group), width = 0.6, position = position_dodge(0.6), size = 0.2) +
  # Group SEM
  geom_errorbar(data = group_summary, aes(y = group_mean, ymin = group_mean - sem, ymax = group_mean + sem), width = 0.2, position = position_dodge(0.6)) +
  # Horizontal line at y=0 for each subplot
  geom_hline(yintercept = 0, linetype="solid", color = "black", size = 0.25) +
  scale_fill_manual(values = colors) +  # Set custom colors
  scale_color_manual(values = colors) + # Also set the jitter points to the same colors
  facet_wrap(~variable, scales = "free") +
  labs(x = "Group", y = "Value") +
  theme_bw(base_size = 16) + 
  geom_vline(xintercept = vline_positions, color = "grey", size = 0.25) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),  
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 10),
    axis.line.y = element_line(color = "black", size = 0.25),  # Adjusted size for y-axis line
    axis.text.y = element_text(color = "black", size = 10),
    axis.ticks.y = element_line(color = "black", size = 0.25), # Adjusted size for y-axis ticks
    panel.spacing.x = unit(1.5, "lines")   # Adjust the space between subplots
  )
dev.off()


# create an empty list to store p-values
t_test_results_list <- list()

# Loop through pairs of groups
for(i in seq(1, nlevels(animal_means$group), by = 2)){
  group1 <- levels(animal_means$group)[i]
  group2 <- levels(animal_means$group)[i + 1]
  
  subset_data <- animal_means %>% filter(group %in% c(group1, group2))
  
  if(length(unique(subset_data$group)) == 2) { # Check if the subset contains both groups
    current_test <- subset_data %>%
      group_by(variable) %>%
      do(tidy(t.test(animal_mean ~ group, data = ., var.equal = TRUE)))
    # Add a column to specify the groups being compared
    current_test$comparison <- paste0(group1, " vs ", group2)
    
    t_test_results_list[[paste0(group1, "_vs_", group2)]] <- current_test
  } else {
    warning(paste("Skipped t-test for", group1, "and", group2, "due to insufficient data."))
  }
}


# Bind all results into a single data frame
t_test_results <- bind_rows(t_test_results_list)

# Show the results
t_test_results

#function  for p-values into asterisks
p_to_asterisk <- function(p){
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  return("NS")
}

#apply function and create a column named significance
t_test_results$significance <- sapply(t_test_results$p.value, p_to_asterisk)

# Save t-test results to Excel
openxlsx::write.xlsx(t_test_results, "t_test_results.xlsx", rowNames=FALSE)


#Counts the numbers of NMJs that are below 10% overlap and then computes a percentage for each animal
data_long_overlap_percentage <- data_long %>%
  filter(variable == "overlap") %>%  # filter only rows where variable is 'overlap'
  group_by(Animal.ID, group) %>%
  summarise(percentage = mean(value <= 10, na.rm = TRUE) * 100,  # compute the percentage
            .groups = "drop")

# Show the result
data_long_overlap_percentage
# Save df_overlap_percentage to an Excel file
openxlsx::write.xlsx(data_long_overlap_percentage, "overlap10_percentage.xlsx", rowNames=FALSE)



#Radar plots. Figure 2h

data <- read.csv("nmj_data.csv", header = TRUE)
data$group <- factor(data$group, levels = c("1mo_wt", "1mo_ko", "2mo_wt", "2mo_ko", "3mo_wt", "3mo_ko"))
data <- na.omit(data)
data <- data %>%
  bind_cols(cluster = umap_data[,"cluster"])

data <- data %>% select(-ID)

# Set the factor levels for 'variable' in the long dataframe
data_long$variable <- factor(data_long$variable, levels = variable_order)



data_long <- data %>%
  pivot_longer(cols = -c(Animal.ID, group,cluster),
               names_to = "variable",
               values_to = "value")
# Set the factor levels for 'variable' in the long dataframe
data_long$variable <- factor(data_long$variable, levels = variable_order)

# 1. Data Transformation
data_wide <- data_long %>%
  group_by(cluster, variable) %>%
  summarise(avg_value = mean(value, na.rm = TRUE)) %>%
  spread(key = variable, value = avg_value) %>%
  as.data.frame()

ordered_vars <- colnames(data_wide)[-1]
max_values <- apply(data_wide[,-1], 2, max)
min_values <- apply(data_wide[,-1], 2, min)

# Add max and min rows for plotting (they will be constant across clusters)
data_for_radar <- rbind(data.frame(cluster = "max", t(max_values)), 
                        data.frame(cluster = "min", t(min_values)), 
                        data_wide)


colors = c("#003F5C", "#58508D", "#BC5090", "#FF6361")

# 2. Radar Chart Plotting

pdf("radar.pdf", width = 4, height = 4)
for(cluster_name in sort(unique(data_long$cluster))) {
  
  # Subset data for the cluster and add max/min rows
  subset_data <- data_for_radar %>%
    filter(cluster %in% c(cluster_name, "max", "min"))
  
  radarchart(subset_data[-1],
             axistype = 1,
             pcol = colors[which(cluster_name == sort(unique(data_long$cluster)))],
             pfcol = adjustcolor(colors[which(cluster_name == sort(unique(data_long$cluster)))], alpha.f = 0.2),
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             cglwd = 0.8,
             vlcex = 0.8,
             title = paste("Radar Chart for", cluster_name))
}
dev.off()

