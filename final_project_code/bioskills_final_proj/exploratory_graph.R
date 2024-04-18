#clear enviorment
rm(list=ls())

#load packages
library(stats)
library(ggplot2)
library(dplyr)
library(ggExtra)
library(RColorBrewer)

#load data
raw_data <- read.delim("./merged_all_species.tsv", header = TRUE, row.names = 1, sep = "\t")
#load meta-data
meta_data<- read.delim("./my_metadata_species_level.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA", "N/A"), skip = 1)

#transpose data the table- each line become single sample, each column is a bacteria
raw_data <- t(raw_data)

#verify each sample sum to 100
is_approx_100 <- function(x, tolerance = 0.1) {
  abs(x - 100) <= tolerance
}
#sum each row
row_sum <- rowSums(raw_data)
row_sum <- as.data.frame(row_sum)
#apply the function
all_approx_100 <- sapply(row_sum, is_approx_100)
#check if the whole vector is true - print warning message if not
if (all(all_approx_100)) {
  print("Each sample sum to approximately 100")
} else {
  cat("\033[1;31mWarning:Not all of the samples sum to approximately 100[0m\n")
}

#meta-data preprocess
meta_data <- meta_data[!duplicated(meta_data$X), ]
meta_data <- meta_data[complete.cases(meta_data$X), ]
rownames(meta_data) <- meta_data$X

#filter the samples' sources
meta_data <- meta_data %>% filter(source == "mother_rectal" | source == "child_stool" 
                                  |source == "mother_stool")
#replace all NA at time point to 0
meta_data <- meta_data %>% mutate(time_point = ifelse(is.na(time_point), 0, time_point))

#take only the first sample (the minimum of the time_point)
filtered_table <- meta_data %>%
  group_by(source, family_id) %>%
  filter(time_point == min(time_point))

meta_data <- ungroup(filtered_table)

sorted_data <- meta_data %>%
  arrange(family_id, source, desc(case_when(grepl("^mother_rectal", source,) ~ 1, TRUE ~ 0)))

#filter if there are two samples from the same mother (rectal and stool)
filtered_data <- sorted_data %>%
  distinct(family_id, source, .keep_all = TRUE)

# verify there are not duplicates
group_counts <- filtered_data %>%
  group_by(family_id) %>%
  summarize(count = n())

# Filter out groups with a count greater than 2
groups_larger_than_2 <- group_counts %>%
  filter(count > 2)

# Print the duplicates, if there are
print(groups_larger_than_2)

meta_data <- filtered_data

#PCA manipulation
pca<- prcomp(raw_data, center =TRUE)
pca_data<- pca$x

#set empty  label to each sample
raw_data <- as.data.frame(raw_data)
meta_data$X <- gsub("\\.", "_", meta_data$X)
meta_data$X <- gsub("-", "_", meta_data$X)


current_row_names <- rownames(raw_data)

# Replace dots with dashes in the row names
new_row_names <- gsub("\\.", "_", current_row_names)
new_row_names <- gsub("-", "_", new_row_names)

# Set the new row names to the data frame
rownames(raw_data) <- new_row_names

names(meta_data)[names(meta_data)=="X"] <- "RowName"
data <- raw_data %>% filter(rownames(raw_data) %in% meta_data$RowName)

#PCA manipulation
pca<- prcomp(data, center =TRUE)
pca_data<- pca$x

# Create a data frame with the first two principal components and row names as labels
two_dim_pca <- data.frame(PC1 = pca_data[, 1], PC2 = pca_data[, 2], RowName = rownames(data))

#set the row name as the sample id
names(two_dim_pca)[names(two_dim_pca)=="Label"] <- "RowName"
names(meta_data)[names(meta_data)=="X"] <- "RowName"
selected_meta_data <- select(meta_data, RowName, family_id, source, cohort, ShannonIndex)

#joins the useful data
final_data<- inner_join(two_dim_pca, selected_meta_data, by = "RowName")
rownames(final_data) <- final_data$RowName

final_data$source <- gsub("\\bchild_stool\\b", "child", final_data$source, ignore.case = TRUE)
final_data$source <- gsub("\\bmother_rectal\\b", "mother", final_data$source, ignore.case = TRUE)
final_data$source <- gsub("\\bmother_stool\\b", "mother", final_data$source, ignore.case = TRUE)

final_data$cohort <- ifelse(final_data$cohort == "edia", "Finland",
                                 ifelse(final_data$cohort == "origin", "Boston", "Beduian"))
# Plot the 2D PCA plot
pca_biplot <- ggplot(data = final_data, aes(x = PC1, y = PC2, label = RowName, color = cohort, shape = source)) +
  geom_point(alpha = 0.7, size = 4) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_blank()) +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA over mothers and infants samples") +
  scale_colour_manual(values = c("Finland" = '#ea661d', "Beduian" = '#e1a516', "Boston" = '#13b6b6'))
pca_biplot + 
  xlim(range(final_data$PC1) + c(-1, 1) * 2) +
  ylim(range(final_data$PC2) + c(-1, 1) * 2)

#saving the plot to current local directory
ggsave("pca_biplot.png", plot = pca_biplot, width = 8, height =6,dpi=300)

