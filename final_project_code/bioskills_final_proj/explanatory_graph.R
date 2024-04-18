#clear enviorment
rm(list=ls())

#load packages
library(stats)
library(ggplot2)
library(dplyr)
library(ggExtra)
library(RColorBrewer)

#install.packages("vegan")
library(vegan)

#load data
raw_data <- read.delim("./merged_all_species.tsv", header = TRUE, row.names = 1, sep = "\t") %>% t()
#load meta-data
meta_data<- read.delim("./my_metadata_species_level.tsv", header = TRUE, sep = "\t", na.strings = c("", "NA", "N/A"), skip = 1)

meta_data <- meta_data %>% distinct() %>% filter(source == "mother_rectal" | source == "child_stool" 
                                  |source == "mother_stool") 



unique(meta_data$source)
row.names(raw_data) <- gsub("[.-]", "_", row.names(raw_data))
meta_data$X <- gsub("[.-]", "_", meta_data$X)
rownames(meta_data) <- meta_data$X

meta_data <- meta_data %>% group_by(family_id) %>% filter(any(grepl("^mother_rectal|^mother_stool", source)) & any(grepl("^child_stool", source)))

match_samples <- function(group){
  matched_samples <- raw_data[row.names(raw_data) %in% group$X]
  return(matched_samples)
}

calculate_avg_distance <- function(group1, group2) {
  dist_matrix <- vegdist(rbind(group1, group2), method = "euclidean")
  avg_distance <- mean(as.matrix(dist_matrix)[lower.tri(dist_matrix)])
  return(avg_distance)
}

vector_gen <- function(cohort){
  uniqe_family <- unique(cohort$family_id)
  family_group <- cohort %>% group_by(family_id)
  
  avg_dist <- vector("numeric", length(uniqe_family))
  for (i in 1:length(uniqe_family)){
    cur_id <- uniqe_family[i] 
    samples <- family_group %>% ungroup %>% filter(family_id == cur_id) %>% group_by(source) 
    row.names(samples) <- samples$X
    unique_sources <- unique(samples$source)
    g1 <- samples %>% filter(source == unique_sources[1])
    
    group1 <- match_samples(samples %>% filter(source == unique_sources[1]))
    group2 <- match_samples(samples %>% filter(source == unique_sources[2]))
    
    avg_dist[i] <- calculate_avg_distance(group1,group2)
    
  }
  return(avg_dist)
}

min_max_normalize <- function(x, factor) {
  x / (factor)
}

unique(meta_data$cohort)
Bedouin <- vector_gen(filter(meta_data, cohort == "mobile"))
Boston <- vector_gen(filter(meta_data, cohort == "origin"))
Finland <- vector_gen(filter(meta_data, cohort == "edia"))

factor <-max(max(Bedouin),max(Boston),max(Finland))
Bedouin <- min_max_normalize(Bedouin, factor)
Boston <- min_max_normalize(Boston, factor)
Finland <- min_max_normalize(Finland, factor)



# creat the graph
data <- data.frame(
  group = c(rep(paste("Finland \n (n = " ,length(Finland),")", sep=""), length(Finland)),
            rep(paste("Bedouin \n(n = " ,length(Bedouin),")",sep=""), length(Bedouin)),
            rep(paste("Boston \n(n  = " ,length(Boston),")",sep = ""), length(Boston))),
  value = c(Finland, Bedouin, Boston)
)
assign_colors <- function(group_name) {
  if (grepl("Finland", group_name)) {
    return('#FF5733')  # Finland color
  } else if (grepl("Bedouin", group_name)) {
    return('#e1a516')  # Bedouin color
  } else if (grepl("Boston", group_name)) {
    return('#13b6b6')  # Boston color
  } else {
    return('#000000')  # Default black color for other cases
  }
}

all_groups <- data$group %>% unique() 
my_comparisons <- list(c(all_groups[1], all_groups[2]), c(all_groups[2], all_groups[3]), c(all_groups[1], all_groups[3]))
graph <- ggplot(data, aes(x = group, y = value, fill = group, color = group)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 16))  +
  guides(fill = "none", color = "none") +
  geom_boxplot(alpha = 0.5) +
  labs(x = "",y = "Euclidean Distance", title = "Maternal-Infant Microbiome Distances")+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = data, aes(x = group, y = value) ,position = position_jitter(width = 0.05)) +
  scale_fill_manual(values = sapply(data$group, assign_colors)) +
  scale_color_manual(values = sapply(data$group, assign_colors)) + 
  stat_compare_means(comparisons = my_comparisons) 


ggsave("explanatory graph.png", plot = graph, width = 8, height =6,dpi=300)




