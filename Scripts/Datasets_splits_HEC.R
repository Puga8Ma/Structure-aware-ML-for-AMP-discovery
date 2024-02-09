# Split the training and external validation sets into 7 structural regions according to their 3-state structural values %H, %E and %C
## Working directory
setwd("/Users/fabienplisson/Desktop/AF2_predictions/")

## Packages
library(ggplot2)
library(ggtern)

## Load datasets
# training set
## Load datasets
mdps_db <- read.csv('./STRIDE_results/MDPs_stride_results.csv', header=TRUE)
row.names(mdps_db) <- mdps_db$PDB_Name
mpps_db <- read.csv('./STRIDE_results/MPPs_stride_results.csv', header=TRUE)
row.names(mpps_db) <- mpps_db$PDB_Name
paps_db <- read.csv('./STRIDE_results/PAPs_stride_results.csv', header=TRUE)
row.names(paps_db) <- paps_db$PDB_Name

# external validation set
ev_db <- read.csv('./STRIDE_results/EV_stride_results.csv', header=TRUE)
row.names(ev_db) <- ev_db$PDB_Name

## Transform datasets
subset_dataset <- function(dataset){
  #Subsetting
  subset <- dataset[,3:5]
  colnames(subset) <- c('helix_H', 'sheet_E', 'coil_C')
  
  # Change null values from 0.0 to 2 (so they can be counted in the density plot)
  subset[subset == 0.0] <- 2
  return(subset)
}

# training set
mdps_sub <- subset_dataset(mdps_db)
mpps_sub <- subset_dataset(mpps_db)
paps_sub <- subset_dataset(paps_db)
training_sub <- rbind(mdps_sub, mpps_sub, paps_sub)
training_sub$Grp <- c(rep("MDPs", 412), rep("MPPs", 327), rep("PAPs", 307))

# external validation set
ev_sub <- subset_dataset(ev_db)
ev_sub$Grp <- c(rep("MDPs", 72), rep("MPPs", 57), rep("PAPs", 133))  

## Filter datasets according to structural regions
dataset = training_sub
#dataset = ev_sub

condition_1 <- (dataset$sheet_E <= 10) & (dataset$helix_H > 10) & (dataset$helix_H <= 90)
region1 <- subset(dataset, condition_1)
  
condition_2 <- (dataset$sheet_E <= 10) & (dataset$helix_H > 90)
region2 <- subset(dataset, condition_2)
  
condition_3 <- (dataset$sheet_E > 10) & (dataset$sheet_E <= 90) & (dataset$coil_C <= 10)
region3 <- subset(dataset, condition_3)
  
condition_4 <- (dataset$sheet_E > 90) & (dataset$coil_C <= 10)
region4 <- subset(dataset, condition_4)
  
condition_5 <- (dataset$helix_H <= 10) & (dataset$coil_C > 10) & (dataset$coil_C <= 90) 
region5 <- subset(dataset, condition_5)
  
condition_6 <- (dataset$helix_H <= 10) & (dataset$coil_C > 90)
region6 <- subset(dataset, condition_6)
  
condition_7 <- (dataset$helix_H > 10) & (dataset$sheet_E > 10) & (dataset$coil_C > 10)
region7 <- subset(dataset, condition_7)
  
# Looking at the distributions in each "region", the number of peptides might be too small for training model. So, we gathered the regions as follows:
# cluster I: regions 1+2 (mainly helices), cluster II: regions 3+4 (mainly beta-sheets), cluster III: regions 5+6 (mainly coils), cluster IV: region 7 (mixed structures)

tr_cluster_I <- rbind(region1, region2)
tr_cluster_II <- rbind(region3, region4)
tr_cluster_III <- rbind(region5, region6)
tr_cluster_IV <- rbind(region7)

write.csv(training_sub, './Structural_subsets/model0_all.csv')
write.csv(tr_cluster_I, './Structural_subsets/model1_helices.csv')
write.csv(tr_cluster_III, './Structural_subsets/model3_coils.csv')
write.csv(tr_cluster_IV, './Structural_subsets/model4_mixes.csv')

ev_cluster_I <- rbind(region1, region2)
ev_cluster_II <- rbind(region3, region4)
ev_cluster_III <- rbind(region5, region6)
ev_cluster_IV <- rbind(region7)

write.csv(ev_sub, './Structural_subsets/validation0_all.csv')
write.csv(ev_cluster_I, './Structural_subsets/validation1_helices.csv')
write.csv(ev_cluster_III, './Structural_subsets/validation3_coils.csv')
write.csv(ev_cluster_IV, './Structural_subsets/validation4_mixes.csv')


# Plot regions

RegDist <- ggtern(NULL, aes(x = sheet_E, y = helix_H,z = coil_C, color='black')) +  
  geom_point(data = training_sub, size=2, shape=21, fill='#8073ac', alpha=0.7) +
  geom_point(data = ev_sub, size=2, shape=21, fill='#e08214', alpha=0.7) +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  scale_color_manual('Datasets', limits=c('model', 'validation'), values = c('#8073ac','#e08214')) +
  guides(colour = guide_legend(override.aes = list(pch = c(21, 21), fill = c('#8073ac','#e08214'))))

RegDist
ggsave('./Ternary_plots/overlay_datasets.pdf')



# Identify clusters by colored zones;
#clusterI.zone <- matrix(c(100, 0, 0, 90, 10, 0, 10, 10, 80, 10, 0, 90),nrow=4,byrow=TRUE)
#colnames(clusterI.zone) <- c('helix_H', 'sheet_E', 'coil_C')

RegDist2 <- ggtern(NULL, aes(x = sheet_E, y = helix_H,z = coil_C)) +  
  geom_point(data = ev_sub, size=1.4, colour='darkgrey', alpha=0.7) +
  #geom_point(data = ev_cluster_I, size=1.4, colour='#7a0177', alpha=0.7) +
  #geom_point(data = ev_cluster_II, size=1.4, colour='#c51b8a', alpha=0.7) +
  #geom_point(data = ev_cluster_III, size=1.4, colour='#f768a1', alpha=0.7) +
  #geom_point(data = ev_cluster_IV, size=1.4, colour='#fbb4b9', alpha=0.7) +
  #geom_polygon(mapping = aes(fill=clusterI.zone),  alpha = 0.75, size = 0.5, color = "black") +
  theme_bw() +
  theme_showarrows() +
  theme_anticlockwise() +
  #guides(ill = guide_legend(title = "Clusters", title.position = "left")) +
  ggtitle('Clustering external validation set into 4 regions')
RegDist2
