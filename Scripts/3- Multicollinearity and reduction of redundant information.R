#### Library ####
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(corrplot)
library(pheatmap)

#### Loading data ####

Ds  <-read.csv("Dataset_2_118_PFQ_MDPs_vs_MPPs.csv")
Dataset$Label <- as.character(Dataset$Label)
Ds.v <- Dataset[,3:120] #only variables

#Normalize variables
norm.fun <- function(x) {(x  - min(x , na.rm=FALSE))/(max (x ,na.rm=TRUE) - min(x , na.rm=FALSE))}
Ds.norm <- as.data.frame(lapply(Ds.v, norm.fun))

#Correlation

Ds_cor    <-cor(Ds.norm, method = c("pearson")) #Pearson
Ds_cor_df <-as.data.frame(Ds_cor)
Ds.cor.mtx<-as.matrix(Ds_cor)

#Cutoff 0.95 

Ds_cor_pfq.95 <-findCorrelation(Ds.cor.mtx,cutoff = .95,verbose = F,names =T)
length(Ds_cor_pfq.95) #29
print(Ds_cor_pfq.95)
PFQ_cor_p.95 <-as.factor(c(Ds_cor_pfq.95))
PFQ_sn_cor.95<-Ds[ , !(names(Ds) %in% PFQ_cor_p.95)]#81 

#Saving data
write.csv(PFQ_sn_cor.95,"~/TRABAJO EN R/Articulo para titulaci?n/Datasets/MDPs_vs_MPPs/Cutoffs/PFQ_corr(0.95)_MDPs_vs_MPPs.csv",row.names = F)

#Cutoff 0.90

Ds_cor_pfq.90 <-findCorrelation(Ds.cor.mtx,cutoff = .90,verbose = F,names =T)
length(Ds_cor_pfq.90) #55
print(Ds_cor_pfq.90)
PFQ_cor_p.90 <-as.factor(c(Ds_cor_pfq.90))
PFQ_sn_cor.90<-Ds[ , !(names(Ds) %in% PFQ_cor_p.90)]#59 

#Saving data
write.csv(PFQ_sn_cor.90,"~/TRABAJO EN R/Articulo para titulaci?n/Datasets/MDPs_vs_MPPs/Cutoffs/PFQ_corr(0.90)_MDPs_vs_MPPs.csv",row.names = F)


#Cutoff 0.85

Ds_cor_pfq.85 <-findCorrelation(Ds.cor.mtx,cutoff = .85,verbose = F,names =T)
length(Ds_cor_pfq.85) #70
print(Ds_cor_pfq.85)
PFQ_cor_p.85 <-as.factor(c(Ds_cor_pfq.85))
PFQ_sn_cor.85<-Ds[ , !(names(Ds) %in% PFQ_cor_p.85)]#48

#Saving data
write.csv(PFQ_sn_cor.85,"~/TRABAJO EN R/Articulo para titulaci?n/Datasets/MDPs_vs_MPPs/Cutoffs/PFQ_corr(0.85)_MDPs_vs_MPPs.csv",row.names = F)

#Cutoff 0.80

Ds_cor_pfq.80 <-findCorrelation(Ds.cor.mtx,cutoff = .80,verbose = F,names =T)
length(Ds_cor_pfq.80) #83
print(Ds_cor_pfq.80)
PFQ_cor_p.80 <-as.factor(c(Ds_cor_pfq.80))
PFQ_sn_cor.80<-Ds[ , !(names(Ds) %in% PFQ_cor_p.80)]#35 

#Saving data
write.csv(PFQ_sn_cor.80,"~/TRABAJO EN R/Articulo para titulaci?n/Datasets/MDPs_vs_MPPs/Cutoffs/PFQ_corr(0.80)_MDPs_vs_MPPs.csv",row.names = F)

#Cutoff 0.75

Ds_cor_pfq.75 <-findCorrelation(Ds.cor.mtx,cutoff = .75,verbose = F,names =T)
length(Ds_cor_pfq.75) #91
print(Ds_cor_pfq.75)
PFQ_cor_p.75 <-as.factor(c(Ds_cor_pfq.75))
PFQ_sn_cor.75<-Ds[ , !(names(Ds) %in% PFQ_cor_p.75)]#27

#Saving data
write.csv(PFQ_sn_cor.75,"~/TRABAJO EN R/Articulo para titulaci?n/Datasets/MDPs_vs_MPPs/Cutoffs/PFQ_corr(0.75)_MDPs_vs_MPPs.csv",row.names = F)

#Cutoff 0.70

Ds_cor_pfq.70 <-findCorrelation(Ds.cor.mtx,cutoff = .70,verbose = F,names =T)
length(Ds_cor_pfq.70)#95
print(Ds_cor_pfq.70)
PFQ_cor_p.70 <-as.factor(c(Ds_cor_pfq.70))
PFQ_sn_cor.70<-Ds[ , !(names(Ds) %in% PFQ_cor_p.70)]#23

#Saving data
write.csv(PFQ_sn_cor.70,"~/TRABAJO EN R/Articulo para titulaci?n/Datasets/MDPs_vs_MPPs/Cutoffs/PFQ_corr(0.70)_MDPs_vs_MPPs.csv",row.names = F)
