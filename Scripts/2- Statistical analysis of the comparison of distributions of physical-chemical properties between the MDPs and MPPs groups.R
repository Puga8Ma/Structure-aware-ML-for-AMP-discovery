#### Library ####
library(dplyr)
library(arm)
library(caret)
library(sjPlot)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ggsignif)
library(nortest) #Lilliefors test
library(schoolmath)
library(grid) 
library(qvalue) 
library(data.table)

# Loading file
PCPs.raw <- read.csv("Training_dataset.csv")

PCPs <- PCPs.raw[, c(2:ncol(PCPs.raw))] #drop sequences
PCPs$Label <- as.character(PCPs.raw$Label) #label as character

# Divide two classes (0 y 1)
pfq.0  <- pfq[pfq$Label == "MDPs",]#Class 0
pfq.1  <- pfq[pfq$Label == "MPPs",]#Class 1

PCPs <- rbind(pfq.0, pfq.1) #bind two classes

#### BOXPLOTS ####
PCPs.all <- PCPs

boxplots.list <- combn(names(PCPs.all)[1:ncol(PCPs.all)], 2, simplify=FALSE)
boxplots.plot <- list()
for (i in 2:length(PCPs.all)-1) {
  p <- ggplot(PCPs.all,
              aes_string(x = boxplots.list[[i]][1], y = boxplots.list[[i]][2])) +
    geom_boxplot(aes(fill = Label)) +
    theme_bw() +
    scale_fill_manual(values = c("#43a2ca", "#de2d26"), 
                      labels = c("MDPs", "MPPs")) + 
    ggtitle(boxplots.list[[i]][2]) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) + 
    labs(y = boxplots.list[[i]][2], fill = "Class") 
  
  boxplots.plot[[i]] <- p
}

# Saving all plots in pdf (4*4) per sheet
boxplots.multi.pages <- ggarrange(plotlist = boxplots.plot, 
                                  nrow = 4, 
                                  ncol = 4, 
                                  common.legend = T, 
                                  
                                  egend = "right")  

ggexport(boxplots.multi.pages, filename = "Boxplots_General_MDPs_vs_MPPs(0.90).pdf", width = 18, height = 20)

# Density plot
pdf("Density_plots-General_MDPs_vs_MPPs(0.90).pdf")
for (i in 2:ncol(PCPs)) {
  p.0 <- ggplot(data = PCPs.0, aes(x = PCPs.0[, i])) + 
    geom_density(colour = "#43a2ca", fill = "#43a2ca") +
    geom_vline(aes(xintercept = mean(PCPs.0[, i])),
               color="black", linetype="dashed", size=1) +
    stat_function(fun = dnorm, colour = "firebrick",
                  args = list(mean = mean(PCPs.0[, i]),
                              sd = sd(PCPs.0[, i]))) +
    ggtitle("MDPs") +
    labs(x = colnames(PCPs.0)[i]) +
    theme_bw()
  
  print(p.0)
  
  p.1 <- ggplot(data = PCPs.1, aes(x = PCPs.1[, i])) + 
    geom_density(colour = "#de2d26", fill = "#de2d26") +
    geom_vline(aes(xintercept = mean(PCPs.1[, i])),
               color="black", linetype="dashed", size=1) +
    stat_function(fun = dnorm, colour = "firebrick",
                  args = list(mean = mean(PCPs.1[, i]),
                              sd = sd(PCPs.1[, i]))) +
    ggtitle("MPPs density plot") +
    labs(x = colnames(PCPs.1)[i]) +
    theme_bw() 
  
  print(p.1)
}
dev.off()

# Q-plot
pdf("QQ-plot-General_MDPss_vs_MPPss(0.90).pdf")
for (i in 2:ncol(PCPs)) {
  p.q<- ggplot(PCPs, aes(sample = PCPs[, i], colour = Label)) +
    stat_qq(aes(color = Label)) + 
    stat_qq_line(aes(color = Label)) +
    ggtitle(colnames(PCPs)[i]) +
    theme_bw() + 
    scale_color_manual(values = c("#43a2ca", "#de2d26")) + 
    labs(y = colnames(PCPs)[i])
  print(p.q)
}
dev.off()

#### Statistics ####

#  The p-value > 0.05 implying that the distribution of the data are not significantly different from 
#    normal distribution. In other words, we can assume the normality.

# Lilliefors y Shapiro test for each group
NT.0 <- apply(pfq.0[,-c(1)], 2, lillie.test) 
NT.1 <- apply(pfq.1[,-c(1)], 2, lillie.test) # > n=50

# new object for the Lilliefors test p-values of the groups separately
NT.pvalue <- data.frame(matrix(nrow = length(NT.1), ncol = 3))
colnames(NT.pvalue) <- c("Variable" , "0", "1")

# variables name (PCPs)
NT.pvalue$Variable <- colnames(PCPs[, 2:ncol(PCPs)])

# "MDPs" values
for (f in 1:length(NT.0)) {
  p.0 <- NT.0[[f]][2]
  NT.pvalue[f, 2] <- p.0
}

# "MPPs" values
for (f in 1:length(NT.1)) {
  p.1 <- NT.1[[f]][2]
  NT.pvalue[f, 3] <- p.1
}

# Checking how many variables have p-value > 0.05 from the two groups?
sum(NT.pvalue$"0" > 0.05) # 6 normally distributed
sum(NT.pvalue$"1" > 0.05) # 4 normally distributed

# Normally distributed in each group
var.normal.0 <- which(NT.pvalue$"0" > 0.05) ; length(var.normal.0) # 6
var.normal.1 <- which(NT.pvalue$"1" > 0.05) ; length(var.normal.1) # 4

length(intersect(var.normal.0,var.normal.1)) #coincide with a normal distribution in both groups?

var.normal.0.1 <- intersect(var.normal.0,var.normal.1); length(var.normal.0.1) #Wich PCPs coincide?

# Saving variable names with a normally distributed
var.name.normal.0.1  <- NT.pvalue[var.normal.0.1, 1]

# Only columns with normally distributed
pep.0.1.normal <- as.data.frame(PCPs[, var.name.normal.0.1]) ; length(pep.0.1.normal) # 0
pep.0.1.normal <- cbind(PCPs[, 1], pep.0.1.normal) # add label
length(pep.0.1.normal) # 2 (1)
colnames(pep.0.1.normal)[1] <- "Label"
colnames(pep.0.1.normal)[2] <- "Z3"

# Only columns with non-normally distributed data

pep.0.1.NO.normal  <-(PCPs[, -c(var.normal.0.1 + 1)])
pep.0.1.NO.normal  <-(PCPs) #when there is no normal distribution?
length(pep.0.1.NO.normal) # 382 

#### Variance ####

# F-test
#Variance comparison for each PCPs

var.test.normal <- var.test(Z3~Label, 
                            data = pep.0.1.normal,
                            alternative = "two.sided")#one variable

#Saving p-values for variance test

VT.normal.pvalue <- data.frame(matrix(nrow = length(var.test.normal), ncol = 2))

VT.normal.pvalue <- data.frame(matrix(nrow = 1, ncol = 2))#PCPs

colnames(VT.normal.pvalue) <- c("Variable" , "p_value")

#variable names (PCPs)

VT.normal.pvalue$Variable <- colnames(pep.0.1.normal)[2:ncol(pep.0.1.normal)] # 2 porque 1 son Label

VT.normal.pvalue$p_value <- var.test.normal$p.value

# >0.05  (equal variance)

sum(VT.normal.pvalue$p_value > 0.05) # 0 no significant difference 

# <0.05 (different variance)
sum(VT.normal.pvalue$p_value <= 0.05) # 2 significant difference

# equal variance
VT.normal.igual  <- which(VT.normal.pvalue$p_value > 0.05) ; length(VT.normal.igual) # 
VT.name.normal.igual  <- VT.normal.pvalue[VT.normal.igual, 1]

# Saving data with normally distributed and equal variance
pep.0.1.normal.VT.igual <- pep.0.1.normal[, VT.name.normal.igual] ; length(pep.0.1.normal.VT.igual)
pep.0.1.normal.VT.igual <- cbind(pep.0.1.normal[, 1], pep.0.1.normal.VT.igual)
length(pep.0.1.normal.VT.igual) # 0 (0)
colnames(pep.0.1.normal.VT.igual)[1] <- "Label"

# different variance
VT.normal.diferente  <- which(VT.normal.pvalue$p_value < 0.05) ; length(VT.normal.diferente)
VT.name.normal.diferente  <- VT.normal.pvalue[VT.normal.diferente, 1]

# Saving data with normally distributed and different variance
pep.0.1.normal.VT.diferente  <- pep.0.1.normal[, -c(VT.normal.igual + 1)] 
pep.0.1.normal.VT.diferente  <- pep.0.1.normal
length(pep.0.1.normal.VT.diferente) # 2 (1)

#---- Non-normally distributed (Fligner-Killeen test) ----#

var.test.NO.normal <- lapply(pep.0.1.NO.normal[, c(2:ncol(pep.0.1.NO.normal))], 
                             function(x) fligner.test(x~pep.0.1.NO.normal$Label))

# Saving p-values from test variance
VT.NO.normal.pvalue <- data.frame(matrix(nrow = length(var.test.NO.normal), ncol = 2))
colnames(VT.NO.normal.pvalue) <- c("Variable" , "p_value")

# PCPs colnames
VT.NO.normal.pvalue$Variable <- colnames(pep.0.1.NO.normal[, 2:ncol(pep.0.1.NO.normal)]) # 2  

for (i in 1:length(var.test.NO.normal)) {
  p.nnormal <- var.test.NO.normal[[i]][3]
  VT.NO.normal.pvalue[i, 2] <- p.nnormal
}

# > 0.05 (equal variance)
sum(VT.NO.normal.pvalue$p_value > 0.05) # 236 

# < 0.05 (different variance)
sum(VT.NO.normal.pvalue$p_value <= 0.05) # 145 
VT.NO.normal.igual  <- which(VT.NO.normal.pvalue$p_value > 0.05) ; length(VT.NO.normal.igual) # 3

# Saving data with non-normally distributed and equal variance
VT.name.NO.normal.igual  <- VT.NO.normal.pvalue[VT.NO.normal.igual, 1]
pep.0.1.NO.normal.VT.igual <- pep.0.1.NO.normal[, VT.name.NO.normal.igual] ; length(pep.0.1.NO.normal.VT.igual)
pep.0.1.NO.normal.VT.igual <- cbind(pep.0.1.NO.normal[, 1], pep.0.1.NO.normal.VT.igual) 
length(pep.0.1.NO.normal.VT.igual) # 237 (236)
colnames(pep.0.1.NO.normal.VT.igual)[1] <- "Label"

# Saving data with non-normally distributed and different variance
pep.0.1.NO.normal.VT.diferente  <- (pep.0.1.NO.normal[, -c(VT.NO.normal.igual + 1)]) # 1 es porque son las columnas que se eliminan en NT.pvalue
length(pep.0.1.NO.normal.VT.diferente) # 146 (145)

#### Welch test ####

#all PCPs

welch.test <- t.test(Z3~Label,#for only one
                     data=pep.0.1.normal.VT.diferente,
                     var.equal = F, 
                     conf.level = 0.99, 
                     paired = FALSE, 
                     alternative = "two.sided")


# Saving p-values for Welch test
welch.test.pvalue <- data.frame(matrix(nrow = 1, ncol = 2))
colnames(welch.test.pvalue) <- c("Variable" , "p_value")
welch.test.pvalue$Variable <- colnames(pep.0.1.normal.VT.diferente[2:ncol(pep.0.1.normal.VT.diferente)]) #colnames
welch.test.pvalue$p_value <- welch.test$p.value # p-value for one

# N for p-value
sum(welch.test.pvalue$p_value > 0.001) # 0 NS
sum(welch.test.pvalue$p_value <= 0.001) # 1 *

welch.test.1  <- welch.test.pvalue

#### Wilcoxon test ####

# all PCPs

wilcox.test <- lapply(pep.0.1.NO.normal.VT.igual[, c(2:ncol(pep.0.1.NO.normal.VT.igual))],# 3 porque 1 y 2 son Label y ID
                      function(x) wilcox.test(x~pep.0.1.NO.normal.VT.igual$Label,
                                              var.equal   = T,
                                              conf.level  = 0.99,
                                              paired      = FALSE, 
                                              alternative = "two.sided"))

#  p-values for Welch test
wilcox.test.pvalue <- data.frame(matrix(nrow = length(pep.0.1.NO.normal.VT.igual)-1, ncol = 2)) 
colnames(wilcox.test.pvalue) <- c("Variable" , "p_value")
wilcox.test.pvalue$Variable <- colnames(pep.0.1.NO.normal.VT.igual[, 2:ncol(pep.0.1.NO.normal.VT.igual)]) #PCPs names

for (i in 1:(length(pep.0.1.NO.normal.VT.igual)-1)) {
  p.w <- wilcox.test[[i]][3]
  wilcox.test.pvalue[i, 2] <- p.w
}

# N for p-value
sum(wilcox.test.pvalue$p_value > 0.01) # 1 NS
sum(wilcox.test.pvalue$p_value <= 0.01) # 9 *

wilcox.test.NS <- wilcox.test.pvalue[which(wilcox.test.pvalue$p_value > 0.01), ]
wilcox.test.1  <- wilcox.test.pvalue[which(wilcox.test.pvalue$p_value <= 0.01), ]

#### Kolmogorov-smirnov test ####

# Divide "pep.0.1.NO.normal.VT.diferente" btw 0 y 1
pep.0.NO.normal.VT.diferente <- pep.0.1.NO.normal.VT.diferente[pep.0.1.NO.normal.VT.diferente$Label == "MDPs",]
pep.1.NO.normal.VT.diferente <- pep.0.1.NO.normal.VT.diferente[pep.0.1.NO.normal.VT.diferente$Label == "MPPs",]

# all PCPs
ks.test.pvalue <- data.frame(matrix(nrow = length(pep.0.1.NO.normal.VT.diferente), ncol = 2))
colnames(ks.test.pvalue) <- c("Variable" , "p_value")
ks.test.pvalue$Variable <- colnames(pep.0.1.NO.normal.VT.diferente[, 1:ncol(pep.0.1.NO.normal.VT.diferente)])#PCPs names

for (i in 2:ncol(pep.0.NO.normal.VT.diferente)) {
  ks <- ks.test(x = pep.0.NO.normal.VT.diferente[, i], y = pep.1.NO.normal.VT.diferente[, i],
                var.equal   = F,
                conf.level  = 0.99,
                paired      = FALSE, 
                alternative = "two.sided")
  ks.test.pvalue[i,2] <- ks$p.value
}

ks.test.pvalue <- ks.test.pvalue[-c(1),] 

# N for p-value
sum(ks.test.pvalue$p_value > 0.01) # 4 NS
sum(ks.test.pvalue$p_value <= 0.01) #38 *

ks.test.NS <- ks.test.pvalue[which(ks.test.pvalue$p_value > 0.01), ]
ks.test.1  <- ks.test.pvalue[which(ks.test.pvalue$p_value <= 0.01), ]

#### saving PCPs statistically significant ####

welch.test.1$Prueba  <- "welch.test"
wilcox.test.1$Prueba <- "wilcox.test"
ks.test.1$Prueba     <- "ks.test"

PFQ.sig.dif <- rbind(welch.test.1,wilcox.test.1, ks.test.1)
PFQ.sig.dif <- rbind(wilcox.test.1, ks.test.1)

#Save data

write.csv(PFQ.sig.dif,"~/TRABAJO EN R/Artículo ordenado/3. Estadística/Modelo 2/PFQ.sig.did_AAA_MDPs_vs_MPPs.csv",row.names = F)

####
PFQ.sig.dif # Bringing out the ranges of significantly different PCPs
PCPs.1.PFQ.sig.dif  <- pfq.1[, c("Label", PFQ.sig.dif$Variable)]
PCPs.0.PFQ.sig.dif  <- pfq.0[, c("Label", PFQ.sig.dif$Variable)]

# Re-checking

all.equal(colnames(PCPs.1.PFQ.sig.dif[, 2:ncol(PCPs.1.PFQ.sig.dif)]), PFQ.sig.dif$Variable) # checked
all.equal(colnames(PCPs.0.PFQ.sig.dif [, 2:ncol(PCPs.0.PFQ.sig.dif)]) , PFQ.sig.dif$Variable) # checked

#--- PCPs.1 ---#

#saving data
PCPs.1.PFQ.sig.dif.summary <- as.data.frame(matrix(ncol = 6, nrow = ncol(PCPs.1.PFQ.sig.dif)-1))

colnames(PCPs.1.PFQ.sig.dif.summary) <- c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max")

# statistically significant PCPs summary
for (i in 1:(ncol(PCPs.1.PFQ.sig.dif)-1)) {
  print("Summary")
  s <- summary(PCPs.1.PFQ.sig.dif[, -c(1)][[i]])
  print(s)
  # Saving data
  PCPs.1.PFQ.sig.dif.summary[i,] <- s
}

# add colnames for each PCPs
PCPs.1.PFQ.sig.dif.summary$PFQ <- colnames(PCPs.1.PFQ.sig.dif[, -c(1)]) 
PCPs.1.PFQ.sig.dif.summary <- PCPs.1.PFQ.sig.dif.summary[, c(ncol(PCPs.1.PFQ.sig.dif.summary), 1:(ncol(PCPs.1.PFQ.sig.dif.summary)-1))]

# Saving data

write.csv(PCPs.1.PFQ.sig.dif.summary,"~/TRABAJO EN R/Artículo ordenado/3. Estadística/Modelo 2/PCPs.1.PFQ.sig.dif.summary_AAA_MDPs_vs_MPPs.csv",row.names = F)

#--- PCPs.0 ---#

# Saving data
PCPs.0.PFQ.sig.dif.summary <- as.data.frame(matrix(ncol = 6, nrow = ncol(PCPs.0.PFQ.sig.dif)-1))

colnames(PCPs.0.PFQ.sig.dif.summary) <- c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max")

# statistically significant PCPs summary
for (i in 1:(ncol(PCPs.0.PFQ.sig.dif)-1)) {
  print("Summary")
  s <- summary(PCPs.0.PFQ.sig.dif[, -c(1)][[i]])
  print(s)
  # save data
  PCPs.0.PFQ.sig.dif.summary[i,] <- s
}

# add colnames for each PCPs
PCPs.0.PFQ.sig.dif.summary$PFQ <- colnames(PCPs.0.PFQ.sig.dif[, -c(1)]) 
PCPs.0.PFQ.sig.dif.summary <- PCPs.0.PFQ.sig.dif.summary[, c(ncol(PCPs.0.PFQ.sig.dif.summary), 1:(ncol(PCPs.0.PFQ.sig.dif.summary)-1))]

# save data
write.csv(PCPs.0.PFQ.sig.dif.summary,"~/TRABAJO EN R/Artículo ordenado/3. Estadística/Modelo 2/PCPs.0.PFQ.sig.dif.summary_AAA_MDPs_vs_MPPs.csv",row.names = F)


#### Comparing the averages of each PCPs ####

PCPs.0.1.PFQ.sig.dif.prom <- as.data.frame(matrix(ncol = 2, nrow = nrow(PCPs.0.PFQ.sig.dif.summary)))
colnames(PCPs.0.1.PFQ.sig.dif.prom) <- c("PCPs", "Comparacion promedios")
PCPs.0.1.PFQ.sig.dif.prom$PCPs <- PCPs.0.PFQ.sig.dif.summary$PFQ

for (i in 1:nrow(PCPs.0.PFQ.sig.dif.summary)) {
  print("Promedio MDPs")
  prom.S <- PCPs.0.PFQ.sig.dif.summary[i,5]
  print(prom.S)
  print("Promedio MPPs")
  prom.NS <- PCPs.1.PFQ.sig.dif.summary[i,5]
  print(prom.NS)
  if (prom.S > prom.NS) {
    print("PCPs.0 > PCPs.1")
    PCPs.0.1.PFQ.sig.dif.prom[i, 2] <- "PCPs.0 > PCPs.1"
  } else {
    print("PCPs.0 < PCPs.1")
    PCPs.0.1.PFQ.sig.dif.prom[i, 2] <- "PCPs.0 < PCPs.1"
  }
}

#add averages and quartiles
PCPs.0.1.PFQ.sig.dif.prom$PCPs.0.1Qu  <- PCPs.0.PFQ.sig.dif.summary$`1st Qu`
PCPs.0.1.PFQ.sig.dif.prom$PCPs.0.Mean <- PCPs.0.PFQ.sig.dif.summary$Mean
PCPs.0.1.PFQ.sig.dif.prom$PCPs.0.3Qu  <- PCPs.0.PFQ.sig.dif.summary$`3rd Qu`

PCPs.0.1.PFQ.sig.dif.prom$PCPs.1.1Qu  <- PCPs.1.PFQ.sig.dif.summary$`1st Qu`
PCPs.0.1.PFQ.sig.dif.prom$PCPs.1.Mean <- PCPs.1.PFQ.sig.dif.summary$Mean
PCPs.0.1.PFQ.sig.dif.prom$PCPs.1.3Qu  <- PCPs.1.PFQ.sig.dif.summary$`3rd Qu`

# saving data

write.csv(PCPs.0.1.PFQ.sig.dif.prom,"~/TRABAJO EN R/Artículo ordenado/3. Estadística/Modelo 2/PCPs.0.1.PFQ.sig.dif.prom_AAA_MDPs_vs_MPPs.csv",row.names = F)

#### Boxplots ####

pep.0.1.normal.VT.diferente$Label    <- as.character(pep.0.1.normal.VT.diferente$Label)
pep.0.1.NO.normal.VT.igual$Label     <- as.character(pep.0.1.NO.normal.VT.igual$Label)
pep.0.1.NO.normal.VT.diferente$Label <- as.character(pep.0.1.NO.normal.VT.diferente$Label)

# Welch test
welch.list <- combn(names(pep.0.1.normal.VT.diferente)[1:ncol(pep.0.1.normal.VT.diferente)], 2, simplify=FALSE)
welch.plot <- list()
for (i in 2:length(pep.0.1.normal.VT.diferente)-1) {
  p <- ggplot(pep.0.1.normal.VT.diferente,
              aes_string(x = welch.list[[i]][1], y = welch.list[[i]][2])) +
    geom_boxplot(aes(fill = Label)) +
    theme_bw() +
    scale_fill_manual(values = c("#43a2ca", "#de2d26"), labels = c("MDPs", "MPPs")) +
    ggtitle(welch.list[[i]][2], "Welch Test") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    labs(y = welch.list[[i]][2], fill = "Class") +
    geom_signif(comparisons = list(c("MDPs", "MPPs")),
                map_signif_level = c("*"=0.01),
                test = "t.test")
  welch.plot[[i]] <- p
}

# Wilcoxon
wilcox.list <- combn(names(pep.0.1.NO.normal.VT.igual)[1:ncol(pep.0.1.NO.normal.VT.igual)], 2, simplify=FALSE)
wilcox.plot <- list()
for (i in 2:length(pep.0.1.NO.normal.VT.igual)-1) {
  p <- ggplot(pep.0.1.NO.normal.VT.igual,
              aes_string(x = wilcox.list[[i]][1], y = wilcox.list[[i]][2])) +
    geom_boxplot(aes(fill = Label)) +
    theme_bw() +
    scale_fill_manual(values = c("#43a2ca", "#de2d26"), labels = c("MDPs", "MPPs")) + 
    ggtitle(wilcox.list[[i]][2], "Wilcoxon Test") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
    labs(y = wilcox.list[[i]][2], fill = "Class") + 
    geom_signif(comparisons = list(c("MDPs", "MPPs")), 
                map_signif_level = c("*"=0.01),
                test = "wilcox.test")
  
  wilcox.plot[[i]] <- p
}

# Kolmogorov-Smirnov
KS.list <- combn(names(pep.0.1.NO.normal.VT.diferente)[1:ncol(pep.0.1.NO.normal.VT.diferente)], 2, simplify=FALSE)
KS.plot <- list()
for (i in 2:length(pep.0.1.NO.normal.VT.diferente)-1) {
  p <- ggplot(pep.0.1.NO.normal.VT.diferente,
              aes_string(x = KS.list[[i]][1], y = KS.list[[i]][2])) +
    geom_boxplot(aes(fill = Label)) +
    theme_bw() +
    scale_fill_manual(values = c("#43a2ca", "#de2d26"), labels = c("MDPs", "MPPs")) + 
    ggtitle(KS.list[[i]][2], "Kolmogorov-Smirnov Test") +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    labs(y = KS.list[[i]][2], fill = "Class") + 
    geom_signif(comparisons = list(c("MDPs", "MPPs")), 
                map_signif_level = c("*"=0.01),
                test = "ks.test")
  
  KS.plot[[i]] <- p
}

#### saving all plots in pdf ####

all_plot <- do.call(c, list(welch.plot[1:length(welch.plot)],wilcox.plot[1:length(wilcox.plot)], KS.plot[1:length(KS.plot)]))

all_plot <- do.call(c, list(wilcox.plot[1:length(wilcox.plot)], KS.plot[1:length(KS.plot)]))

multi.pages <- ggarrange(plotlist = all_plot, 
                         nrow = 4, 
                         ncol = 4, 
                         common.legend = T, 
                         legend = "right")

ggexport(multi.pages, filename = "Boxplots_MDPs_vs_MPPs_General_difsig_AAA.pdf", width = 18, height = 20)



welch.pvalue.plot <-
  ggplot(data = welch.test.pvalue, 
         mapping = aes(elch.test.pvalue$p_value)) +
  geom_histogram(fill = "cyan3") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  theme_bw() +
  ggtitle("Welch Test") +
  labs(x = "p-value")

#### FDR (False Discovery Rate) ####

# p.adjust

all.p.adjust <- as.data.frame(matrix(ncol = 8, 
                                     nrow = (nrow(welch.test.pvalue)+nrow(wilcox.test.pvalue) + nrow(ks.test.pvalue))))

colnames(all.p.adjust) <- c("Variable", "Prueba", "p.raw", "p.raw<=0.01", "p.BH", "p.BH<=0.01", "p.BY", "p.BY<=0.01")

all.p.adjust[c(1:(nrow(welch.test.pvalue))), c(1,3)]   <- welch.test.pvalue
all.p.adjust[c((nrow(welch.test.pvalue)+1):(nrow(wilcox.test.pvalue)+nrow(welch.test.pvalue))), c(1,3)] <- wilcox.test.pvalue
all.p.adjust[c((nrow(wilcox.test.pvalue)+nrow(welch.test.pvalue)+1):(nrow(ks.test.pvalue)+nrow(welch.test.pvalue)+nrow(wilcox.test.pvalue))), c(1,3)] <- ks.test.pvalue

rep(x = "Welch",              times = nrow(welch.test.pvalue)),
rep(x = "Wilcoxon",           times = nrow(wilcox.test.pvalue)),
rep(x = "Kolmogorov-Smirnov", times = nrow(ks.test.pvalue)))


all.p.adjust$Prueba <- c(rep(x = "Welch",              times = nrow(welch.test.pvalue)),
                         rep(x = "Wilcoxon",           times = nrow(wilcox.test.pvalue)),
                         rep(x = "Kolmogorov-Smirnov", times = nrow(ks.test.pvalue)))


# BH and BY methods
all.p.adjust$p.BH <- p.adjust(p = all.p.adjust$p.raw, method = "BH")
all.p.adjust$p.BY <- p.adjust(p = all.p.adjust$p.raw, method = "BY")

# Results
all.p.adjust$`p.raw<=0.01` <- all.p.adjust$p.raw <= 0.01
all.p.adjust$`p.BH<=0.01`  <- all.p.adjust$p.BH  <= 0.01
all.p.adjust$`p.BY<=0.01`  <- all.p.adjust$p.BY  <= 0.01

table(all.p.adjust$p.raw <= 0.01) # 83
table(all.p.adjust$p.BH  <= 0.01) # 83
table(all.p.adjust$p.BY  <= 0.01) # 82

#### Saving FDR results ####

write.csv(all.p.adjust,"~/TRABAJO EN R/Artículo ordenado/3. Estadística/Modelo 2/all.p.adjust_AAA_MDPs_vs_MPPs.csv",row.names = F)
