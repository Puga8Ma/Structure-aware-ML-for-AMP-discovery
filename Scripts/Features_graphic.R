#library
library(dplyr)
library(ggplot2)
library(ggpubr)


#loading data
FI <- read.csv("features.csv")
FI.v<-data.frame(FI[,-1],row.names = FI[,1])
colnames(FI.v)[1] <- "Model 1.0"
colnames(FI.v)[2] <- "Model 1.I"
colnames(FI.v)[3] <- "Model 1.III"
colnames(FI.v)[4] <- "Model 1.IV"
colnames(FI.v)[5] <- "Model 1.V"
columnas_a_redondear <- c("Model 1.0", "Model 1.I","Model 1.III","Model 1.IV","Model 1.V")
FI.v[columnas_a_redondear] <- lapply(FI.v[columnas_a_redondear], function(x) round(x, digits = 3))


# plot data
pdf("features_again.pdf")

my_cols_rev <- c("#F0F921FF","#FCA636FF","#E16462FF","#B12A90FF", "#6A00A8FF","#0D0887FF")

ggballoonplot(FI.v, fill = "value")+
  scale_fill_gradientn(colors = my_cols_rev)

dev.off()


#Transformación de los datos 
str(FI)
summary(FI)

FI.t<- log(FI.v+1)

#### Ternary models ####

#loading data
FI.ternary <- read.csv("features_ternary.csv")
FI.t.v<-data.frame(FI.ternary[,-1],row.names = FI.ternary[,1])
colnames(FI.t.v)[1] <- "Model 2.0"
colnames(FI.t.v)[2] <- "Model 2.I"
colnames(FI.t.v)[3] <- "Model 2.III"
colnames(FI.t.v)[4] <- "Model 2.IV"
colnames(FI.t.v)[5] <- "Model 2.V"

# plot data
pdf("features_ternary_again.pdf")
                        
my_cols_rev <- c("#F0F921FF","#FCA636FF","#E16462FF","#B12A90FF", "#6A00A8FF","#0D0887FF")

ggballoonplot(FI.t.v, fill = "value")+
  scale_fill_gradientn(colors = my_cols_rev)

dev.off()

