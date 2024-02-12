#### Packages ####

install.packages("Peptides", dependencies=TRUE)
install_github("dosorio/Peptides")

library(devtools)
library(Peptides)
library(dplyr)
library(usethis)

#### Loading file ####

Training.ds <-read.csv("Trainig_dataset.csv")

#### PCPs with R package ####

#Isoelectric point
Training.ds$pI<-pI(Training.ds$Sequence)

#Hydrophobic scales
Training.ds$hydrophobicity_Eisenberg     <- hydrophobicity(Training.ds$Sequence, scale = "Eisenberg")
Training.ds$hydrophobicity_KyteDoolittle <- hydrophobicity(Training.ds$Sequence, scale = "KyteDoolittle")
Training.ds$hydrophobicity_Cowan7.5      <- hydrophobicity(Training.ds$Sequence, scale = "Cowan7.5")
Training.ds$hydrophobicity_Cowan3.4      <- hydrophobicity(Training.ds$Sequence, scale = "Cowan3.4")
Training.ds$hydrophobicity_Argos         <- hydrophobicity(Training.ds$Sequence, scale = "Argos")
Training.ds$hydrophobicity_Aboderin      <- hydrophobicity(Training.ds$Sequence, scale = "Aboderin")
Training.ds$hydrophobicity_AbrahamLeo    <- hydrophobicity(Training.ds$Sequence, scale = "AbrahamLeo")
Training.ds$hydrophobicity_BlackMould    <- hydrophobicity(Training.ds$Sequence, scale = "BlackMould")
Training.ds$hydrophobicity_BullBreese    <- hydrophobicity(Training.ds$Sequence, scale = "BullBreese")
Training.ds$hydrophobicity_Casari        <- hydrophobicity(Training.ds$Sequence, scale = "Casari")
Training.ds$hydrophobicity_Chothia       <- hydrophobicity(Training.ds$Sequence, scale = "Chothia")
Training.ds$hydrophobicity_Cid           <- hydrophobicity(Training.ds$Sequence, scale = "Cid")
Training.ds$hydrophobicity_Engelman      <- hydrophobicity(Training.ds$Sequence, scale = "Engelman")
Training.ds$hydrophobicity_Fasman        <- hydrophobicity(Training.ds$Sequence, scale = "Fasman")
Training.ds$hydrophobicity_Fauchere      <- hydrophobicity(Training.ds$Sequence, scale = "Fauchere")
Training.ds$hydrophobicity_Goldsack      <- hydrophobicity(Training.ds$Sequence, scale = "Goldsack")
Training.ds$hydrophobicity_Guy           <- hydrophobicity(Training.ds$Sequence, scale = "Guy")
Training.ds$hydrophobicity_HoppWoods     <- hydrophobicity(Training.ds$Sequence, scale = "HoppWoods")
Training.ds$hydrophobicity_Janin         <- hydrophobicity(Training.ds$Sequence, scale = "Janin")
Training.ds$hydrophobicity_Jones         <- hydrophobicity(Training.ds$Sequence, scale = "Jones")
Training.ds$hydrophobicity_Juretic       <- hydrophobicity(Training.ds$Sequence, scale = "Juretic")
Training.ds$hydrophobicity_Kidera        <- hydrophobicity(Training.ds$Sequence, scale = "Kidera")
Training.ds$hydrophobicity_Kuhn          <- hydrophobicity(Training.ds$Sequence, scale = "Kuhn")
Training.ds$hydrophobicity_Levitt        <- hydrophobicity(Training.ds$Sequence, scale = "Levitt")
Training.ds$hydrophobicity_Manavalan     <- hydrophobicity(Training.ds$Sequence, scale = "Manavalan")
Training.ds$hydrophobicity_Miyazawa      <- hydrophobicity(Training.ds$Sequence, scale = "Miyazawa")
Training.ds$hydrophobicity_Parker        <- hydrophobicity(Training.ds$Sequence, scale = "Parker")
Training.ds$hydrophobicity_Ponnuswamy    <- hydrophobicity(Training.ds$Sequence, scale = "Ponnuswamy")
Training.ds$hydrophobicity_Prabhakaran   <- hydrophobicity(Training.ds$Sequence, scale = "Prabhakaran")
Training.ds$hydrophobicity_Rao           <- hydrophobicity(Training.ds$Sequence, scale = "Rao")
Training.ds$hydrophobicity_Rose          <- hydrophobicity(Training.ds$Sequence, scale = "Rose")
Training.ds$hydrophobicity_Roseman       <- hydrophobicity(Training.ds$Sequence, scale = "Roseman")
Training.ds$hydrophobicity_Sweet         <- hydrophobicity(Training.ds$Sequence, scale = "Sweet")
Training.ds$hydrophobicity_Tanford       <- hydrophobicity(Training.ds$Sequence, scale = "Tanford")
Training.ds$hydrophobicity_Welling       <- hydrophobicity(Training.ds$Sequence, scale = "Welling")
Training.ds$hydrophobicity_Wilson        <- hydrophobicity(Training.ds$Sequence, scale = "Wilson")
Training.ds$hydrophobicity_Wolfenden     <- hydrophobicity(Training.ds$Sequence, scale = "Wolfenden")
Training.ds$hydrophobicity_Zimmerman     <- hydrophobicity(Training.ds$Sequence, scale = "Zimmerman")

#General PCPs
Training.ds$hmoment<-hmoment(Training.ds$Sequence)
Training.ds$charge<-charge(Training.ds$Sequence)
Training.ds$mw<-mw(Training.ds$Sequence)
Training.ds$lengthpep<-lengthpep(Training.ds$Sequence)
Training.ds$boman<-boman(Training.ds$Sequence)
Training.ds$instaIndex<-instaIndex(Training.ds$Sequence)

# aaComponents
aaComp  <- aaComp(Training.ds$Sequence)
aaComp  <- as.data.frame(aaComp)

aaComp1 <- aaComp[, c(seq(1, 524, by = 2))]
aaComp1 <- as.data.frame(t(aaComp1), row.names = rownames(Training.ds))


aaComp2 <- aaComp[, c(seq(2, 524, by = 2))]
aaComp2 <- as.data.frame(t(aaComp2), row.names = rownames(Training.ds))

#Rename columns
colnames(aaComp2)[1] <- "%Tiny"
colnames(aaComp2)[2] <- "%Small"
colnames(aaComp2)[3] <- "%Aliphatic"
colnames(aaComp2)[4] <- "%Aromatic"
colnames(aaComp2)[5] <- "%NonPolar"
colnames(aaComp2)[6] <- "%Polar"
colnames(aaComp2)[7] <- "%Charged"
colnames(aaComp2)[8] <- "%Basic"
colnames(aaComp2)[9] <- "%Acidic"

aaComp <- cbind(aaComp1,aaComp2)

Training.ds <- cbind(Training.ds, aaComp)

# zScales

zScales <- zScales(Training.ds$Sequence)


zScales <- as.data.frame(zScales, col.names = rownames(Training.ds))
zScales <- as.data.frame(t(zScales))

colnames(zScales)[1] <- "Lipophilicity"
colnames(zScales)[2] <- "Steric_Properties"
colnames(zScales)[3] <- "Electronic_Properties"
colnames(zScales)[4] <- "Electrophilicity"
colnames(zScales)[5] <- "Hardness"

Training.ds   <- cbind(Training.ds, zScales)

# Cruciani properties
CP <- crucianiProperties(Training.ds$Sequence)
CP <- as.data.frame(CP, col.names = rownames(Training.ds))
CP <- as.data.frame(t(CP))

colnames(CP)[1] <- "Polarity"
colnames(CP)[2] <- "Hydrophobicity_CP"
colnames(CP)[3] <- "H_bonding"

Training.ds<- cbind(Training.ds, CP)

# FASGAI vectors 
FV <- fasgaiVectors(Training.ds$Sequence)
FV <- as.data.frame(FV, col.names = rownames(Training.ds))
FV <- as.data.frame(t(FV))

colnames(FV)[1] <- "Hydrophobicity_index"
colnames(FV)[2] <- "Alpha_and_turn_propensities"
colnames(FV)[3] <- "Bulky_properties,"
colnames(FV)[4] <- "Compositional_characteristic_index"
colnames(FV)[5] <- "Local_flexibility"
colnames(FV)[6] <- "Electronic_properties"

Training.ds<- cbind(Training.ds, FV)

Training.ds$Alpha_and_turn_propensities<-NULL

#Guardar archivos csv

write.csv(Training.ds,"~/TRABAJO EN R/Artículo ordenado/PCPs_training_dataset.csv",row.names = F)
