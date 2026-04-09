#####################################################

# PCA genetic clustering of transcriptome and CAPTUS
# code developed by Abby Pearse, Jessie Pelosi
# Last updated: 04/01/2026

#####################################################

### Libraries

library(vcfR)
library(ggplot2)
library(readr)
library(dplyr)
library(viridis)
library(paletteer)
library(gghighlight)

### Reading in files

transcriptome.vcf <- read.vcfR("files/txm.50missing.CTmarked.recalc.maf0.05.thinned.vcf")
captus.vcf <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")
popmap <- read.csv("files/popmap_updatedcoords.csv")

####################################### Creating genind objects ##########################################

### Transcriptome

transcriptome.genind <- vcfR2genind(transcriptome.vcf, ploidy = 4)
pop(transcriptome.genind) <- popmap$Temporal_Group
transcriptome.genind.nomissing <- missingno(transcriptome.genind, type = "mean")

### CAPTUS

captus.genind <- vcfR2genind(captus.vcf, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group
captus.genind.nomissing <- missingno(captus.genind, type = "mean")


############################################ Generating PCA dataframes ###########################################

### Transcriptome

pca.transcriptome <- dudi.pca(transcriptome.genind.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.transcriptome <- (pca.transcriptome$eig / sum(pca.transcriptome$eig))*100
pve.transcriptome <- round(pve.transcriptome, digits = 2)

pca.transcriptome.df <- pca.transcriptome$li
pca.transcriptome.df$Ind <- rownames(pca.transcriptome.df)
pca.transcriptome.df.popmap <- inner_join(pca.transcriptome.df, popmap)


### CAPTUS

pca.captus <- dudi.pca(captus.genind.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.captus <- (pca.captus$eig / sum(pca.captus$eig))*100
pve.captus <- round(pve.captus, digits = 2)

pca.captus.df <- pca.captus$li
pca.captus.df$Ind <- rownames(pca.captus.df)
pca.captus.df.popmap <- inner_join(pca.captus.df, popmap)


######################################### Plotting PCAs #####################################

native.locs <- c("Japan", "China", "Indonesia", "Palau", "Philippines", "Taiwan", "Vietnam")
invaded.locs <- c("Florida", "Alabama", "Arkansas", "Georgia", "Louisiana", "Mississippi", "South Carolina", "Texas")

### Transcriptome

# Collection year - collection year legend disappeared :((

pca.transcriptome.collectionyear <- ggplot(pca.transcriptome.df.popmap, 
                                           aes(Axis1, Axis2, color = as.numeric(Collection.Year), shape = Invaded)) +
  geom_point(size = 3.5) +
  geom_point(data = subset(pca.transcriptome.df.popmap, Location %in% invaded.locs), aes(color = as.numeric(Collection.Year))) +
  gghighlight(Location %in% invaded.locs, unhighlighted_params = list(color = "black", alpha = 1)) +
  scale_color_viridis_c(limits = range(pca.transcriptome.df.popmap$Collection.Year)) +
  xlab(paste0("PC 1 (", pve.transcriptome[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome[2], "% variation explained)")) +
  theme_bw() +
  labs(color = "Collection Year")

#ggsave("transcriptome.collectionyear.pca.pdf", height = 6, width = 8)

# Location

pca.transcriptome.location <- ggplot(pca.transcriptome.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.transcriptome[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome[2], "% variation explained)")) +
  ggtitle("Transcriptome") +
  theme_bw()

#ggsave("transcriptome.location.pca.pdf", height = 6, width = 8)

### CAPTUS

# Collection year

pca.captus.collectionyear <- ggplot(pca.captus.df.popmap, aes(Axis1, Axis2, color = Collection.Year, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  ggtitle("CAPTUS") +
  theme_bw()

#ggsave("captus.collectionyear.pca.pdf", height = 6, width = 8)

# Location

pca.captus.location <- ggplot(pca.captus.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  ggtitle("CAPTUS") +
  theme_bw()

#ggsave("captus.location.pca.pdf", height = 6, width = 8)



##################### Creating PCAs exlcuding Pacific Islands (for higher resolution of invaded?) ################

pacific.islands <- c("SRR29127777", "SRR29127778", "SRR29127780", "SRR29127781", "SRR29127782", "SRR29127782",
                     "SRR29127775", "SRR29127776") # removes Taiwan, Phillippines, and Palau

### Creating and subsetting new geninds

# Transcriptome

transcriptome.genind.noPI <- transcriptome.genind[!indNames(transcriptome.genind) %in% pacific.islands]
transcriptome.genind.noPI.nomissing <- missingno(transcriptome.genind.noPI, type = "mean")

# CAPTUS

captus.genind.noPI <- captus.genind[!indNames(captus.genind) %in% pacific.islands]
captus.genind.noPI.nomissing <- missingno(captus.genind.noPI, type = "mean")

#################################### Creating PCA dataframes ###########################################

### Transcriptome

pca.transcriptome.noPI <- dudi.pca(transcriptome.genind.noPI.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.transcriptome.noPI <- (pca.transcriptome.noPI$eig / sum(pca.transcriptome.noPI$eig))*100
pve.transcriptome.noPI <- round(pve.transcriptome.noPI, digits = 2)

pca.transcriptome.noPI.df <- pca.transcriptome.noPI$li
pca.transcriptome.noPI.df$Ind <- rownames(pca.transcriptome.noPI.df)
pca.transcriptome.noPI.df.popmap <- inner_join(pca.transcriptome.noPI.df, popmap)

### CAPTUS

pca.captus.noPI <- dudi.pca(captus.genind.noPI.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.captus.noPI <- (pca.captus.noPI$eig / sum(pca.captus.noPI$eig))*100
pve.captus.noPI <- round(pve.captus.noPI, digits = 2)

pca.captus.noPI.df <- pca.captus.noPI$li
pca.captus.noPI.df$Ind <- rownames(pca.captus.noPI.df)
pca.captus.noPI.df.popmap <- inner_join(pca.captus.noPI.df, popmap)


####################################### Plotting PCAs ################################

### Transcriptome

# Collection Year

pca.transcriptome.noPI.collectionYear <- ggplot(pca.transcriptome.noPI.df.popmap, aes(Axis1, Axis2, color = Collection.Year, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  xlab(paste0("PC 1 (", pve.transcriptome.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome.noPI[2], "% variation explained)")) +
  ggtitle("Transcriptome") +
  theme_bw()

#ggsave("transcriptome.noPI.collectionyear.pca.pdf", width = 8, height = 6)

# Location

pca.transcriptome.noPI.location <- ggplot(pca.transcriptome.noPI.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.transcriptome.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome.noPI[2], "% variation explained)")) +
  ggtitle("Transcriptome") +
  theme_bw()

#ggsave("transcriptome.noPI.location.pca.pdf", width = 8, height = 6)

### CAPTUS

# Collection year

pca.captus.noPI.location.pca <- ggplot(pca.captus.noPI.df.popmap, aes(Axis1, Axis2, color = Collection.Year, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_viridis() +
  xlab(paste0("PC 1 (", pve.captus.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus.noPI[2], "% variation explained)")) +
  ggtitle("CAPTUS") +
  theme_bw()

#ggsave("captus.noPI.collectionyear.pca.pdf", width = 8, height = 6)

# Location

pca.captus.noPI.pca <- ggplot(pca.captus.noPI.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 3) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.captus.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus.noPI[2], "% variation explained)")) +
  ggtitle("CAPTUS") +
  theme_bw()

#ggsave("captus.noPI.location.pca.pdf", width = 8, height = 6)




################################### temporal groups - excluded from analyses ########################

### Transcriptome

pca.transcriptome.tempgroups <- ggplot(pca.transcriptome.df.popmap, aes(Axis1, Axis2, color = Temporal_Group)) + geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC 1 (", pve.transcriptome[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome[2], "% variation explained)")) +
  scale_color_manual(values = c ("native" = "dodgerblue2", "early invasion" = "pink",
                                 "mid invasion" = "deeppink2", "late invasion" = "firebrick4"))+
  ggtitle("TRANSCRIPTOME PCA") 

#ggsave("transcriptome.temporalgroups.pdf", width = 8, height = 6)

### CAPTUS

pca.captus.tempgroups <- ggplot(pca.captus.df.popmap, aes(Axis1, Axis2, color = Temporal_Group)) + geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  scale_color_manual(values = c ("native" = "dodgerblue2", "early invasion" = "pink",
                                 "mid invasion" = "deeppink2", "late invasion" = "firebrick4")) +
  ggtitle("CAPTUS PCA") 

#ggsave("captus.temporalgroups.pdf", height = 6, width = 8)
