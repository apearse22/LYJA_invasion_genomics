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
library(adegenet)
library(poppr)

### Reading in files

captus.vcf <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")
popmap <- read.csv("files/popmap_updatedcoords.csv")

####################################### Creating genind objects ##########################################

### CAPTUS

captus.genind <- vcfR2genind(captus.vcf, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group
captus.genind.nomissing <- missingno(captus.genind, type = "mean")


############################################ Generating PCA dataframes ###########################################

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

### CAPTUS

# Collection year
pca.captus.df.inv <- filter(pca.captus.df.popmap, Invaded == "Y")
pca.captus.df.nat <- filter(pca.captus.df.popmap, Invaded == "N")

pca.captus.collectionyear <- ggplot() +
  geom_point(pca.captus.df.inv, mapping = aes(x = Axis1, y = Axis2, color = as.numeric(Collection.Year)), size = 3, alpha = 0.85) +
  geom_point(pca.captus.df.nat, mapping = aes(x = Axis1, y = Axis2), color = "black", shape = 17, size = 2) +
  scale_color_viridis_c(limits = range(pca.captus.df.popmap$Collection.Year)) +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +

  labs(color = "Collection Year") +
  ggtitle("CAPTUS Reference")

ggsave("captus.collectionyear.pca.pdf", height = 6, width = 8)
ggsave("captus.collectionyear.pca.png", height = 6, width = 8, dpi = 300)



# pca colored by native geographic location
pca.captus.collectionyear <- ggplot() +
  geom_point(pca.captus.df.inv, mapping = aes(x = Axis1, y = Axis2, color = "black", size = 3, alpha = 0.85)) +
  geom_point(pca.captus.df.nat, mapping = aes(x = Axis1, y = Axis2, color = Location), shape = 17, size = 2) +
  scale_color_viridis_c(limits = range(pca.captus.df.popmap$Collection.Year)) +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  
  labs(color = "Collection Year") +
  ggtitle("CAPTUS Reference")

ggsave("captus.collectionyear.pca.pdf", height = 6, width = 8)
ggsave("captus.collectionyear.pca.png", height = 6, width = 8, dpi = 300)


























#ggsave("captus.collectionyear.pca.pdf", height = 6, width = 8)

# Location

pca.captus.location <- ggplot(pca.captus.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  ggtitle("CAPTUS Reference") +
  theme_bw()

ggsave("captus.location.pca.pdf", height = 6, width = 8)
ggsave("captus.location.pca.png", height = 6, width = 8, dpi = 300)



##################### Creating PCAs exlcuding Pacific Islands (for higher resolution) ################

pacific.islands <- c("SRR29127777", "SRR29127778", "SRR29127780", "SRR29127781", "SRR29127782", "SRR29127782",
                     "SRR29127775", "SRR29127776") # removes Taiwan, Philippines, and Palau

### Creating and subsetting new geninds

# CAPTUS

captus.genind.noPI <- captus.genind[!indNames(captus.genind) %in% pacific.islands]
captus.genind.noPI.nomissing <- missingno(captus.genind.noPI, type = "mean")

#################################### Creating PCA dataframes ###########################################

### CAPTUS

pca.captus.noPI <- dudi.pca(captus.genind.noPI.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.captus.noPI <- (pca.captus.noPI$eig / sum(pca.captus.noPI$eig))*100
pve.captus.noPI <- round(pve.captus.noPI, digits = 2)

pca.captus.noPI.df <- pca.captus.noPI$li
pca.captus.noPI.df$Ind <- rownames(pca.captus.noPI.df)
pca.captus.noPI.df.popmap <- inner_join(pca.captus.noPI.df, popmap)


####################################### Plotting PCAs ################################

### CAPTUS

# Collection year

#pca.captus.noPI.location.pca <- ggplot(pca.captus.noPI.df.popmap, aes(Axis1, Axis2, color = Collection.Year, shape = Invaded)) +
#  geom_point(size = 3) +
#  scale_color_viridis() +
#  xlab(paste0("PC 1 (", pve.captus.noPI[1], "% variation explained)")) +
#  ylab(paste0("PC 2 (", pve.captus.noPI[2], "% variation explained)")) +
#  ggtitle("CAPTUS") +
#  theme_bw()

#ggsave("captus.noPI.collectionyear.pca.pdf", width = 8, height = 6)

# Location


pca.captus.noPI.inv.df <- pca.captus.noPI.df.popmap %>% 
  filter(Invaded == "Y")

pca.captus.noPI.nat.df <- pca.captus.noPI.df.popmap %>% 
  filter(Invaded == "N")

# invaded total
pca.captus.collectionyear.pca <- ggplot() +
  geom_point(pca.captus.df.inv, mapping = aes(x = Axis1, y = Axis2, color = as.numeric(Collection.Year)), size = 3, alpha = 0.85) +
  geom_point(pca.captus.df.nat, mapping = aes(x = Axis1, y = Axis2), color = "black", shape = 17, size = 3) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_viridis_c(limits = range(pca.captus.noPI.df.popmap$Collection.Year)) +
  xlab(paste0("PC 1 (", pve.captus.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus.noPI[2], "% variation explained)")) +
  theme_bw() +
  labs(color = "Collection Year") +
  theme(text = element_text(size = 16)) 

ggsave("captus.collectionyear.inv.pca.pdf", width = 8, height = 6)
ggsave("captus.collectionyear.inv.pca.png", width = 8, height = 6, dpi = 300)



# no pacific islands by invaded collection year
pca.captus.noPI.collectionyear.pca <- ggplot() +
  geom_point(pca.captus.noPI.inv.df, mapping = aes(x = Axis1, y = Axis2, color = as.numeric(Collection.Year)), size = 3, alpha = 0.85) +
  geom_point(pca.captus.noPI.nat.df, mapping = aes(x = Axis1, y = Axis2), color = "black", shape = 17, size = 3) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_viridis_c(limits = range(pca.captus.noPI.df.popmap$Collection.Year)) +
  xlab(paste0("PC 1 (", pve.captus.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus.noPI[2], "% variation explained)")) +
  theme_bw() +
  labs(color = "Collection Year") +
  theme(text = element_text(size = 16)) 
  
ggsave("captus.collectionyear.noPI.pca.pdf", width = 8, height = 6)
ggsave("captus.collectionyear.noPI.pca.png", width = 8, height = 6, dpi = 300)


# native location
pca.captus.collectionyear <- ggplot() +
  geom_point(pca.captus.df.inv, mapping = aes(x = Axis1, y = Axis2), color = "black", size = 3, alpha = 0.85) +
  geom_point(pca.captus.df.nat, mapping = aes(x = Axis1, y = Axis2, color = Location), shape = 17, size = 3) +
  #scale_color_viridis_c(limits = range(pca.captus.df.popmap$Collection.Year)) +
  scale_color_manual(values = c("China" ="chartreuse3",
                                "Japan" = "darkorchid3",
                                "Indonesia" = "firebrick3",
                                "Philippines" = "dodgerblue2",
                                "Vietnam" = "gold1",
                                "Taiwan" = "lightblue",
                                "Palau" = "pink")) +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  theme_bw() +
  labs(color = "Native Location") +
  theme(text = element_text(size = 16)) 
  

ggsave("captus.nativelocation.pca.pdf", width = 8, height = 6)
ggsave("captus.nativelocation.pca.png", width = 8, height = 6, dpi = 300)


# no Pacific islands native location
pca.captus.noPI.collectionyearnative.pca <- ggplot() +
  geom_point(pca.captus.noPI.nat.df, mapping = aes(x = Axis1, y = Axis2, color = Location), size = 3, alpha = 0.85, shape = 17) +
  geom_point(pca.captus.noPI.inv.df, mapping = aes(x = Axis1, y = Axis2), color = "black", size = 3) +
  scale_shape_manual(values = c(17,16)) +
  #scale_color_viridis_c(limits = range(pca.captus.noPI.df.popmap$Collection.Year)) +
  scale_color_manual(values = c("China" ="chartreuse3",
                                "Japan" = "darkorchid3",
                                "Indonesia" = "firebrick3",
                                "Philippines" = "dodgerblue2",
                                "Vietnam" = "gold1")) +
  xlab(paste0("PC 1 (", pve.captus.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus.noPI[2], "% variation explained)")) +
  theme_bw() +
  labs(color = "Native Location") +
  theme(text = element_text(size = 16)) 

ggsave("captus.nativelocation.noPI.pca.pdf", width = 8, height = 6)
ggsave("captus.nativelocation.noPI.pca.png", width = 8, height = 6, dpi = 300)


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




################################################# Transcriptome analyses ###############################################

transcriptome.vcf <- read.vcfR("files/txm.50missing.CTmarked.recalc.maf0.05.thinned.vcf")


### Transcriptome

transcriptome.genind <- vcfR2genind(transcriptome.vcf, ploidy = 4)
pop(transcriptome.genind) <- popmap$Temporal_Group
transcriptome.genind.nomissing <- missingno(transcriptome.genind, type = "mean")

### Transcriptome

pca.transcriptome <- dudi.pca(transcriptome.genind.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.transcriptome <- (pca.transcriptome$eig / sum(pca.transcriptome$eig))*100
pve.transcriptome <- round(pve.transcriptome, digits = 2)

pca.transcriptome.df <- pca.transcriptome$li
pca.transcriptome.df$Ind <- rownames(pca.transcriptome.df)
pca.transcriptome.df.popmap <- inner_join(pca.transcriptome.df, popmap)


### Transcriptome

# Collection year 

pca.txm.df.inv <- filter(pca.transcriptome.df.popmap, Invaded == "Y")
pca.txm.df.nat <- filter(pca.transcriptome.df.popmap, Invaded == "N")

pca.transcriptome.collectionyear <- ggplot() +
  geom_point(pca.txm.df.inv, mapping = aes(x = Axis1, y = Axis2, color = as.numeric(Collection.Year)), size = 3, alpha = 0.75) +
  geom_point(pca.txm.df.nat, mapping = aes(x = Axis1, y = Axis2), color = "black", shape = 17, size = 2) +
  scale_color_viridis_c(limits = range(pca.transcriptome.df.popmap$Collection.Year)) +
  xlab(paste0("PC 1 (", pve.transcriptome[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome[2], "% variation explained)")) +
  theme_bw() +
  labs(color = "Collection Year") +
  ggtitle("Transcriptome Reference")

ggsave("transcriptome.collectionyear.pca.pdf", height = 6, width = 8)
ggsave("transcriptome.collectionyear.pca.png", height = 6, width = 8, dpi = 300)

# Location

pca.transcriptome.location <- ggplot(pca.transcriptome.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(17, 16)) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.transcriptome[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome[2], "% variation explained)")) +
  ggtitle("Transcriptome Reference") +
  theme_bw()

ggsave("transcriptome.location.pca.pdf", height = 6, width = 8)
ggsave("transcriptome.location.pca.png", height = 6, width = 8, dpi = 300)



# Transcriptome

transcriptome.genind.noPI <- transcriptome.genind[!indNames(transcriptome.genind) %in% pacific.islands]
transcriptome.genind.noPI.nomissing <- missingno(transcriptome.genind.noPI, type = "mean")

### Transcriptome

pca.transcriptome.noPI <- dudi.pca(transcriptome.genind.noPI.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.transcriptome.noPI <- (pca.transcriptome.noPI$eig / sum(pca.transcriptome.noPI$eig))*100
pve.transcriptome.noPI <- round(pve.transcriptome.noPI, digits = 2)

pca.transcriptome.noPI.df <- pca.transcriptome.noPI$li
pca.transcriptome.noPI.df$Ind <- rownames(pca.transcriptome.noPI.df)
pca.transcriptome.noPI.df.popmap <- inner_join(pca.transcriptome.noPI.df, popmap)

### Transcriptome

# Collection Year

#pca.transcriptome.noPI.collectionYear <- ggplot(pca.transcriptome.noPI.df.popmap, aes(Axis1, Axis2, color = Collection.Year, shape = Invaded)) +
#  geom_point(size = 3) +
#  scale_color_viridis() +
#  xlab(paste0("PC 1 (", pve.transcriptome.noPI[1], "% variation explained)")) +
#  ylab(paste0("PC 2 (", pve.transcriptome.noPI[2], "% variation explained)")) +
#  ggtitle("Transcriptome") +
#  theme_bw()

#ggsave("transcriptome.noPI.collectionyear.pca.pdf", width = 8, height = 6)

# Location

pca.transcriptome.noPI.location <- ggplot(pca.transcriptome.noPI.df.popmap, aes(Axis1, Axis2, color = Location, shape = Invaded)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_paletteer_d("colorBlindness::paletteMartin") +
  xlab(paste0("PC 1 (", pve.transcriptome.noPI[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome.noPI[2], "% variation explained)")) +
  ggtitle("Transcriptome Reference") +
  theme_bw()

ggsave("transcriptome.noPI.location.pca.pdf", width = 8, height = 6)
ggsave("transcriptome.noPI.location.pca.png", width = 8, height = 6, dpi = 300)

