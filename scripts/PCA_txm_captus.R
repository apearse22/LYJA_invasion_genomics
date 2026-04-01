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


############################################ Generating PCAs ###########################################

### Transcriptome

pca.transcriptome <- dudi.pca(transcriptome.genind.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.transcriptome <- (pca.transcriptome$eig / sum(pca.transcriptome$eig))*100
pve.transcriptome <- round(pve.transcriptome, digits = 2)

pca.transcriptome.df <- pca.transcriptome$li
pca.transcriptome.df$Ind <- rownames(pca.transcriptome.df)
pca.transcriptome.df.popmap <- inner_join(pca.transcriptome.df, popmap)


pca.transcriptome.tempgroups <- ggplot(pca.transcriptome.df.popmap, aes(Axis1, Axis2, color = Temporal_Group)) + geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC 1 (", pve.transcriptome[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.transcriptome[2], "% variation explained)")) +
  scale_color_manual(values = c ("native" = "dodgerblue2", "early invasion" = "pink",
                                 "mid invasion" = "deeppink2", "late invasion" = "firebrick4"))+
  ggtitle("TRANSCRIPTOME PCA") 

#ggsave("transcriptome.temporalgroups.pdf", width = 8, height = 6)

### CAPTUS

pca.captus <- dudi.pca(captus.genind.nomissing, scale = FALSE, scannf = FALSE, nf = 50)

pve.captus <- (pca.captus$eig / sum(pca.captus$eig))*100
pve.captus <- round(pve.captus, digits = 2)

pca.captus.df <- pca.captus$li
pca.captus.df$Ind <- rownames(pca.captus.df)
pca.captus.df.popmap <- inner_join(pca.captus.df, popmap)

pca.captus.tempgroups <- ggplot(pca.captus.df.popmap, aes(Axis1, Axis2, color = Temporal_Group)) + geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC 1 (", pve.captus[1], "% variation explained)")) +
  ylab(paste0("PC 2 (", pve.captus[2], "% variation explained)")) +
  scale_color_manual(values = c ("native" = "dodgerblue2", "early invasion" = "pink",
                                 "mid invasion" = "deeppink2", "late invasion" = "firebrick4")) +
  ggtitle("CAPTUS PCA") 

#ggsave("captus.pca.pdf", height = 6, width = 8)