###########################################################

# Calculating individual and temporal group heterozygosity
# code developed: Abby Pearse, Jessie Pelosi
# Last updated: 03/25/2026

###########################################################

### Libraries

library(ggplot2)
library(dplyr)
library(vcfR)
library(adegenet)
library(poppr)
library(snpStats)
library(dartR.base)
library(ggpubr)

################################################### Reading in files and shared variables #########################################

popmap <- read.csv("data/popmap_updatedcoords.csv")

transcriptome.vcf <- read.vcfR("data/one_snp_per_locus_thinned.vcf")
captus.vcf <- read.vcfR("data/captus_0.5mising.maf0.5.thinned.vcf")

### Creating genind objects

# Transcriptome

transcriptome.genind <- vcfR2genind(transcriptome.vcf, ploidy = 4)
pop(transcriptome.genind) <- popmap$Temporal_Group

# CAPTUS

captus.genind <- vcfR2genind(captus.vcf, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group


#################################################### Calculating individual level heterozygosity #######################################

### Transcriptome

transcriptome.gl <- gi2gl(transcriptome.genind) # this changes the ploidy to 2 rather than keep 4 from genind object
ploidy(transcriptome.gl) <- 4 # set the ploidy to 4

gl.compliance.check(transcriptome.gl)

transcriptome.ho.ind <-  gl.report.polyploid_heterozygosity(transcriptome.gl, method = "ind", error.bar = "SE") 

### CAPTUS

captus.gl <- gi2gl(captus.genind)
ploidy(captus.gl) <- 4

captus.ho.ind <- gl.report.polyploid_heterozygosity(captus.gl, method = "ind", error.bar = "SE")

######################### Plotting individual level heterozygosity

### Transcriptome

colnames(transcriptome.ho.ind) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

transcriptome.ho.ind.popmap <- inner_join(transcriptome.ho.ind, popmap)

transcriptome.ho.ind.popmap.plot <- ggplot(transcriptome.ho.ind.popmap, aes(x = Collection.Year, y = Ho, color = Invaded)) + geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  ggtitle("Transcriptome Invididual Heterozygosity")

### CAPTUS

colnames(captus.ho.ind) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

captus.ho.ind.popmap <- inner_join(captus.ho.ind, popmap)

captus.ho.ind.popmap.plot <- ggplot(captus.ho.ind.popmap, aes(Collection.Year, Ho, color = Invaded)) + geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  stat_regline_equation(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  ggtitle("CAPTUS Individual Heterozygosity")

############################################# Calculating population (temporal group) level heterozygosity ###########################

transcriptome.ho.popmetric <- gl.report.polyploid_heterozygosity(transcriptome.gl, method = "pop", error.bar = "SE") # transcriptome
captus.ho.popmetric <- gl.report.polyploid_heterozygosity(captus.gl, method = "pop", error.bar = "SE") # CAPTUS

############### Plotting population level heterozygosity 

### Transcriptome

transcriptome.ho.popmetric.df <- data.frame(Temporal_Group = c("Early Invasion", "Mid Invasion", "Late Invasion", "Native"),
                                            Obs_Ho = transcriptome.ho.popmetric$Ho)

transcriptome.popmetric.plot.df <-  ggplot(transcriptome.ho.popmetric.df, aes(Temporal_Group, y = Obs_Ho, fill = Temporal_Group)) + 
  geom_col() +
  theme_bw() +
  ggtitle("Transcriptome Population (Temporal Group) Observed Heterozygosity") +
  scale_fill_manual(values = c ("Native" = "dodgerblue2", "Early Invasion" = "pink",
                                 "Mid Invasion" = "deeppink2", "Late Invasion" = "firebrick4")) +
  xlab("Temporal Group") +
  ylab("Observed Heterozygosity") 
  

### CAPTUS

captus.ho.popmetric.df <- data.frame(Temporal_Group = c("Early Invasion", "Mid Invasion", "Late Invasion", "Native"),
                                            Obs_Ho = captus.ho.popmetric$Ho)

captus.popmetric.plot.df <-  ggplot(captus.ho.popmetric.df, aes(Temporal_Group, y = Obs_Ho, fill = Temporal_Group)) + 
  geom_col() +
  theme_bw() +
  ggtitle("CAPTUS Population (Temporal Group) Observed Heterozygosity") +
  scale_fill_manual(values = c ("Native" = "dodgerblue2", "Early Invasion" = "pink",
                                "Mid Invasion" = "deeppink2", "Late Invasion" = "firebrick4")) +
  xlab("Temporal Group") +
  ylab("Observed Heterozygosity")


############################################# Conducting statistical analyses ###########################################

### Transcriptome

# Model with Interactions

transcriptome.Ho.lm.int <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = transcriptome.ho.ind.popmap)
summary(transcriptome.Ho.lm.int)


# Model without Interactions

transcriptome.Ho.lm.noint <- lm(Ho ~ Collection.Year + Invaded, data = transcriptome.ho.ind.popmap)
summary(transcriptome.Ho.lm.noint) 

### CAPTUS

# Model with Interactions

captus.Ho.lm.int <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = captus.ho.ind.popmap)
summary(captus.Ho.lm.int)

# Model without Interactions

captus.Ho.lm.noint <- lm(Ho ~ Collection.Year + Invaded, data = captus.ho.ind.popmap)
summary(captus.Ho.lm.noint) 

