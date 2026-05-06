###########################################################

# Calculating coding and non-coding heterozygosity
# code developed: Abby Pearse, Jessie Pelosi
# Last updated: 04/27/2026

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

captus <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")
popmap <- read.csv("files/popmap_updatedcoords.csv")

txm.captus.blast <- read.delim("files/LYJA.CAPTUS.blast.out")
colnames(txm.captus.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart",
                                "qend", "sstart", "send", "evalue", "bitscore")

txm.captus.blast.05eval <- txm.captus.blast %>% 
  filter(evalue <= 1e-5) # forget to set evalue threshold when creating file in HPC


################################# creating / subsetting genind objects ########################

# creating captus genind object

captus.genind <- vcfR2genind(captus, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group

# subsetting captus genind based on hits to transcriptome

locNames(captus.genind) <- gsub("_[0-9]..", "", locNames(captus.genind))
locNames(captus.genind) <- gsub("_[0-9].", "", locNames(captus.genind))

captus.coding.genind <- captus.genind[loc = locNames(captus.genind) %in% txm.captus.blast.05eval$qseqid]
captus.noncoding.genind <- captus.genind[loc = !locNames(captus.genind) %in% txm.captus.blast.05eval$qseqid]



######################################### Calculating individual level heterozygosity #######################################

### coding

coding.gl <- gi2gl(captus.coding.genind)
ploidy(coding.gl) <- 4

captus.coding.ho <- gl.report.polyploid_heterozygosity(coding.gl, method = "ind", error.bar = "SE")


### non-coding

noncoding.gl <- gi2gl(captus.noncoding.genind)
ploidy(noncoding.gl) <- 4

captus.noncoding.ho <- gl.report.polyploid_heterozygosity(noncoding.gl, method = "ind", error.bar = "SE")


################################## Plotting individual level heterozygosity ###############################

### coding
colnames(captus.coding.ho) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

captus.coding.ho.popmap <- inner_join(captus.coding.ho, popmap)

captus.coding.ho.popmap.plot <- ggplot(captus.coding.ho.popmap, aes(Collection.Year, Ho, color = Invaded, shape = Invaded)) + geom_point(size = 2) +
  theme_bw() +
  geom_smooth(method = "lm", linetype = "dashed") +
  stat_regline_equation(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  scale_shape_manual(values = c(17, 16)) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  theme(text = element_text(size = 16))

ggsave("captus.coding.ho.plot.pdf", captus.coding.ho.popmap.plot, height = 6, width = 8)
ggsave("captus.coding.ho.plot.png", captus.coding.ho.popmap.plot, height = 6, width = 8, dpi = 300)


### non-coding

colnames(captus.noncoding.ho) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

captus.noncoding.ho.popmap <- inner_join(captus.noncoding.ho, popmap)

captus.noncoding.ho.popmap.plot <- ggplot(captus.noncoding.ho.popmap, aes(Collection.Year, Ho, color = Invaded, shape = Invaded, linetype = Invaded)) + geom_point(size = 2) +
  theme_bw() +
  geom_smooth(method = "lm") +
  scale_linetype_manual(values = c("Y" = "solid",
                                   "N" = "dashed")) +
  scale_shape_manual(values = c(17, 16)) +
  stat_regline_equation(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  theme(text = element_text(size = 16)) #+ xlim(1940, 2015)

ggsave("captus.noncoding.ho.plot.pdf", captus.noncoding.ho.popmap.plot, height = 6, width = 8)
ggsave("captus.noncoding.ho.plot.png", captus.noncoding.ho.popmap.plot, height = 6, width = 8, dpi = 300)





############################################## Conducting statistical analyses ##############################

# standardize Collection Year to start at the beginning of the dataset

captus.coding.ho.popmap$std.Year <- captus.coding.ho.popmap$Collection.Year - min(captus.coding.ho.popmap$Collection.Year)

### coding loci

# interaction

coding.lm.interaction <- lm(Ho ~ std.Year + Invaded + std.Year:Invaded, data = captus.coding.ho.popmap)
summary(coding.lm.interaction)

# no interaction

coding.lm.nointeraction <- lm(Ho ~ Collection.Year + Invaded, data = captus.coding.ho.popmap)
summary(coding.lm.nointeraction)

### non-coding

# standardize Collection Year to start at the beginning of the dataset

captus.noncoding.ho.popmap$std.Year <- captus.noncoding.ho.popmap$Collection.Year - min(captus.noncoding.ho.popmap$Collection.Year)


# interaction
noncoding.lm.interaction <- lm(Ho ~ std.Year + Invaded + std.Year:Invaded, data = captus.noncoding.ho.popmap)
summary(noncoding.lm.interaction)

# standardize: year - minimum year 
# run a model for just the invaded range samples lm(Ho~Collection.Year, data = invaded.samples)

# no effect of year on Ho in the native samples
# across all years, there is an effect of invasion status on Ho --> invaded samples have lower Ho than native samples 
#   at the year 0, this is the y-intercept  (5.91e-01)
# and there is an additional effect of collection year only in the invaded samples 
# maginally significant interaction between collection year and invasion status, increased Ho as year increases 

t <- filter(captus.noncoding.ho.popmap, Invaded == "Y")
test <- lm(Ho~Collection.Year, data = t)
summary(test)

# no interaction

noncoding.lm.nointeraction <- lm(Ho ~ Collection.Year + Invaded, data = captus.noncoding.ho.popmap)
summary(noncoding.lm.nointeraction)












############################################# old analyses after this point! ###################################


####################### Transcriptome analyses #########################

transcriptome.vcf <- read.vcfR("files/txm.50missing.CTmarked.recalc.maf0.05.thinned.vcf")

# Transcriptome

transcriptome.genind <- vcfR2genind(transcriptome.vcf, ploidy = 4)
pop(transcriptome.genind) <- popmap$Temporal_Group

transcriptome.gl <- gi2gl(transcriptome.genind) # this changes the ploidy to 2 rather than keep 4 from genind object
ploidy(transcriptome.gl) <- 4 # set the ploidy to 4

#gl.compliance.check(transcriptome.gl)

transcriptome.ho.ind <-  gl.report.polyploid_heterozygosity(transcriptome.gl, method = "ind", error.bar = "SE") 

### Transcriptome

colnames(transcriptome.ho.ind) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

transcriptome.ho.ind.popmap <- inner_join(transcriptome.ho.ind, popmap)

transcriptome.ho.ind.popmap.plot <- ggplot(transcriptome.ho.ind.popmap, aes(x = Collection.Year, y = Ho, color = Invaded)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  ggtitle("Transcriptome Reference - All Loci") +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154"))

ggsave("transcriptome.ho.ind.popmap.plot.pdf", transcriptome.ho.ind.popmap.plot, height = 6, width = 8)
ggsave("transcriptome.ho.ind.popmap.plot.png", transcriptome.ho.ind.popmap.plot, height = 6, width = 8, dpi = 300)

transcriptome.ho.popmetric <- gl.report.polyploid_heterozygosity(transcriptome.gl, method = "pop", error.bar = "SE") # transcriptome

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

### Transcriptome

# Model with Interactions

transcriptome.Ho.lm.int <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = transcriptome.ho.ind.popmap)
summary(transcriptome.Ho.lm.int)


# Model without Interactions

transcriptome.Ho.lm.noint <- lm(Ho ~ Collection.Year + Invaded, data = transcriptome.ho.ind.popmap)
summary(transcriptome.Ho.lm.noint) 


############################################# Calculating population (temporal group) level heterozygosity ###########################

captus.ho.popmetric <- gl.report.polyploid_heterozygosity(captus.gl, method = "pop", error.bar = "SE") # CAPTUS

############### Plotting population level heterozygosity 

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


### calculating captus statistical models


##### Conducting statistical analyses

### CAPTUS

# Model with Interactions

captus.Ho.lm.int <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = captus.ho.ind.popmap)
summary(captus.Ho.lm.int)

# Model without Interactions

captus.Ho.lm.noint <- lm(Ho ~ Collection.Year + Invaded, data = captus.ho.ind.popmap)
summary(captus.Ho.lm.noint) 
