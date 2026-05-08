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






############################################## Conducting statistical analyses ##############################

# standardize Collection Year to start at the beginning of the dataset

captus.coding.ho.popmap$std.Year <- captus.coding.ho.popmap$Collection.Year - min(captus.coding.ho.popmap$Collection.Year)

### coding loci

# interaction

coding.lm.interaction <- lm(Ho ~ std.Year + Invaded + std.Year:Invaded, data = captus.coding.ho.popmap)
summary(coding.lm.interaction)

# no interaction

coding.lm.nointeraction <- lm(Ho ~ std.Year + Invaded, data = captus.coding.ho.popmap)
summary(coding.lm.nointeraction)



r2_valu <- summary(coding.lm.nointeraction)$r.squared
f <- summary(coding.lm.nointeraction)$fstatistic

p_valu <- pf(f[1], f[2], f[3], lower.tail = F)

stats_label <- paste0("Adj.R^2 == ", round(r2_valu, 5),
                      "~~italic(P) == ", round(p_valu, 3))


ggplot(data = captus.coding.ho.popmap, mapping = aes(x = Collection.Year, y = Ho, color = Invaded, shape = Invaded)) + geom_point(size = 2) + 
  geom_smooth(method = "lm", linetype = "solid") + theme_bw() +
  #stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  annotate("text", 
           x = 1947, y = 0.23, 
           label = stats_label, 
           parse = T, 
           hjust = 1.1, vjust=1.5, size = 5) +
  scale_shape_manual(values = c(17, 16)) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  theme(text = element_text(size = 16), legend.position = "none") +
  ggtitle("B) Coding Loci")


ggsave("captus.coding.ho.plot.pdf", height = 6, width = 8)
ggsave("captus.coding.ho.plot.png", height = 6, width = 8, dpi = 300)

### non-coding

# standardize Collection Year to start at the beginning of the dataset

captus.noncoding.ho.popmap$std.Year <- captus.noncoding.ho.popmap$Collection.Year - min(captus.noncoding.ho.popmap$Collection.Year)


# interaction
noncoding.lm.interaction <- lm(Ho ~ std.Year + Invaded + std.Year:Invaded, data = captus.noncoding.ho.popmap)
summary(noncoding.lm.interaction)

captus.noncoding.ho.popmap.invaded <- filter(captus.noncoding.ho.popmap, Invaded == "Y")
captus.noncoding.ho.popmap.invaded.lm <- lm(Ho~std.Year, data = captus.noncoding.ho.popmap.invaded)
summary(captus.noncoding.ho.popmap.invaded.lm)

yr = data.frame(std.Year = 0)

predict(captus.noncoding.ho.popmap.invaded.lm, newdata = yr)

captus.noncoding.ho.popmap.native <- filter(captus.noncoding.ho.popmap, Invaded == "N")
captus.noncoding.ho.popmap.native.lm <- lm(Ho~std.Year, data = captus.noncoding.ho.popmap.native)
summary(captus.noncoding.ho.popmap.native.lm)

yr2 = data.frame(std.Year = 35)

predict(captus.noncoding.ho.popmap.native.lm, newdata = yr2)



stats_df <- captus.noncoding.ho.popmap %>% 
  group_by(Invaded) %>% 
  do(glance(lm(Ho ~ std.Year, data = .))) %>% 
  mutate(
    label = paste0("Adj.R^2 == ", round(r.squared, 5),
                   "~~italic(P) == ", round(p.value, 3)),
    y_pos = ifelse(Invaded == "N", 0.2, 0.195), 
    x_pos = 1945
  )
   

ggplot(data = captus.noncoding.ho.popmap, mapping = aes(x = Collection.Year, y = Ho, color = Invaded, shape = Invaded, linetype = Invaded)) + 
  geom_point(size = 2) + 
  geom_smooth(method = "lm") + theme_bw() +
  scale_linetype_manual(values = c("Y" = "solid",
                                   "N" = "dashed")) +
  #stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  geom_text(data = stats_df,
            aes(x = x_pos, y = y_pos, label = label, color = Invaded), 
            parse = T, hjust=1, size = 5, show.legend = F)+
  scale_shape_manual(values = c(17, 16)) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  theme(text = element_text(size = 16), legend.position = "none") +
  ggtitle("A) Non-Coding Loci")

ggsave("captus.noncoding.ho.plot.pdf", height = 6, width = 8)
ggsave("captus.noncoding.ho.plot.png",  height = 6, width = 8, dpi = 300)



# GEA loci -- see script with GEA analysis 





############################################# old analyses after this point! ###################################


####################### Transcriptome analyses #########################

#transcriptome.vcf <- read.vcfR("files/txm.50missing.CTmarked.recalc.maf0.05.thinned.vcf")

# Transcriptome

#transcriptome.genind <- vcfR2genind(transcriptome.vcf, ploidy = 4)
#pop(transcriptome.genind) <- popmap$Temporal_Group

#transcriptome.gl <- gi2gl(transcriptome.genind) # this changes the ploidy to 2 rather than keep 4 from genind object
#ploidy(transcriptome.gl) <- 4 # set the ploidy to 4

#gl.compliance.check(transcriptome.gl)

#transcriptome.ho.ind <-  gl.report.polyploid_heterozygosity(transcriptome.gl, method = "ind", error.bar = "SE") 

### Transcriptome

#colnames(transcriptome.ho.ind) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

#transcriptome.ho.ind.popmap <- inner_join(transcriptome.ho.ind, popmap)

#transcriptome.ho.ind.popmap.plot <- ggplot(transcriptome.ho.ind.popmap, aes(x = Collection.Year, y = Ho, color = Invaded)) + 
#  geom_point() + 
#  geom_smooth(method = "lm") + 
#  theme_bw() +
#  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
#  ggtitle("Transcriptome Reference - All Loci") +
#  xlab("Collection Year") +
#  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154"))

#ggsave("transcriptome.ho.ind.popmap.plot.pdf", transcriptome.ho.ind.popmap.plot, height = 6, width = 8)
#ggsave("transcriptome.ho.ind.popmap.plot.png", transcriptome.ho.ind.popmap.plot, height = 6, width = 8, dpi = 300)

#transcriptome.ho.popmetric <- gl.report.polyploid_heterozygosity(transcriptome.gl, method = "pop", error.bar = "SE") # transcriptome

### Transcriptome

#transcriptome.ho.popmetric.df <- data.frame(Temporal_Group = c("Early Invasion", "Mid Invasion", "Late Invasion", "Native"),
#                                            Obs_Ho = transcriptome.ho.popmetric$Ho)

#transcriptome.popmetric.plot.df <-  ggplot(transcriptome.ho.popmetric.df, aes(Temporal_Group, y = Obs_Ho, fill = Temporal_Group)) + 
#  geom_col() +
#  theme_bw() +
#  ggtitle("Transcriptome Population (Temporal Group) Observed Heterozygosity") +
#  scale_fill_manual(values = c ("Native" = "dodgerblue2", "Early Invasion" = "pink",
#                                "Mid Invasion" = "deeppink2", "Late Invasion" = "firebrick4")) +
#  xlab("Temporal Group") +
#  ylab("Observed Heterozygosity") 

### Transcriptome

# Model with Interactions

#transcriptome.Ho.lm.int <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = transcriptome.ho.ind.popmap)
#summary(transcriptome.Ho.lm.int)


# Model without Interactions

#transcriptome.Ho.lm.noint <- lm(Ho ~ Collection.Year + Invaded, data = transcriptome.ho.ind.popmap)
#summary(transcriptome.Ho.lm.noint) 


############################################# Calculating population (temporal group) level heterozygosity ###########################

#captus.ho.popmetric <- gl.report.polyploid_heterozygosity(captus.gl, method = "pop", error.bar = "SE") # CAPTUS

############### Plotting population level heterozygosity 

#captus.ho.popmetric.df <- data.frame(Temporal_Group = c("Early Invasion", "Mid Invasion", "Late Invasion", "Native"),
#                                     Obs_Ho = captus.ho.popmetric$Ho)

#captus.popmetric.plot.df <-  ggplot(captus.ho.popmetric.df, aes(Temporal_Group, y = Obs_Ho, fill = Temporal_Group)) + 
#  geom_col() +
#  theme_bw() +
#  ggtitle("CAPTUS Population (Temporal Group) Observed Heterozygosity") +
#  scale_fill_manual(values = c ("Native" = "dodgerblue2", "Early Invasion" = "pink",
#                                "Mid Invasion" = "deeppink2", "Late Invasion" = "firebrick4")) +
#  xlab("Temporal Group") +
#  ylab("Observed Heterozygosity")


### calculating captus statistical models


##### Conducting statistical analyses

### CAPTUS

# Model with Interactions

#captus.Ho.lm.int <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = captus.ho.ind.popmap)
#summary(captus.Ho.lm.int)

# Model without Interactions

#captus.Ho.lm.noint <- lm(Ho ~ Collection.Year + Invaded, data = captus.ho.ind.popmap)
#summary(captus.Ho.lm.noint) 
