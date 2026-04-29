######################################################################################################

# Subsetting climate-associated genes (Whiting et al. 2024) from the CAPTUS coding and non-coding loci
# code developed by: Abby Pearse, Jessie Pelosi
# Last updated: 04/13/2026

######################################################################################################

### Libraries

library(vcfR)
library(ggplot2)
library(dplyr)
library(dartR.base)
library(adegenet)

### Reading in files

captus <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")

popmap <- read.csv("files/popmap_updatedcoords.csv", sep = ",")

Whiting.blast.hits <- read.delim("files/LYJA_Whiting_blast.out", header = F, sep = "")

txm.captus.blast <- read.delim("files/LYJA.CAPTUS.blast.out")
colnames(txm.captus.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart",
                                "qend", "sstart", "send", "evalue", "bitscore")

txm.captus.blast.05eval <- txm.captus.blast %>% 
  filter(evalue <= 1e-5) # forget to set evalue threshold when creating file in HPC

txm.arabidopsis.blast <- read.delim("files/Arabidopsis.LYJA.blast")
colnames(txm.arabidopsis.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                     "qstart",
                                     "qend", "sstart", "send", "evalue", "bitscore")
txm.arabidopsis.blast.05eval <- txm.arabidopsis.blast %>% 
  filter(evalue <= 1e-5) # forget to set evalue threshold when creating file in HPC


###################################### subsetting CAPTUS for coding and non-coding loci ################

captus.genind <- vcfR2genind(captus, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group

# subsetting captus genind based on hits to transcriptome

locNames(captus.genind) <- gsub("_[0-9]..", "", locNames(captus.genind))
locNames(captus.genind) <- gsub("_[0-9].", "", locNames(captus.genind))

captus.coding.genind <- captus.genind[loc = locNames(captus.genind) %in% txm.captus.blast.05eval$qseqid]
captus.noncoding.genind <- captus.genind[loc = !locNames(captus.genind) %in% txm.captus.blast.05eval$qseqid]


#################################### joining Arabidopsis genes with captus genes ########################

### finding top blast hit per query

best.hit.arabidopsis.txm <- txm.arabidopsis.blast %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(row_number() == 1)


colnames(best.hit.arabidopsis.txm) <- c("ID", "sseqid", "pident", "length", "mismatch", "gapopen",
                                        "qstart", "qend", "sstart", "send", "evalue", "bitscore")


best.hit.captus.txm <- txm.captus.blast.05eval %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(row_number() == 1)

colnames(best.hit.captus.txm) <- c("qseqid", "ID", "pident", "length", "mismatch", "gapopen",
                                   "qstart", "qend", "sstart", "send", "evalue", "bitscore")

######################### subsetting geninds for coding and whiting genes ##################################

captus.arabidopsis.blast <- inner_join(best.hit.arabidopsis.txm, best.hit.captus.txm, by = "ID", suffix = c(".arabidopsis", ".captus"))


coding.arabidopsis <- captus.coding.genind[loc = locNames(captus.coding.genind) %in% captus.arabidopsis.blast$qseqid]
noncoding.arabidopsis <- captus.noncoding.genind[loc = locNames(captus.noncoding.genind) %in% captus.arabidopsis.blast$qseqid] # 0 - think this is just bc there are no noncoding hits


arabiopsis.in.captus.df <- captus.arabidopsis.blast %>% 
  filter(qseqid %in% locNames(coding.arabidopsis))


coding.in.whiting.genind <- coding.arabidopsis[loc = locNames(coding.arabidopsis) %in% captus.arabidopsis.blast$qseqid]



#colnames(Whiting.blast.hits) <- c("query", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "ssend", 
 #                                 "evalue", "bitscore")

#coding.in.whiting <- inner_join(Whiting.blast.hits, arabiopsis.in.captus.df, by = "sseqid")


################################ calculating heterozygosity for whiting genes in coding loci ###################

gl.coding.in.whiting <- gi2gl(coding.in.whiting.genind)
ploidy(gl.coding.in.whiting) <- 4
gl.compliance.check(gl.coding.in.whiting) 


coding.in.whiting.ho <- gl.report.polyploid_heterozygosity(gl.coding.in.whiting, method = "ind")
colnames(coding.in.whiting.ho) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

coding.in.whiting.ho.popmap <- inner_join(coding.in.whiting.ho, popmap)


### plotting
coding.in.whiting.plot <- ggplot(coding.in.whiting.ho.popmap, aes(Collection.Year, Ho, color = Invaded, shape = Invaded)) + geom_point(size = 2) +
  theme_bw() +
  geom_smooth(method = "lm") +
  scale_shape_manual(values = c(17, 16)) +
  stat_regline_equation(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  theme(text = element_text(size = 16))



############################################## statistical tests! #############################

# interaction

coding.in.whiting.interaction <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = coding.in.whiting.ho.popmap)
summary(coding.in.whiting.interaction)

# no interaction

coding.in.whiting.nointeraction <- lm(Ho ~ Collection.Year + Invaded, data = coding.in.whiting.ho.popmap)
summary(coding.in.whiting.nointeraction)







############################## Subsetting transcriptome genind for Whiting et al. genes #########################

all.loci <- as.data.frame(transcriptome.genind@loc.fac)
all.loci.uniq <- distinct(all.loci)

all.loci.uniq$query <- gsub("_[0-9]+","", all.loci.uniq$`transcriptome.genind@loc.fac`)

colnames(blast.hits) <- c("query", "subject", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "ssend", 
                     "evalue", "bitscore")

blast.loci<- inner_join(all.loci.uniq, blast.hits)
blast.loci.uniq <- distinct(blast.loci, `transcriptome.genind@loc.fac`) 


############################# Creating genind of climate genes in the transcriptome ############################

genind.climate <- transcriptome.genind[loc = blast.loci.uniq$`transcriptome.genind@loc.fac`]

gl.climate <- dartR.base::gi2gl(genind.climate)

ploidy(gl.climate) <- 4 # set the ploidy to 4
gl.compliance.check(gl.climate)


climate.ho <- gl.report.polyploid_heterozygosity(gl.climate, method = "ind")
colnames(climate.ho) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

climate.ho.popmap <- inner_join(climate.ho, popmap)

ggplot(data = climate.ho.popmap, mapping = aes(x = Collection.Year, y = Ho, color = Invaded)) + geom_point() + 
  geom_smooth(method = "lm") + theme_bw() +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")))




########################################## calculating heterozygosity ####################################

Ho.climate <- gl.report.polyploid_heterozygosity(gl.climate, method = "ind")
Ho.climate.popmetric <- gl.report.polyploid_heterozygosity(gl.climate, method = "pop") # temporal group heterozygosity



colnames(Ho.climate) <- c("Ind", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")

Ho.climate.popmap <- inner_join(Ho.climate, popmap)

