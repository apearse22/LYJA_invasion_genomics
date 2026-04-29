library(vcfR)
library(ggplot2)
library(adegenet)
library(dplyr)
library(readr)

### reading in files

captus <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")
popmap <- read.csv("files/popmap_updatedcoords.csv")

txm.captus.blast <- read.delim("files/LYJA.CAPTUS.blast.out")
colnames(txm.captus.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart",
                                "qend", "sstart", "send", "evalue", "bitscore")

txm.captus.blast.05eval <- txm.captus.blast %>% 
  filter(evalue <= 0.05) # forget to set evalue threshold when creating file in HPC

txm.arabidopsis.blast <- read.delim("files/Arabidopsis.LYJA.blast")
colnames(txm.arabidopsis.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                     "qstart",
                                "qend", "sstart", "send", "evalue", "bitscore")

# at some point, get best hit per locus 


################################# creating / subsetting genind objects ########################

# creating transcriptome / captus genind object

txm.genind <- vcfR2genind(txm, ploidy = 4)
pop(txm.genind) <- popmap$Temporal_Group

captus.genind <- vcfR2genind(captus, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group

# subsetting captus genind based on hits to transcriptome

locNames(captus.genind) <- gsub("_[0-9]..", "", locNames(captus.genind))
locNames(captus.genind) <- gsub("_[0-9].", "", locNames(captus.genind))

captus.coding.genind <- captus.genind[loc = locNames(captus.genind) %in% txm.captus.blast.05eval$qseqid]
captus.noncoding.genind <- captus.genind[loc = !locNames(captus.genind) %in% txm.captus.blast.05eval$qseqid]
