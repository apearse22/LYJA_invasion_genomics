#####################################################################

# Correlating genetic and geographic distance for genomic references
# Code developed by Abby Pearse, Jessie Pelosi
# Last updated: 03/25/2026

#####################################################################

### Libraries

library(vegan)
library(vcfR)
library(dplyr)
library(ggplot2)

################################################ Reading in files and shared variables #######################################

samples.native <- c("SRR29127787", "SRR29127784", "SRR29127767","SRR29127785","SRR29127770","SRR29127789","SRR29127769","SRR29127786","SRR29127774","SRR29127777","SRR29127780",
                    "SRR29127782","SRR29127788","SRR29127783","SRR29127776","SRR29127775","SRR29127771","SRR29127764","SRR29127766",
                    "SRR29127773","SRR29127781","SRR29127778","SRR29127765", "SRR29127793", "SRR29127772")

popmap <- read.csv("files/popmap_updatedcoords.csv")
captus <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")


######################################### Creating population maps of invaded and native samples ####################################

# Invaded

popmap.inv <- popmap[!popmap$Ind %in% samples.native, ]
popmap.inv <- arrange(popmap.inv, Ind)

popmap.inv.mat <- as.data.frame(popmap.inv, row.names = popmap.inv$Ind)

# Native

popmap.nat <- popmap[popmap$Ind %in% samples.native, ]
popmap.nat <- arrange(popmap.nat, Ind)

popmap.nat.mat <- as.data.frame(popmap.nat, row.names = popmap.nat$Ind)


############################################## Obtaining sample coordinates and collection years ###################################

### Coordinates

# Invaded samples

popmap.inv.coords <- popmap.inv.mat %>% 
  dplyr::select(Approx..Latitude, Approx..Longitude) %>% 
  na.omit()

# Native samples

popmap.nat.coords <- popmap.nat.mat %>% 
  dplyr::select(Approx..Latitude, Approx..Longitude) %>% 
  na.omit()

### Collection year

# Invaded samples

popmap.inv.years <- popmap.inv.mat %>% 
  dplyr::select(Collection.Year) 

# Native samples

popmap.nat.years <- popmap.nat.mat %>% 
  dplyr::select(Collection.Year)


############################################################ Creating genind objects #############################################

captus.genind <- vcfR2genind(captus, ploidy = 4)

captus.genind.inv <- captus.genind[!indNames(captus.genind) %in% samples.native]
captus.genind.nat <- captus.genind[indNames(captus.genind) %in% samples.native & indNames(transcriptome.genind) != "SRR29127772"]

# same logic as above
captus.genind.nat.temp <- captus.genind[indNames(captus.genind) %in% samples.native]

################################################# Calculating genetic distance matrices ###############################

### CAPTUS invaded samples

captus.gen.dist.inv <- prevosti.dist(captus.genind.inv)

### CAPTUS native samples

captus.gen.dist.nat <- prevosti.dist(captus.genind.nat)
captus.gen.dist.nat.temp <- prevosti.dist(captus.genind.nat.temp) # temporal distance, including SRR29127772


###################################################### Isolation by distance ##############################################

### Creating geographic distance matrices - these will be the same for the transcriptome and cAPTUS

geo.dist.inv <- dist(popmap.inv.coords, method = "euclidean") # Invaded samples

geo.dist.nat <- dist(popmap.nat.coords, method = "euclidean") # Native samples


########### Conducting isolation by distance Mantel tests

# CAPTUS invaded samples

mantel.captus.inv <- mantel(captus.gen.dist.inv, geo.dist.inv, method = "pearson", permutations = 999)
summary(mantel.captus.inv)

# CAPTUS native samples

mantel.captus.nat <- mantel(captus.gen.dist.nat, geo.dist.nat, method = "pearson", permutations = 999)
summary(mantel.captus.nat)


########### Plotting isolation by distance

geo.dist.inv.vec <- as.vector(geo.dist.inv)
geo.dist.nat.vec <- as.vector(geo.dist.nat)


### CAPTUS invaded samples

captus.gen.dist.inv.vec <- as.vector(captus.gen.dist.inv)

captus.inv.geo.plot.df <- data.frame(
  GenDist = captus.gen.dist.inv.vec, 
  GeoDist = geo.dist.inv.vec)

captus.inv.geo.plot <- ggplot(captus.inv.geo.plot.df, aes(GeoDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Geographic Distance") +
  ylab("Genetic Distance") +
  ggtitle("Invaded CAPTUS Isolation by Distance") +
  theme_bw()


### CAPTUS native samples

captus.gen.dist.nat.vec <- as.vector(captus.gen.dist.nat)

captus.nat.geo.plot.df <- data.frame(
  GenDist = captus.gen.dist.nat.vec, 
  GeoDist = geo.dist.nat.vec)

captus.nat.geo.plot <- ggplot(captus.nat.geo.plot.df, aes(GeoDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Geographic Distance") +
  ylab("Genetic Distance") +
  ggtitle("Native CAPTUS Isolation by Distance") +
  theme_bw()


# Saving plot images

#ggsave("invaded.transcriptome.geo.pdf", plot = transcriptome.inv.geo.plot, height = 7, width = 7)
#ggsave("native.transcriptome.geo.pdf", plot = transcriptome.nat.geo.plot, height = 7, width = 7)
#ggsave("invaded.captus.geo.pdf", plot = captus.inv.geo.plot, height = 7, width = 7)
#ggsave("native.captus.geo.pdf", plot = captus.nat.geo.plot, height = 7, width = 7)

################################################### Isolation by time ##################################


temp.dist.inv <- dist(popmap.inv.years$Collection.Year, method = "euclidean") # Invasive temporal matrix
temp.dist.nat <- dist(popmap.nat.years$Collection.Year, method = "euclidean") # Native temporal matrix

temp.dist.inv.vec <- as.vector(temp.dist.inv)
temp.dist.nat.vec <- as.vector(temp.dist.nat)

########### Conducting Mantel tests

# CAPTUS invaded samples

mantel.captus.inv.temp <- mantel(captus.gen.dist.inv, temp.dist.inv, method = "pearson", permutations = 999)
summary(mantel.captus.inv.temp)


# CAPTUS native samples

mantel.captus.nat.temp <- mantel(captus.gen.dist.nat.temp, temp.dist.nat, method = "pearson", permutations = 999)
summary(mantel.captus.nat.temp)


########### Plotting temporal isolation by distance

### Invaded CAPTUS

captus.inv.temp.plot.df <- data.frame(
  GenDist = captus.gen.dist.inv.vec, 
  TempDist = temp.dist.inv.vec)

captus.inv.temp.plot <- ggplot(captus.inv.temp.plot.df, aes(TempDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Temporal Distance") +
  ylab("Genetic Distance") +
  ggtitle("Invaded CAPTUS Isolation by Time") +
  theme_bw()


### Native CAPTUS

captus.gen.dist.nat.temp.vec <- as.vector(captus.gen.dist.nat.temp)

captus.nat.temp.plot.df <- data.frame(
  GenDist = captus.gen.dist.nat.temp.vec, 
  TempDist = temp.dist.nat.vec)

captus.nat.temp.plot <- ggplot(captus.nat.temp.plot.df, aes(TempDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Temporal Distance") +
  ylab("Genetic Distance") +
  ggtitle("Native CAPTUS Isolation by Time") +
  theme_bw()


# Saving plot images

#ggsave("invaded.transcriptome.temp.pdf", plot = transcriptome.inv.temp.plot, height = 7, width = 7)
#ggsave("native.transcriptome.temp.pdf", plot = transcriptome.nat.temp.plot, height = 7, width = 7)
#ggsave("invaded.captus.temp.pdf", plot = captus.inv.temp.plot, height = 7, width = 7)
#ggsave("native.captus.temp.pdf", plot = captus.nat.temp.plot, height = 7, width = 7)


######################################### Genetic distance of plastid haplotypes #############################################

cp1 <- c("SRR29127726", "SRR29127727", "SRR29127728", "SRR29127729", "SRR29127730", "SRR29127731", 
         "SRR29127732", "SRR29127733", "SRR29127734", "SRR29127735", "SRR29127736", "SRR29127737", 
         "SRR29127738", "SRR29127739", "SRR29127740", "SRR29127741", "SRR29127779", "SRR29127743", 
         "SRR29127744", "SRR29127748", "SRR29127742", "SRR29127747", "SRR29127749", "SRR29127750", 
         "SRR29127751", "SRR29127752", "SRR29127753", "SRR29127754", "SRR29127755", "SRR29127756", 
         "SRR29127757", "SRR29127758", "SRR29127759", "SRR29127760", "SRR29127761", "SRR29127762", 
         "SRR29127763", "SRR29127790", "SRR29127791", "SRR29127792", "SRR29127793", "SRR29127794", 
         "SRR29127795", "SRR29127796", "SRR29127797", "SRR29127798", "SRR29127799")


cp2 <- c("SRR29127768")
cp3 <- c("SRR29127745", "SRR29127746")


### CAPTUS

captus.inv <- as.matrix(captus.gen.dist.inv)
captus.inv[upper.tri(captus.inv, diag = FALSE)] <- NA
dimnames(captus.inv) <- list(rownames(captus.inv), colnames(captus.inv))

captus.inv.df <- as.data.frame(as.table(captus.inv))
captus.inv.df <- na.omit(captus.inv.df)

captus.inv.df <- captus.inv.df %>% 
  mutate(plastid_haplotype_1 = case_when(
    Var1 %in% cp1 ~ "1",
    Var1 %in% cp2 ~ "2",
    Var1 %in% cp3 ~ "3",
    TRUE            ~ "Not Found" 
  )) %>% 
  mutate(plastid_haplotype_2 = case_when (
    Var2 %in% cp1 ~ "1", 
    Var2 %in% cp2 ~ "2",
    Var2 %in% cp3 ~ "3",
    TRUE            ~ "Not Found" 
  )) %>% 
  mutate(haplotype_div = case_when (
    plastid_haplotype_1 == plastid_haplotype_2 ~ "same", 
    plastid_haplotype_1 != plastid_haplotype_2 ~ "different"
  )) %>% 
  filter(Freq != 0) %>% 
  group_by(haplotype_div)

captus.cphaps.plot <- ggplot(captus.inv.df, aes(x = haplotype_div, y = Freq, fill = haplotype_div)) + 
  geom_boxplot(alpha = 0.5) +
  xlab("Same/Different Plastid Haplotype") + 
  ylab("Genetic Distance") + 
  ggtitle("CAPTUS") +
  theme_bw()

wilcox.test(Freq ~ haplotype_div, data = captus.inv.df)





############################################################# TRANSCRIPTOME ANALYSES ##############################################



transcriptome <- read.vcfR("files/txm.50missing.CTmarked.recalc.maf0.05.thinned.vcf")
### Transcriptome

transcriptome.genind <- vcfR2genind(transcriptome, ploidy = 4)

transcriptome.genind.inv <- transcriptome.genind[!indNames(transcriptome.genind) %in% samples.native]
transcriptome.genind.nat <- transcriptome.genind[indNames(transcriptome.genind) %in% samples.native & indNames(transcriptome.genind) != "SRR29127772"]

# SRR...72 only needs to be removed for the geographic isolation because it has no coordinates, it's fine to include for temporal isolation. making a new genind

transcriptome.genind.nat.temp <- transcriptome.genind[indNames(transcriptome.genind) %in% samples.native]

### Transcriptome invaded samples

transcriptome.gen.dist.inv <- prevosti.dist(transcriptome.genind.inv)

### Transcriptome native samples

transcriptome.gen.dist.nat <- prevosti.dist(transcriptome.genind.nat)
transcriptome.gen.dist.nat.temp <- prevosti.dist(transcriptome.genind.nat.temp) # temporal distance, including SRR29127772

# Transcriptome invaded samples

mantel.transcriptome.inv <- mantel(transcriptome.gen.dist.inv, geo.dist.inv, method = "pearson", permutations = 999)
summary(mantel.transcriptome.inv)

# Transcriptome native samples

mantel.transcriptome.nat <- mantel(transcriptome.gen.dist.nat, geo.dist.nat, method = "pearson", permutations = 999)
summary(mantel.transcriptome.nat)

### Transcriptome invaded samples

transcriptome.gen.dist.inv.vec <- as.vector(transcriptome.gen.dist.inv)

transcriptome.inv.geo.plot.df <- data.frame(
  GenDist = transcriptome.gen.dist.inv.vec, 
  GeoDist = geo.dist.inv.vec)


transcriptome.inv.geo.plot <- ggplot(transcriptome.inv.geo.plot.df, aes(GeoDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Geographic Distance") +
  ylab("Genetic Distance") +
  ggtitle("Invaded Transcriptome Isolation by Distance") +
  theme_bw()


### Transcriptome native samples

transcriptome.gen.dist.nat.vec <- as.vector(transcriptome.gen.dist.nat)

transcriptome.nat.geo.plot.df <- data.frame(
  GenDist = transcriptome.gen.dist.nat.vec, 
  GeoDist = geo.dist.nat.vec)

transcriptome.nat.geo.plot <- ggplot(transcriptome.nat.geo.plot.df, aes(GeoDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Geographic Distance") +
  ylab("Genetic Distance") +
  ggtitle("Native Transcriptome Isolation by Distance") +
  theme_bw()

# Transcriptome invaded samples

mantel.transcriptome.inv.temp <- mantel(transcriptome.gen.dist.inv, temp.dist.inv, method = "pearson", permutations = 999)
summary(mantel.transcriptome.inv.temp)

# Transcriptome native samples

mantel.transcriptome.nat.temp <- mantel(transcriptome.gen.dist.nat.temp, temp.dist.nat, method = "pearson", permutations = 999)
summary(mantel.transcriptome.nat.temp)


### Invaded transcriptome

transcriptome.inv.temp.plot.df <- data.frame(
  GenDist = transcriptome.gen.dist.inv.vec, 
  TempDist = temp.dist.inv.vec)

transcriptome.inv.temp.plot <- ggplot(transcriptome.inv.temp.plot.df, aes(TempDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Temporal Distance") +
  ylab("Genetic Distance") +
  ggtitle("Invaded Transcriptome Isolation by Time") +
  theme_bw()


### Native transcriptome

transcriptome.gen.dist.nat.temp.vec <- as.vector(transcriptome.gen.dist.nat.temp)

transcriptome.nat.temp.plot.df <- data.frame(
  GenDist = transcriptome.gen.dist.nat.temp.vec, 
  TempDist = temp.dist.nat.vec)

transcriptome.nat.temp.plot <- ggplot(transcriptome.nat.temp.plot.df, aes(TempDist, GenDist)) +
  geom_point(color = "royalblue") + 
  geom_smooth(method = "lm", color = "red") +
  xlab("Temporal Distance") +
  ylab("Genetic Distance") +
  ggtitle("Native Transcriptome Isolation by Time") +
  theme_bw()

### Transcriptome

transcriptome.inv <- as.matrix(transcriptome.gen.dist.inv)
transcriptome.inv[upper.tri(transcriptome.inv, diag = FALSE)] <- NA
dimnames(transcriptome.inv) <- list(rownames(transcriptome.inv), colnames(transcriptome.inv))

transcriptome.inv.df <- as.data.frame(as.table(transcriptome.inv))
transcriptome.inv.df <- na.omit(transcriptome.inv.df)

transcriptome.inv.df <- transcriptome.inv.df %>% 
  mutate(plastid_haplotype_1 = case_when(
    Var1 %in% cp1 ~ "1",
    Var1 %in% cp2 ~ "2",
    Var1 %in% cp3 ~ "3",
    TRUE            ~ "Not Found" 
  )) %>% 
  mutate(plastid_haplotype_2 = case_when (
    Var2 %in% cp1 ~ "1", 
    Var2 %in% cp2 ~ "2",
    Var2 %in% cp3 ~ "3",
    TRUE            ~ "Not Found" 
  )) %>% 
  mutate(haplotype_div = case_when (
    plastid_haplotype_1 == plastid_haplotype_2 ~ "same", 
    plastid_haplotype_1 != plastid_haplotype_2 ~ "different"
  )) %>% 
  filter(Freq != 0) %>% 
  group_by(haplotype_div)

transcriptome.cphaps.plot <- ggplot(transcriptome.inv.df, aes(x = haplotype_div, y = Freq, fill = haplotype_div)) +
  geom_boxplot(alpha = 0.5) +
  xlab("Same/Different Plastid Haplotype") + 
  ylab("Genetic Distance") + 
  ggtitle("Transcriptome") +
  theme_bw()

wilcox.test(Freq ~ haplotype_div, data = transcriptome.inv.df)
