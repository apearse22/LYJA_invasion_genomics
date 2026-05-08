###############################################################################

# Finding candidate alleles/loci associated with environment (GEA) from CAPTUS
# Code developed by Jessie Pelosi, Abby Pearse
# Last updated: 04/27/2026

###############################################################################

### Libraries

library(vcfR)
library(ggplot2)
library(vegan)
library(adegenet)
library(reshape2)
library(gridExtra)
library(dartR)
library(dartR.base)
library(dplyr)
library(adegenet)
library(broom)

### Reading in / manipulating files

captus <- read.vcfR("files/captus.SNPs.0.5missing.CTmarked.recalc.maf0.5.thinned.vcf")
envs <- read.csv("files/sample_envs.csv")
popmap <- read.csv("files/popmap_updatedcoords.csv", sep=",")
si <- read.csv("files/si_table.csv")


### Shared functions and variables

outliers <- function(x,z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}

bg=c("#ff7f00","#1f78b4","#ffff33","#a6cee3","#33a02c","#e31a1c", "purple", "brown4")

################################### Creating geninds  #######################

captus.genind <- vcfR2genind(captus, ploidy = 4)
pop(captus.genind) <- popmap$Temporal_Group

### subsetting captus genind to only include invaded samples

samples.native <- c("SRR29127787", "SRR29127784", "SRR29127767","SRR29127785","SRR29127770","SRR29127789","SRR29127769","SRR29127786","SRR29127774","SRR29127777","SRR29127780",
                 "SRR29127782","SRR29127788","SRR29127783","SRR29127776","SRR29127775","SRR29127772","SRR29127771","SRR29127764","SRR29127766",
                 "SRR29127773","SRR29127781","SRR29127778","SRR29127765", "SRR29127793")


captus.invaded <- captus.genind[!indNames(captus.genind) %in% samples.native] 
captus.impute <- tab(captus.invaded, freq = TRUE, NA.method = "mean") # matrix containing genotypic data


envs <- envs[!envs$Ind %in% samples.native,] # subset environmental data to only be invaded range
envs2 <- dplyr::select(envs, BIO1, BIO2, BIO13, BIO14, BIO15, Ind)
rownames(envs2) <- envs$Ind
colnames(envs2) <- c("BIO1", "BIO2", "BIO13", "BIO14", "BIO15", "Ind")


predictors <- envs2[,1:5]

envs3 <- inner_join(popmap, envs2)
levels(envs3$Location) <- c("Florida", "Texas", "Alabama", "Arkansas", "South_Carolina", "Georgia", "Mississippi", "Louisiana")
state <- envs3$Location



###################################################### PCA #################################################

impute.PCA <- dudi.pca(captus.impute, center=TRUE, scale=FALSE, nf=50, scannf = F)
impute.PCA.df <- impute.PCA$li
impute.PCA.df$Ind <- row.names(impute.PCA.df)
impute.PCA.df.popmap <- inner_join(impute.PCA.df, popmap)

PCA.plot <- ggplot(impute.PCA.df.popmap, aes(Axis1, Axis2, color = Location)) +
  geom_point() +
  theme_bw()


##################################################### RDA ################## - ask jessie about this part

rda <- rda(captus.impute ~., envs2, scale = T) # Calculating redundancy analysis of genotypic data against bioclim variables
summary(eigenvals(rda, model = "constrained")) # summarizing the eigenvalues from our redundancy analysis
RsquareAdj(rda) # calculating R squared values
screeplot(rda) # majority of our variation is seen in RDA1, RDA2, and RDA3


################################################ pRDA #######################################


#envs2$sample <- rownames(envs2)

envs.add <- inner_join(envs2, impute.PCA.df.popmap)

pRDA <- rda(captus.impute ~ + BIO14 + BIO2 + BIO15 + BIO13 + BIO1 + Condition(Axis1+Axis2), envs.add) 
summary(eigenvals(pRDA, model = "constrained"))
screeplot(pRDA) # majority of variation is seen in RDA1, RDA2, and RDA3

plot(pRDA, scaling = 3) # plotting pRDA results, blue vectors are the environmental predictors 


############################## plotting pRDA results based on environmental predictors #######################

points(pRDA, display = "species", pch = 20, cex=0.7, col="gray32", scaling = 3) #SNPs
points(pRDA, display = "sites", pch = 21, cex = 1.3, col = "gray32", bg=bg, scaling = 3)
text(pRDA, scaling = 3, display = "bp", col = "#0868ac", cex = 1)

#legend("bottomright", legend = levels(state), bty= 'n', col = "gray32", pch = 21, pt.bg = bg)


########## plotting directionality of SNPs explained by bioclim variables, according to RDA1/RDA2 & RDA1 / RDA3 ########

## RDA 1 & RDA 2
plot(pRDA, type = 'n', scaling = 3)
points(pRDA, display = "sites", pch =21, cex=1.75, col = 'black', scaling = 3)
text(pRDA, scaling = 3, display = "bp", col = "black", cex=1.2)

## RDA 2 & RDA 3
plot(pRDA, type = 'n', scaling = 3, choices =c(1,3))
points(pRDA, display = "sites", pch=21, cex=1.75, col="black", scaling = 3, choices = c(1,3))
text(pRDA, scaling = 3, display = "bp", col="black", cex=1.2, choices = c(1,3))


############################### Calculating SNPs explained by each environmental predictor ####################

load.rda <- scores(pRDA, choices = c(1:3), display = "species")
load.rda.unique <- unique(load.rda)

#### exploration of distribution of unique SNPs
hist(load.rda.unique[,1])
hist(load.rda.unique[,2])
hist(load.rda.unique[,3])

#### Identifying SNPs that are outliers 
cand1 <- outliers(load.rda.unique[,1], 3) # 65
cand2 <- outliers(load.rda.unique[,2], 3) # 56
cand3 <- outliers(load.rda.unique[,3], 3) # 71

ncand <- length(cand1) + length(cand2) + length(cand3) 

### Creating a dataframe of the loading value of each SNP
cand1.df <- cbind.data.frame(rep(1, times=length(cand1)), names(cand1), unname(cand1))
cand2.df <- cbind.data.frame(rep(1, times=length(cand2)), names(cand2), unname(cand2))
cand3.df <- cbind.data.frame(rep(1, times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1.df) <- c("axis", "snp", "loading")
colnames(cand2.df) <- c("axis", "snp", "loading")
colnames(cand3.df) <- c("axis", "snp", "loading")


cand.total.df <- rbind(cand1.df, cand2.df, cand3.df)
cand.total.df$snp <- as.character(cand.total.df$snp)

foo <- matrix(nrow=(ncand), ncol=5)
colnames(foo) <- c("bio14", "bio2", "bio15", "bio13", "bio1")

for(i in 1:length(cand.total.df$snp)) {
  nam <- cand.total.df[i,2]
  snp.gen <- captus.impute[,nam]
  foo[i,] <- apply(predictors,2,function(x) cor(x, snp.gen))
}

cand.total.df.foo <- cbind.data.frame(cand.total.df,foo)

for(i in 1:length(cand.total.df.foo$snp)) {
  bar <- cand.total.df.foo[i,]
  cand.total.df.foo[i,9] <- names(which.max(abs(bar[4:8])))
  cand.total.df.foo[i,10] <- max(abs(bar[4:8]))
}

colnames(cand.total.df.foo)[9] <- "predictor"
colnames(cand.total.df.foo)[10] <- "correlation"



########################################### separating number of loci in coding vs noncoding ##############################

txm.captus.blast <- read.delim("files/LYJA.CAPTUS.blast.out")
colnames(txm.captus.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                "qstart",
                                "qend", "sstart", "send", "evalue", "bitscore")

txm.captus.blast.05eval <- txm.captus.blast %>% 
  filter(evalue <= 1e-05) # forget to set evalue threshold when creating file in HPC

best.hit.txm.captus.blast.05eval <- txm.captus.blast.05eval %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(row_number() == 1)



txm.arabidopsis.blast <- read.delim("files/Arabidopsis.LYJA.blast")
colnames(txm.arabidopsis.blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                                               "qend", "sstart", "send", "evalue", "bitscore")

best.hit.txm.arabidopsis.blast <- txm.arabidopsis.blast %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(row_number() == 1)

# at some point, get best hit per locus 


################################# creating / subsetting genind objects ########################

# subsetting captus genind based on hits to transcriptome

#locNames(captus.genind) <- gsub("_[0-9]..", "", locNames(captus.genind))
#locNames(captus.genind) <- gsub("_[0-9].", "", locNames(captus.genind))


locNames(captus.invaded) <- gsub("_[0-9\\.].*", "", locNames(captus.invaded))


captus.coding.genind <- captus.invaded[loc = locNames(captus.invaded) %in% best.hit.txm.captus.blast.05eval$qseqid.cleaned]
captus.noncoding.genind <- captus.invaded[loc = !locNames(captus.invaded) %in% best.hit.txm.captus.blast.05eval$qseqid.cleaned]


cand.total.df.foo$snp.cleaned <- gsub("_[0-9].*", "", cand.total.df.foo$snp)
cand.total.df.foo.unique <- cand.total.df.foo %>% 
  distinct(cand.total.df.foo$snp.cleaned, .keep_all = TRUE)


gea.loci.in.coding <- captus.coding.genind[loc = locNames(captus.coding.genind) %in% cand.total.df.foo.unique$snp.cleaned]
gea.loci.in.noncoding <- captus.noncoding.genind[loc = locNames(captus.noncoding.genind) %in% cand.total.df.foo.unique$snp.cleaned]

envs.predictors.unique <- levels(as.factor(cand.total.df.foo.unique$predictor))
no.snps.unique <- as.numeric(table(cand.total.df.foo.unique$predictor))
piedata.unique <- data.frame(envs.predictors.unique, no.snps.unique)


ggplot(piedata.unique, aes(x="", y=no.snps.unique, fill=envs.predictors.unique)) +
  geom_col(color = "black") +
  geom_text(aes(label = no.snps.unique), 
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta="y") +
  ggtitle("Nunmber of Candidate SNPs Per Environmental Predictor") +
  theme_void()

# descending bar plot of the 114 loci

ggplot(data = piedata.unique, aes(x = reorder(envs.predictors.unique, -no.snps.unique), y = no.snps.unique, fill = envs.predictors.unique)) + 
  geom_col() + 
  xlab("Bioclimatic Variable") +
  ylab("Number of Candidate SNPs") +
  scale_fill_manual(values = c("#0074D9", "#B10DC9", "#85144b", "#FF4136", "#FF851B")) + 
  theme_bw()+ 
  theme(legend.position = "none")


#ggsave("GEA_descending_bars_updated.pdf", height = 4, width = 6)
#ggsave("GEA_descending_bars_updated.png", height = 4, width = 6, dpi = 300)


################################## doing a chi square test between coding and noncoding loci

cand.total.df.foo.unique.coding <- cand.total.df.foo.unique %>% 
  filter(snp.cleaned %in% locNames(gea.loci.in.coding)) %>% 
  dplyr::select(bio14, bio2, bio15, bio13, bio1, correlation, predictor, snp.cleaned)


envs.predictors.unique.coding <- levels(as.factor(cand.total.df.foo.unique.coding$predictor))
no.snps.unique.coding <- as.numeric(table(cand.total.df.foo.unique.coding$predictor))
piedata.unique.coding <- data.frame(envs.predictors.unique.coding, no.snps.unique.coding)
colnames(piedata.unique.coding) <- c("bioclim", "no. coding snps")



cand.total.df.foo.unique.noncoding <- cand.total.df.foo.unique %>% 
  filter(snp.cleaned %in% locNames(gea.loci.in.noncoding)) %>% 
  dplyr::select(bio14, bio2, bio15, bio13, bio1, correlation, predictor, snp.cleaned)


envs.predictors.unique.noncoding <- levels(as.factor(cand.total.df.foo.unique.noncoding$predictor))
no.snps.unique.noncoding <- as.numeric(table(cand.total.df.foo.unique.noncoding$predictor))
piedata.unique.noncoding <- data.frame(envs.predictors.unique.noncoding, no.snps.unique.noncoding)
colnames(piedata.unique.noncoding) <- c("bioclim", "no. noncoding snps")


chi.df <- inner_join(piedata.unique.coding, piedata.unique.noncoding)
row.names(chi.df) <- chi.df$bioclim
chi.df <- chi.df %>% 
  dplyr::select(-bioclim)

chi.df.results <- chisq.test(chi.df) # not significant


########################## Plotting number of SNPs explained by environmental predictor

envs.predictors <- levels(as.factor(cand.total.df.foo$predictor))
no.snps <- as.numeric(table(cand.total.df.foo$predictor))
piedata <- data.frame(envs.predictors, no.snps)

### Pie chart

ggplot(piedata, aes(x="", y=no.snps, fill=envs.predictors)) +
  geom_col(color = "black") +
  geom_text(aes(label = no.snps), 
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta="y") +
  ggtitle("Nunmber of Candidate SNPs Per Environmental Predictor") +
  theme_void()


### Descending bar plot
ggplot(data = piedata, aes(x = reorder(envs.predictors, -no.snps), y = no.snps, fill = envs.predictors)) + 
  geom_col() + 
  xlab("Bioclimatic Variable") +
  ylab("Number of Candidate SNPs") +
  scale_fill_manual(values = c("#0074D9", "#B10DC9", "#85144b", "#FF4136", "#FF851B")) + 
  theme_bw()+ 
  theme(legend.position = "none")

ggsave("GEA_descending_bars.pdf", height = 4, width = 6)
ggsave("GEA_descending_bars.png", height = 4, width = 6, dpi = 300)


###################################### finding functions of arabidopsis genes in coding environmental loci

coding.gea.subsetted <- cand.total.df.foo.unique.coding %>% 
  dplyr::select(predictor, snp.cleaned)
colnames(coding.gea.subsetted) <- c("predictor", "qseqid")

coding.gea.in.txm.captus.blast <- inner_join(coding.gea.subsetted, best.hit.txm.captus.blast.05eval)
coding.gea.in.txm.captus.blast.subsetted <- coding.gea.in.txm.captus.blast %>% 
  dplyr::select(predictor, qseqid, sseqid, evalue)
colnames(coding.gea.in.txm.captus.blast.subsetted) <- c("predictor", "qseqid", "ID", "evalue")


best.hit.txm.arabidopsis.blast.subset <- best.hit.txm.arabidopsis.blast %>% 
  dplyr::select(qseqid, sseqid)

colnames(best.hit.txm.arabidopsis.blast.subset) <- c("ID", "AT.gene")

coding.gea.loci.in.arabidopsis <- inner_join(coding.gea.in.txm.captus.blast.subsetted, best.hit.txm.arabidopsis.blast.subset)



############################################## Calculating heterozygosity of candidate loci ##############################################


cand.total.df.foo$snp.cleaned <- gsub("\\.[0-9]", "", cand.total.df.foo$snp) 

GEA.loci <- captus.genind[loc = cand.total.df.foo$snp.cleaned]

gl.GEA <- gi2gl(GEA.loci)
ploidy(gl.GEA) <- 4
gl.compliance.check(gl.GEA)

gl.report.polyploid_heterozygosity(gl.GEA, error.bar = "SE") 
Ho.GEA <- gl.report.polyploid_heterozygosity(gl.GEA, method = "ind")

colnames(Ho.GEA) <- c("SRR", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")
Ho.GEA.popmap.year <- inner_join(Ho.GEA, si)


ggplot(data = Ho.GEA.popmap.year, mapping = aes(x = Collection.Year, y = Ho, color = Invaded, shape = Invaded)) + geom_point(size = 2) + 
  geom_smooth(method = "lm") + theme_bw() +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  scale_shape_manual(values = c(17, 16)) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  ggtitle("CAPTUS Reference - GEA Loci") +
  theme(text = element_text(size = 16))



### Running statistical models

Ho.GEA.popmap.year$std.Year <- Ho.GEA.popmap.year$Collection.Year - min(Ho.GEA.popmap.year$Collection.Year)

Ho.GEA.lm <- lm(Ho ~ std.Year + Invaded + std.Year:Invaded, data = Ho.GEA.popmap.year)
summary(Ho.GEA.lm) 

Ho.GEA.noint.lm <- lm(Ho ~ std.Year + Invaded, data = Ho.GEA.popmap.year)
summary(Ho.GEA.noint.lm)


r2_valu <- summary(Ho.GEA.noint.lm)$r.squared
f <- summary(Ho.GEA.noint.lm)$fstatistic

p_valu <- pf(f[1], f[2], f[3], lower.tail = F)

stats_label <- paste0("Adj.R^2 == ", round(r2_valu, 5),
                      "~~italic(P) == ", round(p_valu, 3))

ggplot(data = Ho.GEA.popmap.year, mapping = aes(x = Collection.Year, y = Ho, color = Invaded, shape = Invaded)) + geom_point(size = 2) + 
  geom_smooth(method = "lm", linetype = "dashed") + theme_bw() +
  #stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) +
  annotate("text", 
           x = 1945, y = 0.5, 
           label = stats_label, 
           parse = T, 
           hjust = 1.1, vjust=1.5, size = 5) +
  scale_shape_manual(values = c(17, 16)) +
  xlab("Collection Year") +
  scale_color_manual(values = c("Y" = "#21918c", "N" = "#440154")) +
  theme(text = element_text(size = 16), legend.position = "none") +
  ggtitle("C) GEA Candidate Loci")


ggsave("captus.GEA.ho.ind.popmap.plot.pdf", width = 8, height =6)
ggsave("captus.GEA.ho.ind.popmap.plot.png", width = 8, height =6, dpi = 300)








# Let's look at fixation now
#dosage_matrix <- as.matrix(gl.GEA) 
#dosage_matrix <- as.data.frame(dosage_matrix / 4)
#dosage_matrix$Ind <- rownames(dosage_matrix)

#popmap2 <- read.csv("popmap.csv")

#dosage_locale <- inner_join(dosage_matrix, popmap2)
#dosage_melt_locale <- melt(dosage_locale, id.vars = c("Ind", "Herbarium", "Voucher.Identifier", 
#                                                      "Location", "Approx..Latitude","Approx..Longitude", "Invaded", 
#                                                      "Temporal_Group", "Specimen.Comments", "Specimen.Link", "Collection.Year")) %>% 
#  dplyr::select(Ind, Invaded, Collection.Year, variable, value) %>% 
#  filter(variable != "total...bases")

#ggplot(data = dosage_melt_locale, mapping = aes(x = Collection.Year, y =value, color = Invaded)) + 
#  geom_point() + facet_wrap(~variable)


# Is there a significant trend toward fixation? 
#trend_results <- dosage_melt_locale %>%
#  filter(Invaded == "Y") %>% 
#  group_by(variable) %>%
#  summarize(
#    correlation = cor(Collection.Year, value, method = "spearman", use = "complete.obs"),
#    p_value = cor.test(Collection.Year, value, method = "spearman")$p.value
#  ) %>%
#  filter(!is.na(correlation)) %>%
#  arrange(desc(abs(correlation)))


# Pull out the loci with significant trends: 
#trend_results_sig <- filter(trend_results, p_value < 0.05) 

#allele_freq_changes_sig <- dosage_melt_locale[dosage_melt_locale$variable %in% trend_results_sig$variable, ]
#allele_freq_changes_sig <- na.omit(allele_freq_changes_sig)

#ggplot(data = allele_freq_changes_sig, mapping = aes(x = Collection.Year, y = value, color = Invaded)) + 
#  geom_point() + geom_smooth(se = F, method = "lm") + facet_wrap(~variable) + theme_bw() + 
#  ylim(0,1) + xlab("Collection Year") + ylab("Allele Frequency")



# Let's do this now with the temporal groupings:

#pop(LYJA.gid) <- popmap2$Temporal_Group

#GEA_loci <- LYJA.gid[loc = cand$snp_cleaned]
#ploidy(gl.GEA) <- 4
#gl.GEA <- gi2gl(GEA_loci)
#gl.compliance.check(gl.GEA)

#gl.report.polyploid_heterozygosity(gl.GEA, error.bar = "SE") 
#Ho.GEA <- gl.report.polyploid_heterozygosity(gl.GEA, method = "pop")


########### old analyses 


#sel <- cand$snp
#env <- cand$predictor

#env[env=="bio12"] <- "blue"
#env[env=="bio14"] <- "darkgreen"
#env[env=="bio2"] <- "lightpink"
#env[env=="bio15"] <- "goldenrod"
#env[env=="bio13"] <- "lightblue"
#env[env=="bio1"] <- "lavender"

#col.pred <- rownames(LYJA.pRDA$CCA$v)

#for(i in 1:length(sel)){
#  foo <- match(sel[i], col.pred)
#  col.pred[foo] <- env[i]
#}

#col.pred[grep("L", col.pred)] <- "grey"
#empty <- col.pred

#empty[grep("gray", empty)] <- rgb(0,1,0, alpha =0)

#empty.outline <- ifelse(empty =="gray", "gray", "gray32")
#bg <- c("blue", "darkgreen", "lightpink", "goldenrod", "lightblue", "lavender")

#plot(LYJA.pRDA, type = "n", scaling = 3, xlim = c(-0.75, 0.75), ylim=c(-1,1))
#points(LYJA.pRDA, display = "species", pch=21, cex = 1.3, col = "gray32", bg=col.pred, scaling =3)
#points(LYJA.pRDA, display = "species", pch=21, cex=1.3, col=empty.outline, bg=empty, scaling = 3)
#text(LYJA.pRDA, scaling = 3, display = "bp", col = "black", cex = 1.2)
