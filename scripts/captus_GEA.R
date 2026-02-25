###############################################################################

# Finding candidate alleles/loci associated with environment (GEA) from CAPTUS
# Code developed by Jessie Pelosi, Abby Pearse
# Last updated: 02/23/2026

###############################################################################

### Libraries

library(vcfR)
library(ggplot2)
library(vegan)
library(adegenet)
library(reshape2)

### Reading in files

envs <- read.delim("data/LYJA_bioclim1-2-13-14-15")
popmap <- read.delim("data/popmap_states.txt", sep="", header=FALSE)
LYJA.thinned.vcf <- read.vcfR("data/captus_0.5mising.maf0.5.thinned.vcf")
si <- read.csv("data/si_table.csv")

### Constructing genind oject

LYJA.gid <- vcfR2genind(LYJA.thinned.vcf, ploidy = 4)

# Here, we are subsetting our genind to only include invaded samples
samples.native <- c("SRR29127787", "SRR29127784", "SRR29127767","SRR29127785","SRR29127770","SRR29127789","SRR29127769","SRR29127786","SRR29127774","SRR29127777","SRR29127780",
                 "SRR29127782","SRR29127788","SRR29127783","SRR29127776","SRR29127775","SRR29127772","SRR29127771","SRR29127764","SRR29127766",
                 "SRR29127773","SRR29127781","SRR29127778","SRR29127765", "SRR29127793")

LYJA.gid.inv <- LYJA.gid[!indNames(LYJA.gid) %in% samples.native]
LYJA.gid.impute <- tab(LYJA.gid.inv, freq=TRUE, NA.method="mean") # matrix containing genotypic data


###################################################### PCA ################################################

LYJA.gid.impute.PCA <- dudi.pca(LYJA.gid.impute, center=TRUE, scale=FALSE, nf=NA)
LYJA.gid.impute.PCA.df <- LYJA.gid.impute.PCA$li
LYJA.gid.impute.PCA.df$sample <- row.names(LYJA.gid.impute.PCA.df)
colnames(popmap) <- c("sample", "location")
LYJA.gid.impute.PCA.df.join <- inner_join(LYJA.gid.impute.PCA.df, popmap)

LYJA.PCA.plot <- ggplot(data = LYJA.gid.impute.PCA.df.join) + geom_point(mapping = aes(x = Axis1, y = Axis2, color = location)) +
  theme_bw() 


#################################################### RDA ############################################

LYJA.rda <- rda(LYJA.gid.impute ~., data = envs, scale = T) # Calculating redundancy analysis of genotypic data against bioclim variables
summary(eigenvals(LYJA.rda, model = "constrained")) # Here, we are summarizing the eigenvalues from our redundancy analysis
RsquareAdj(LYJA.rda)

screeplot(LYJA.rda) # We see that the majority of our variation is seen in RDA1, RDA2, and RDA3

envs$sample <- rownames(envs)
envs_add <- inner_join(envs, LYJA.gid.impute.PCA.df.join)

# Here, we perform a redundancy analysis against bioclim variables and variance explained by PCs
LYJA.pRDA <- rda(LYJA.gid.impute ~ + bio12 + bio14 + bio2 + bio15 + bio13 + bio1 + Condition(Axis1+Axis2), envs_add) 
summary(eigenvals(LYJA.pRDA, model = "constrained"))
screeplot(LYJA.pRDA) # We see that majority of variation is seen in RDA1, RDA2, and RDA3

# Here, we look at directionality of SNPs explained by bioclim variables according to RDA1 and RDA2
plot(LYJA.pRDA, type = 'n', scaling = 3)
points(LYJA.pRDA, display = "sites", pch =21, cex=1.75, col = 'black', scaling = 3)
text(LYJA.pRDA, scaling = 3, display = "bp", col = "black", cex=1.2)

# Here, we look at directionality of SNPs explained by bioclim variables according to RDA1 and RDA3
plot(LYJA.pRDA, type = 'n', scaling = 3, choices =c(1,3))
points(LYJA.pRDA, diplay = "sites", pch=21, cex=1.75, col="black", scaling = 3, choices = c(1,3))
text(LYJA.pRDA, scaling = 3, display = "bp", col="black", cex=1.2, choices = c(1,3))



################################## Calculating SNPs explained by each environmental predictor #######################

load.rda <- scores(LYJA.pRDA, choices=c(1:3), display = "species")
load.rda.unique <- unique(load.rda)

# Here, we are obtaining the collective top 5% of SNPs explained by environmental associated. Numeric thresholds are calculated by scaling RDA 
#eigenvalue to summed eigvenvals of top 3 RDA, then multiplied by 222 (5% of total SNPs), and finally divided by 4439 (total number of SNPs)

cand1 <- which(abs(load.rda.unique[,1]) >= quantile(abs(load.rda.unique[,1]), 0.983))
cand2 <- which(abs(load.rda.unique[,2]) >= quantile(abs(load.rda.unique[,2]), 0.984))
cand3 <- which(abs(load.rda.unique[,3]) >= quantile(abs(load.rda.unique[,3]), 0.984))

ncand <- length(cand1) + length(cand2) +length(cand3)

# Here, we create a dataframe of the loading value of each SNP
cand1 <- cbind.data.frame(rep(1, times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(1, times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(1, times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis", "snp", "loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

foo <- matrix(nrow=(ncand), ncol=6)

colnames(foo) <- c("bio12", "bio14", "bio2", "bio15", "bio13", "bio1")

predictors <- envs[,1:6]

for(i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- LYJA.gid.impute[,nam]
  foo[i,] <- apply(predictors,2,function(x) cor(x, snp.gen))
}

cand <- cbind.data.frame(cand,foo)

for(i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,10] <- names(which.max(abs(bar[4:9])))
  cand[i,11] <- max(abs(bar[4:9]))
}

colnames(cand)[10] <- "predictor"
colnames(cand)[11] <- "correlation"

envs_predictors <- levels(as.factor(cand$predictor))
no_snps <- as.numeric(table(cand$predictor))
piedata <- data.frame(envs_predictors, no_snps)


ggplot(piedata, aes(x="", y=no_snps, fill=envs_predictors)) +
  geom_col(color = "black") +
  geom_text(aes(label = no_snps), 
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta="y") +
  ggtitle("Nunmber of Candidate SNPs Per Environmental Predictor") +
  theme_void()

sel <- cand$snp
env <- cand$predictor

env[env=="bio12"] <- "blue"
env[env=="bio14"] <- "darkgreen"
env[env=="bio2"] <- "lightpink"
env[env=="bio15"] <- "goldenrod"
env[env=="bio13"] <- "lightblue"
env[env=="bio1"] <- "lavender"


col.pred <- rownames(LYJA.pRDA$CCA$v)

for(i in 1:length(sel)){
  foo <- match(sel[i], col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("L", col.pred)] <- "grey"
empty <- col.pred

empty[grep("gray", empty)] <- rgb(0,1,0, alpha =0)

empty.outline <- ifelse(empty =="gray", "gray", "gray32")
bg <- c("blue", "darkgreen", "lightpink", "goldenrod", "lightblue", "lavender")

plot(LYJA.pRDA, type = "n", scaling = 3, xlim = c(-0.75, 0.75), ylim=c(-1,1))
points(LYJA.pRDA, display = "species", pch=21, cex = 1.3, col = "gray32", bg=col.pred, scaling =3)
points(LYJA.pRDA, display = "species", pch=21, cex=1.3, col=empty.outline, bg=empty, scaling = 3)
text(LYJA.pRDA, scaling = 3, display = "bp", col = "black", cex = 1.2)

############################################## Calculating heterozygosity of candidate loci ##############################################


cand$snp_cleaned <- gsub("\\.[0-9]", "", cand$snp) 

GEA_loci <- LYJA.gid[loc = cand$snp_cleaned]

gl.GEA <- gi2gl(GEA_loci)
ploidy(gl.GEA) <- 4
gl.compliance.check(gl.GEA)

gl.report.polyploid_heterozygosity(gl.GEA, error.bar = "SE") 
Ho.GEA <- gl.report.polyploid_heterozygosity(gl.GEA, method = "ind")

colnames(Ho.GEA) <- c("SRR", "Ho", "f.hom.ref", "f.hom.alt", "n.Loc")
Ho.GEA.popmap.year <- inner_join(Ho.GEA, si)


ggplot(data = Ho.GEA.popmap.year, mapping = aes(x = Collection.Year, y = Ho, color = Invaded)) + geom_point() + 
  geom_smooth(method = "lm") + theme_bw() +
  stat_regline_equation(aes(label =  paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")))

### Running statistical models

Ho.GEA.lm <- lm(Ho ~ Collection.Year + Invaded + Collection.Year:Invaded, data = Ho.GEA.popmap.year)
summary(Ho.GEA.lm) 

Ho.GEA.noint.lm <- lm(Ho ~ Collection.Year + Invaded, data = Ho.GEA.popmap.year)
summary(Ho.GEA.noint.lm) 

