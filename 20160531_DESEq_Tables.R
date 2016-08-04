#LOADING R PACKAGES#
pac <- c("survival", "phyloseq", "ade4", "vegan", "splines", "Kendall", "lme4", "languageR", 
         "doBy", "MASS", "plyr", "reshape2", "DESeq2")

lapply(pac, library, character.only = T)

#SET UP WORKING DIRECTORY#
setwd("/Users/hayesr03/Documents/Gitstuff/MB_HNC")

#LOADING MAPPING FILE#
map <- read.csv("all0405.csv")

map$age <- as.numeric(as.character(map$AGE))
map$gender[map$sex == "Male"] <- "1"
map$gender[map$sex == "Female"] <- "2"
map$race2 <- as.factor(as.character(map$race2))
map$gender <- as.factor(as.character(map$gender))
map$studytype <- as.factor(as.character(map$studytype))
map$Case[map$case == 1] <- 1
map$Case[map$case == 0] <- 0
map$Case <- as.factor(as.character(map$Case))

#ALCOHOL AND SMOKING STATUS#
map$alc_yn <-0
map$alc_yn[map$ETHA_GRAMS_PER_DAY > 0] <-1
map$alc_yn[is.na(map$ETHA_GRAMS_PER_DAY)] <- 2
map$alc_yn <- as.factor(as.character(map$alc_yn))
map$ETHA_GRAMS_PER_DAY[is.na(map$ETHA_GRAMS_PER_DAY)] <- 0
map$cig_per_day[map$cig_per_day == 999] <- 0
map$cig_per_day[map$cig_per_day == 9999] <- 15
map$matchset <- as.factor(as.character(map$matchset))
map$CIG_STAT <- as.factor(as.character(map$CIG_STAT))
map$smk_yn <- 1
map$smk_yn[map$CIG_STAT == 0] <- 0

# table(map$ETHA_GRAMS_PER_DAY, map$alc_yn)
map <- subset(map, map$Case == 0 | map$Case == 1)  #Call Cases
map.sq <- subset(map, map$Case == 0 | map$sq_cell == 1)  #Sq cell only
row.names(map) <- as.character(map$SampleID)
row.names(map.sq) <- as.character(map.sq$SampleID)

save(map, file = "map.Rdata")
save(map.sq, file = "map.sq.Rdata")


##TABLE 2 DESEQ2##
#ALL CaseS##
load("phy.phylum.abd.href.Rdata")
row.names(map) <- as.character(map$SampleID)
sample_data(phy.abd.phylum) <- map

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 100,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#mean counts#
#Fusobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Proteobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo30924 ~ Case, data=map.all, FUN=c(mean, sd))
#Bacteroidetes#
summaryBy(map.all$denovo82916 ~ Case, data=map.all, FUN=c(mean, sd))
#Firmicutes#
summaryBy(map.all$denovo8807 ~ Case, data=map.all, FUN=c(mean, sd))


#SQ CaseS##
load("phy.phylum.abd.href.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
sample_data(phy.abd.phylum) <- map.sq

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 100,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#mean counts#
#Fusobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Proteobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo30924 ~ Case, data=map.all, FUN=c(mean, sd))
#Bacteroidetes#
summaryBy(map.all$denovo82916 ~ Case, data=map.all, FUN=c(mean, sd))
#Firmicutes#
summaryBy(map.all$denovo8807 ~ Case, data=map.all, FUN=c(mean, sd))


###TABLE 3###
#TEST FOR INTERACTION#
#SMOKING#
design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case, data = map)
dds = DESeq(phylum.deseq, test = "LRT",
            full =~age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case + CIG_STAT*Case,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case )

res <- data.frame(results(dds))
res <-merge(res, taxa, by =0)
res$padj <- p.adjust(res$pvalue, method="fdr")
res

#ALCOHOL#
design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
                                ETHA_GRAMS_PER_DAY + Case + Case*alc_yn, data = map)
dds = DESeq(phylum.deseq, test = "LRT",
            full =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case + Case*alc_yn,
            reduced =~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + 
              ETHA_GRAMS_PER_DAY + Case)

res <- data.frame(results(dds))
res$padj <- p.adjust(res$pvalue, method="fdr")
res<-merge(res, taxa, by =0)
res


##EVER SMOKERS##
map.smk <- subset(map.sq, map.sq$CIG_STAT %in% c("1", "2"))
map.nsmk <- subset(map.sq, map.sq$CIG_STAT == "0")

table(map.smk$Case)
table(map.nsmk$Case)

load("phy.phylum.abd.href.Rdata")
row.names(map.smk) <- as.character(map.smk$SampleID)
sample_data(phy.abd.phylum) <- map.smk

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 100,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map.smk)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)
res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))



#Never-smokers#
load("phy.phylum.abd.href.Rdata")
row.names(map.nsmk) <- as.character(map.nsmk$SampleID)
sample_data(phy.abd.phylum) <- map.nsmk

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 100,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map.nsmk)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))




##EVER ALCOHOL##
map.alc <- subset(map.sq, map.sq$alc_yn == "1")
map.nalc <- subset(map.sq, map.sq$alc_yn == "0")

table(map.alc$Case)
table(map.nalc$Case)

load("phy.phylum.abd.href.Rdata")
row.names(map.alc) <- as.character(map.alc$SampleID)
sample_data(phy.abd.phylum) <- map.alc

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 100,]
phylum.deseq


design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + ETHA_GRAMS_PER_DAY + Case, 
                              data = map.alc)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)
res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))



#Never drinkers#
load("phy.phylum.abd.href.Rdata")
row.names(map.nalc) <- as.character(map.nalc$SampleID)
sample_data(phy.abd.phylum) <- map.nalc

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 1,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + Case, data = map.nalc)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))


##STABLE2##
##CPSII#
load("phy.phylum.abd.href.Rdata")
row.names(map) <- as.character(map$SampleID)

map.a <- subset(map, map$studytype == 1)
sample_data(phy.abd.phylum) <- map.a
table(map.a$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 10,]
phylum.deseq

design(phylum.deseq)<-formula(~ age +  gender + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map.a)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p =round(subs$pvalue2,digits = 4),  p =round(subs$pvalue,digits = 4))
res

count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))
#Fusobacteria#
summaryBy(map.all$denovo30924 ~ Case, data=map.all, FUN=c(mean, sd))
#Bacteroidetes#
summaryBy(map.all$denovo82916 ~ Case, data=map.all, FUN=c(mean, sd))
#Firmicutes#
summaryBy(map.all$denovo8807 ~ Case, data=map.all, FUN=c(mean, sd))


##PLCO##
load("phy.phylum.abd.href.Rdata")
row.names(map) <- as.character(map$SampleID)
map.b <- subset(map, map$studytype == 2)

sample_data(phy.abd.phylum) <- map.b
table(map.b$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 10,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map.b)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             p =round(subs$pvalue2,digits = 4),  p =round(subs$pvalue,digits = 4))
res


count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))
#Fusobacteria#
summaryBy(map.all$denovo30924 ~ Case, data=map.all, FUN=c(mean, sd))
#Bacteroidetes#
summaryBy(map.all$denovo82916 ~ Case, data=map.all, FUN=c(mean, sd))
#Firmicutes#
summaryBy(map.all$denovo8807 ~ Case, data=map.all, FUN=c(mean, sd))


##Supplementary Table 3 SITE##
#SQ CaseS##
#Oral Cavity#
load("phy.phylum.abd.href.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
map.s1 <- subset(map.sq, map.sq$site == 1 | map.sq$Case == 0)

sample_data(phy.abd.phylum) <- map.s1
table(map.s1$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 1,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res


count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.s1, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteri#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteriaa#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))


#Pharynx
load("phy.phylum.abd.href.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
map.s2 <- subset(map.sq, map.sq$site == 2 | map.sq$Case == 0)

sample_data(phy.abd.phylum) <- map.s2
table(map.s2$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 1,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res


count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.s2, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))



#Larynx
load("phy.phylum.abd.href.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
map.s3 <- subset(map.sq, map.sq$site == 3 | map.sq$Case == 0)

sample_data(phy.abd.phylum) <- map.s3
table(map.s3$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 1,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res


count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.s3, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))



#SUPPLEMENTARY TABLE 4 BY TIMETODX#
#year 0-2
load("phy.phylum.abd.href.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
map.sq$TIMETODX <- as.numeric(as.character(map.sq$TIMETODX))
map.t1 <- subset(map.sq, map.sq$TIMETODX <=2 | map.sq$Case == 0)

sample_data(phy.abd.phylum) <- map.t1
table(map.t1$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 1,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res


count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.t1, count_phy, by.x = "SampleID", by.y = "row.names")

#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))



#year 3-10#
load("phy.phylum.abd.href.Rdata")
load("map.sq.Rdata")
row.names(map.sq) <- as.character(map.sq$SampleID)
map.sq$TIMETODX <- as.numeric(as.character(map.sq$TIMETODX))
map.t2 <- subset(map.sq, map.sq$TIMETODX >2 | map.sq$Case == 0)

sample_data(phy.abd.phylum) <- map.t2
table(map.t2$Case)

phy.subs <- subset_taxa(phy.abd.phylum, Phylum %in% c("p__Actinobacteria","p__Bacteroidetes",
                                                      "p__Firmicutes",  "p__Fusobacteria", "p__Proteobacteria"))
taxa <- data.frame(tax_table(phy.abd.phylum))

phylum.deseq=phyloseq_to_deseq2(phy.subs, ~1)
phylum.deseq <- phylum.deseq[rowSums( counts(phylum.deseq) > 2 ) >= 1,]
phylum.deseq

design(phylum.deseq)<-formula(~ age + race2 + gender + studytype + CIG_STAT + cig_per_day + alc_yn + ETHA_GRAMS_PER_DAY + Case, data = map)
dds_phy2<-DESeq(phylum.deseq, minReplicatesForReplace = Inf)

resultsNames(dds_phy2)

res_phy2=results(dds_phy2, contrast = c("Case", "1","0"), cooksCutoff = FALSE,independentFiltering = T)

mcols(dds_phy2)$maxCooks2 <- apply(assays(dds_phy2)[["cooks"]], 1, max)
identical(rownames(dds_phy2),rownames(res_phy2))
res_phy2$maxcook2=mcols(dds_phy2)$maxCooks2
res_phy2$pvalue2=res_phy2$pvalue
res_phy2$pvalue2[res_phy2$maxcook2 > 20] = NA
res_phy2$padj2 <- p.adjust(res_phy2$pvalue2, method="fdr")
res_phy2$lci <-res_phy2$log2FoldChange - 1.96*abs(res_phy2$lfcSE)
res_phy2$uci <-res_phy2$log2FoldChange + 1.96*abs(res_phy2$lfcSE)

res_phy2=data.frame(res_phy2)
res_phy2<-merge(res_phy2, taxa, by =0)
subs <- res_phy2

subs$Phylum <- as.character(subs$Phylum)
subs$Class <- as.character(subs$Class)
subs$Order <- as.character(subs$Order)
subs$Family <- as.character(subs$Family)
subs$Genus <- as.character(subs$Genus)
subs$Species <- as.character(subs$Species)

res <- cbind(taxa = subs$Row.names, phy =subs$Phylum, cla = subs$Class, Ord = subs$Order, fam = subs$Family,
             Est = round(2^(subs$log2FoldChange), digits = 2), 
             paste("(", LL = round(2^(subs$log2FoldChange- 1.96 * subs$lfcSE), digits = 2), ",", " ",  
                   UL = round(2^(subs$log2FoldChange + 1.96 * subs$lfcSE), digits = 2), ")"), 
             qadj  = round(subs$padj2,digits = 4), q  = round(subs$padj,digits = 6))
res


count_phy <- data.frame(t(counts(dds_phy2, normalized = T)), check.names = F)
map.all <- merge(map.t2, count_phy, by.x = "SampleID", by.y = "row.names")
#Proteobacteria#
summaryBy(map.all$denovo39370 ~ Case, data=map.all, FUN=c(mean, sd))
#Actinobacteria#
summaryBy(map.all$denovo71345 ~ Case, data=map.all, FUN=c(mean, sd))
