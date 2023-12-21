# Set working directory and load necessary packages
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
############################################# RUNNING DESEQ2 #############################################
setwd("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/")
##Preparing files for DESeq 
# ff <- list.files( path = "./data/DEG_L1/", pattern = "*L1_star_2pass", full.names = TRUE )
# counts.files <- lapply( ff, read.table, skip = 4 )
# counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) )
# ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
# ff <- gsub( "[.]/counts/", "", ff )
# colnames(counts) <- ff
# row.names(counts) <- counts.files[[1]]$V1
# colnames(counts) <- str_replace_all(colnames(counts),"./all-counts-V2/","")
# colnames(counts) <- str_replace_all(colnames(counts),"ReadsPerGene.out.tab","")
colData <- read.table(file="./data/DEG_L1/expression-design.tsv", header = T)
colData <- arrange(colData, desc(Genotype))
data <- read.table("./data/DEG_L1/L1_star_2pass.tsv", header = T, row.names = 1)
#counts <- counts[,c(7,8,9,10,11,12,1,2,3,4,5,6)]

##Making DESeq Data Set
dds = DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~Genotype + Condition + Genotype:Condition)

#Make sure that dds is properly formatted
as.data.frame(colData(dds))

##Run DESeq
dds <- DESeq(dds)

res <- results(dds)

############################### EXTRACTING AND FILTERING SPECIFIC SAMPLES ################################
## Double check result names
resultsNames(dds)

## Filter differentially expressed genes between L1_coastal and L1_inland
L1IL1C <- results(dds, contrast=c("Condition", "Coastal", "Inland"))
L1 <- as.data.frame(L1IL1C)
L11<-L1[(L1$log2FoldChange>1|L1$log2FoldChange< -1)& !is.na(L1$log2FoldChange),]
L12 <- L11[(L11$padj<0.05)&!is.na(L11$padj),]
L12$gene_model <- row.names(L12)
write.csv(L12, "data/DEG_L1/LMC_aligned_differential_expression_LMC_between_fieldsites.csv")

## Filter differentially expressed genes between S1_coastal and S1_inland
S1IS1C <- results(dds, list(c("Condition_Inland_vs_Coastal","GenotypeS1.ConditionInland")))
S1 <- as.data.frame(S1IS1C)
S11<-S1[(S1$log2FoldChange>1|S1$log2FoldChange< -1)& !is.na(S1$log2FoldChange),]
S12 <- S11[(S11$padj<0.05)&!is.na(S11$padj),]
S12$gene_model <- row.names(S12)
write.csv(S12, "data/DEG_L1/LMC_aligned_differential_expression_SWB_between_fieldsites.csv")

# ## What genes are different between field sites samples AND between genotypes (i.e. are the genes that are different between sites the same across genotypes)?
field_siteres <- results(dds, name="GenotypeS1.ConditionInland")
field_site <- as.data.frame(field_siteres)
field_site1<-field_site[(field_site$log2FoldChange>1|field_site$log2FoldChange< -1)& !is.na(field_site$log2FoldChange),]
field_site2 <- field_site1[(field_site1$padj<0.05)&!is.na(field_site1$padj),]
field_site2$gene_model <- row.names(field_site2)
write.csv(field_site2, "data/DEG_L1/LMC_aligned_differential_expression_genotype_by_fieldsites.csv")

## What genes are different between genotypes?
genores <- results(dds, name="Genotype_S1_vs_L1")
geno <- as.data.frame(genores)
geno1<-geno[(geno$log2FoldChange>1|geno$log2FoldChange< -1)& !is.na(geno$log2FoldChange),]
geno2 <- geno1[(geno1$padj<0.05)&!is.na(geno1$padj),]
geno2$gene_model <- row.names(geno2)
write.csv(geno2, "data/DEG_L1/LMC_aligned_differential_expression_between_genotypes.csv")

## Pulling out intersting genes

# Genes we found interesting because they are ...
  # 1. orthologs with genes found in GOuld et al 2017
  # 2. near breakpoint regions
  # 3. intersting outliers

int_genes <- read.csv("./data/Intersting_genes_for_expression.csv", header = T)

#### Shared between our interesting outliers and the DEGs between field sites between genotypes ####
DE_outliers_int <- subset(int_genes,   int_genes$gene_model %in% field_site2$gene_model)

DE_outliers_int2 <- DE_outliers_int %>% 
  inner_join(field_site2, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(DE_outliers_int2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

## Shared between our interesting outliers and DE in specific genotype

L1_outliers_int <- subset(int_genes,  int_genes$gene_model %in% L12$gene_model)

L1_outliers_int2 <- L1_outliers_int %>% 
  inner_join(L12, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(L1_outliers_int2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

### Checking for DE in all outlier genes ###

# DE between genotypes and between field sites
L1_outliers <- read.csv("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned/genes/data_files/L1_collapsed_alldata.csv")
colnames(L1_outliers) <- c("X","chr","location","gene_model","Function","start", "stop" ,"gstat","fst" , "low_pi" ,"high_pi" )
#L1_outliers <-read.csv("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned/genes/L1_all_data_unique_collapsed.csv")

DE_outliers_L1 <- subset(L1_outliers,   L1_outliers$gene_model %in% field_site2$gene_model)

DE_outliers_L12 <- DE_outliers_L1 %>% 
  inner_join(field_site2, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(DE_outliers_L12) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                               "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

# Shared between all outliers and DE in specific genotype

L1_outliers_all <- subset(L1_outliers,  L1_outliers$gene_model %in% L12$gene_model)

L1_outliers_all2 <- L1_outliers_all %>% 
  inner_join(L12, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(L1_outliers_all2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

# Making plots

gene <- plotCounts(dds,"MgL1_08g16120",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g16120")

gene <- plotCounts(dds,"MgL1_14g10910",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_14g10910")

gene <- plotCounts(dds,"MgL1_05g11230",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_05g11230")

gene <- plotCounts(dds,"MgL1_05g15660",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_05g15660")

gene <- plotCounts(dds,"MgL1_08g06810",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g06810")

gene <- plotCounts(dds,"MgL1_08g23640",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g23640")

gene <- plotCounts(dds,"MgL1_08g02520",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g02520")


# Myb anthocyanin

# amMyb2
# This one has been duplicated.
# Myb1, Myb2, and 3 were tandem duplicates (happened in mg specifically) but evolved seperately. (L1 and S1 both have)
# Common ancester did not have them.
# Expressed higher in L1 than S1 in all locations
# Based on Billies its only happening in the leaf.

gene <- plotCounts(dds,"MgL1_08g07000",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgL1_08g07000")

#### Karyoploter ####
library(karyoploteR)
library(pasilla)
library(apeglm)
library(ballgown)
library(GenomicRanges)
library(rtracklayer)

# Reading in GFF as a Grange object
L1_grange <- gffReadGR("./data/L1-v1.2_gene.gff")
names(L1_grange) <- L1_grange$ID
L1_grange_geno <- L1_grange
L1_grange_fieldgeno <- L1_grange

#Not filterd
mcols(L1_grange) <- res[names(L1_grange), c("log2FoldChange", "stat", "pvalue", "padj")]


# Differences in L1 vs S1
#genores
mcols(L1_grange_geno) <- genores[names(L1_grange_geno), c("log2FoldChange", "stat", "pvalue", "padj")]
L1_grange_geno

# Differences between genotypes and field sites
#field_siteres
mcols(L1_grange_fieldgeno) <- field_siteres[names(L1_grange_fieldgeno), c("log2FoldChange", "stat", "pvalue", "padj")]

#res_geno <- results(dds, alpha=0.05,name="Genotype_S1_vs_L1")
#res_lfc <- lfcShrink(dds, coef=2, res=res, type="apeglm")

mcols(L1_grange) <- res[names(L1_grange), c("log2FoldChange", "stat", "pvalue", "padj")]

# Making Kryotype plot
L1 <- read.table("./data/L1_chromosome.txt")
colnames(L1) <- c("chr", "start", "stop")
L1_cyto <- read.table("./data/L1-v1.2_gene.gff")
# L1_cytobands <- cbind(L1_cyto$V1, L1_cyto$V4, L1_cyto$V5, L1_cyto$V9)
# L1_cytobands[,4] <- gsub("\\.[1]*$", "", L1_cytobands[,4]) #Removing the .1
# L1_cytobands[,4] <- gsub("\\;.*", "", L1_cytobands[,4]) #Revoving the semi colon and should be everything after but thats not working...
# L1_cytobands[,4] <- gsub("Name=.*", "", L1_cytobands[,4]) #Removing Name= and everything after
# L1_cytobands[,4] <- gsub("ID=", "", L1_cytobands[,4]) #Removing "ID="
# colnames(L1_cytobands) <- c("chr", "start", "end", "name")
# L1_cytobands <- as.data.frame(L1_cytobands)
# write.csv(L1_cytobands, "L1_cytobands.csv")
L1_cytobands <- read.csv("./L1_cytobands.csv")
L1_cytobands <- as.data.frame(L1_cytobands[,2:6])

########################################## ALL #######################################################
#Adding the 10 most differentially expressed genes
ordered <- L1_grange[order(L1_grange$padj, na.last = TRUE),]
kp <- plotKaryotype(genome=L1,cex = .5, cytobands = L1_cytobands)
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal", cex = .5)

#Creating plot of all sig p values
filtered.L1_grange <- L1_grange[!is.na(L1_grange$padj)]
log.pval <- -log10(filtered.L1_grange$padj)
mcols(filtered.L1_grange)$log.pval <- log.pval
filtered.L1_grange

sign.genes <- filtered.L1_grange[filtered.L1_grange$padj < 0.05,]
kp <- plotKaryotype(genome=L1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, cex = .5)

# Range of y values
range(sign.genes$log.pval)

# Tame the floating points we will only plot those with padj less than 0.05
sign.genes <- filtered.L1_grange[filtered.L1_grange$padj < 0.05,]
kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))

# Plotting with log2foldchange because it has a smaller range and is easier to view.
fc.ymax <- ceiling(max((range(sign.genes$log2FoldChange)))) #Adjusting ymin and y max since distribution is above 0
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome=L1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)

# Adding y axis labl
kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Using the size of the point to represent significance.
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome=L1, cex =0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Adding names of top significant genes
top.genes <- ordered[1:20]

kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex =0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", cex= 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)

#Adjusting starting point size
points.top <- 0.8
kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

#Putting markers exactly where they are located
kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)


#Colog to differentiate between over and underexpressed genes
col.over <- "#FFBD07AA"
col.under <- "#00A6EDAA"
sign.col <- rep(col.over, length(sign.genes))
sign.col[sign.genes$log2FoldChange<0] <- col.under

kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

# Adding regions to the plot
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 0.5, cytobands = L1_cytobands)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.3, r1=0.55)

#Adding gene desnity plots For chromosome 8
png("./figures/L1_Chromosome8_DESeqPlot_all.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr8")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange, window.size = 10e4, data.panel = 2)
dev.off()
# #Altering plotting parameters
# pp <- getDefaultPlotParams(plot.type = 2)
# pp$data2height <- 75
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome="L1", plot.type=2, plot.params = pp, cex = .5)
# kpAddMainTitle(kp, main = "L1 aligned reciprocal transplant experiment", cex = 1)
# #Data panel 1
# kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
# gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
# kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col="#777777")
# kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, label.dist = 0.008, label.color="#444444", line.color = "#777777")
# 
# #Data panel 2
# kp <- kpPlotDensity(kp, data=L1_grange, window.size = 10e4, data.panel = 2)
# 
#Adding gene desnity plots For chromosome 5
png("./figures/L1_Chromosome5_DESeqPlot_all.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr5")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/L1_Chromosome14_DESeqPlot_all.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr14")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange, window.size = 10e4, data.panel = 2)
dev.off()

########################################## Genotype differences #######################################################
#Adding the 10 most differentially expressed genes
ordered <- L1_grange_geno[order(L1_grange_geno$padj, na.last = TRUE),]
kp <- plotKaryotype(genome=L1,cex = .5, cytobands = L1_cytobands)
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal", cex = .5)

#Creating plot of all sig p values
filtered.L1_grange_geno <- L1_grange_geno[!is.na(L1_grange_geno$padj)]
log.pval <- -log10(filtered.L1_grange_geno$padj)
mcols(filtered.L1_grange_geno)$log.pval <- log.pval
filtered.L1_grange_geno

sign.genes <- filtered.L1_grange_geno[filtered.L1_grange_geno$padj < 0.05,]
kp <- plotKaryotype(genome=L1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, cex = .5)

# Range of y values
range(sign.genes$log.pval)

# Tame the floating points we will only plot those with padj less than 0.05
sign.genes <- filtered.L1_grange_geno[filtered.L1_grange_geno$padj < 0.05,]
kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))

# Plotting with log2foldchange because it has a smaller range and is easier to view.
fc.ymax <- ceiling(max((range(sign.genes$log2FoldChange)))) #Adjusting ymin and y max since distribution is above 0
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome=L1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)

# Adding y axis labl
kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Using the size of the point to represent significance.
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome=L1, cex =0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Adding names of top significant genes
top.genes <- ordered[1:20]

kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex =0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", cex= 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)

#Adjusting starting point size
points.top <- 0.8
kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

#Putting markers exactly where they are located
kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)


#Colog to differentiate between over and underexpressed genes
col.over <- "#FFBD07AA"
col.under <- "#00A6EDAA"
sign.col <- rep(col.over, length(sign.genes))
sign.col[sign.genes$log2FoldChange<0] <- col.under

kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

# Adding regions to the plot
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 0.5, cytobands = L1_cytobands)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.3, r1=0.55)

#Adding gene desnity plots For chromosome 8
png("./figures/L1_Chromosome8_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr8")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_geno, window.size = 10e4, data.panel = 2)
dev.off()
# #Altering plotting parameters
# pp <- getDefaultPlotParams(plot.type = 2)
# pp$data2height <- 75
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome="L1", plot.type=2, plot.params = pp, cex = .5)
# kpAddMainTitle(kp, main = "L1 aligned reciprocal transplant experiment", cex = 1)
# #Data panel 1
# kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
# gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
# kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col="#777777")
# kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, label.dist = 0.008, label.color="#444444", line.color = "#777777")
# 
# #Data panel 2
# kp <- kpPlotDensity(kp, data=L1_grange_geno, window.size = 10e4, data.panel = 2)
# 
#Adding gene desnity plots For chromosome 5
png("./figures/L1_Chromosome5_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr5")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_geno, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/L1_Chromosome14_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr14")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_geno, window.size = 10e4, data.panel = 2)
dev.off()

png("./figures/L1_ChromosomeALL_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands)
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()
################################### Geno x Field site differences ######################################

#Adding the 10 most differentially expressed genes
ordered <- L1_grange_fieldgeno[order(L1_grange_fieldgeno$padj, na.last = TRUE),]
kp <- plotKaryotype(genome=L1,cex = .5, cytobands = L1_cytobands)
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal", cex = .5)

#Creating plot of all sig p values
filtered.L1_grange_fieldgeno <- L1_grange_fieldgeno[!is.na(L1_grange_fieldgeno$padj)]
log.pval <- -log10(filtered.L1_grange_fieldgeno$padj)
mcols(filtered.L1_grange_fieldgeno)$log.pval <- log.pval
filtered.L1_grange_fieldgeno

sign.genes <- filtered.L1_grange_fieldgeno[filtered.L1_grange_fieldgeno$padj < 0.05,]
kp <- plotKaryotype(genome=L1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, cex = .5)

# Range of y values
range(sign.genes$log.pval)

# Tame the floating points we will only plot those with padj less than 0.05
sign.genes <- filtered.L1_grange_fieldgeno[filtered.L1_grange_fieldgeno$padj < 0.05,]
kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))

# Plotting with log2foldchange because it has a smaller range and is easier to view.
fc.ymax <- ceiling(max((range(sign.genes$log2FoldChange)))) #Adjusting ymin and y max since distribution is above 0
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome=L1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)

# Adding y axis labl
kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Using the size of the point to represent significance.
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome=L1, cex =0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Adding names of top significant genes
top.genes <- ordered[1:20]

kp <- plotKaryotype(genome=L1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex =0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", cex= 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)

#Adjusting starting point size
points.top <- .5
kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

#Putting markers exactly where they are located
kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)


#Colog to differentiate between over and underexpressed genes
col.over <- "#FFBD07AA"
col.under <- "#00A6EDAA"
sign.col <- rep(col.over, length(sign.genes))
sign.col[sign.genes$log2FoldChange<0] <- col.under

kp <- plotKaryotype(genome=L1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

# Adding regions to the plot
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 0.5, cytobands = L1_cytobands)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.3, r1=0.55)

#Adding gene desnity plots For chromosome 8
png("./figures/L1_Chromosome8_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr8")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()
# #Altering plotting parameters
# pp <- getDefaultPlotParams(plot.type = 2)
# pp$data2height <- 75
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome="L1", plot.type=2, plot.params = pp, cex = .5)
# kpAddMainTitle(kp, main = "L1 aligned reciprocal transplant experiment", cex = 1)
# #Data panel 1
# kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
# gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
# kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col="#777777")
# kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, label.dist = 0.008, label.color="#444444", line.color = "#777777")
# 
# #Data panel 2
# kp <- kpPlotDensity(kp, data=L1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
# 
#Adding gene desnity plots For chromosome 5
png("./figures/L1_Chromosome5_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr5")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/L1_Chromosome14_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands, chromosomes = "chr14")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/L1_ChromosomeALL_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=L1, plot.type=2, cex= 1, cytobands = L1_cytobands)
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=L1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()




