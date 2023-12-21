# Set working directory and load necessary packages
library(DESeq2)
library(dplyr)
library(stringr)

############################################# RUNNING DESEQ2 #############################################
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
colData <- read.table(file="./data/DEG_S1/expression-design.tsv", header = T)
colData <- arrange(colData, desc(Genotype))
data <- read.table("./data/DEG_S1/S1_star_2pass.tsv", header = T, row.names = 1)
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
write.csv(L12, "data/DEG_S1/SWB_aligned_differential_expression_LMC_between_fieldsites.csv")

## Filter differentially expressed genes between S1_coastal and S1_inland
S1IS1C <- results(dds, list(c("Condition_Inland_vs_Coastal","GenotypeS1.ConditionInland")))
S1 <- as.data.frame(S1IS1C)
S11<-S1[(S1$log2FoldChange>1|S1$log2FoldChange< -1)& !is.na(S1$log2FoldChange),]
S12 <- S11[(S11$padj<0.05)&!is.na(S11$padj),]
S12$gene_model <- row.names(S12)
write.csv(S12, "data/DEG_S1/SWB_aligned_differential_expression_SWB_between_fieldsites.csv")

# ## What genes are different between field sites samples AND between genotypes (i.e. are the genes that are different between sites the same across genotypes)?
field_siteres <- results(dds, name="GenotypeS1.ConditionInland")
field_site <- as.data.frame(field_siteres)
field_site1<-field_site[(field_site$log2FoldChange>1|field_site$log2FoldChange< -1)& !is.na(field_site$log2FoldChange),]
field_site2 <- field_site1[(field_site1$padj<0.05)&!is.na(field_site1$padj),]
field_site2$gene_model <- row.names(field_site2)
write.csv(field_site2, "data/DEG_S1/SWB_aligned_differential_expression_genotype_by_fieldsites.csv")

## What genes are different between genotypes?
genores <- results(dds, name="Genotype_S1_vs_L1")
geno <- as.data.frame(genores)
geno1<-geno[(geno$log2FoldChange>1|geno$log2FoldChange< -1)& !is.na(geno$log2FoldChange),]
geno2 <- geno1[(geno1$padj<0.05)&!is.na(geno1$padj),]
geno2$gene_model <- row.names(geno2)
write.csv(geno2, "data/DEG_S1/SWB_aligned_differential_expression_between_genotypes.csv")

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

S1_outliers_int <- subset(int_genes,  int_genes$gene_model %in% S12$gene_model)

S1_outliers_int2 <- S1_outliers_int %>% 
  inner_join(S12, by ="gene_model", multiple = "all") # Adding S1 outlier data to dataset of syntenic only genes
colnames(S1_outliers_int2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

### Checking for DE in all outlier genes ###

# DE between genotypes and between field sites
S1_outliers <- read.csv("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/S1_aligned/genes/S1_all_data_unique_collapsed.csv")
#L1_outliers <-read.csv("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned/genes/L1_all_data_unique_collapsed.csv")

DE_outliers_S1 <- subset(S1_outliers,   S1_outliers$gene_model %in% field_site2$gene_model)

DE_outliers_S12 <- DE_outliers_S1 %>% 
  inner_join(field_site2, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes
colnames(DE_outliers_S12) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )

# Shared between all outliers and DE in specific genotype

S1_outliers_all <- subset(S1_outliers,  S1_outliers$gene_model %in% S12$gene_model)

S1_outliers_all2 <- S1_outliers_all %>% 
  inner_join(S12, by ="gene_model", multiple = "all") # Adding S1 outlier data to dataset of syntenic only genes
colnames(S1_outliers_all2) <- c("gene_model", "chromosome", "location", "function", "reason", "high_pi", "low_pi", "gstat", "fst", "baseMean",                 
                                "log2FoldChange","lfcSE","stat" ,"pvalue" , "padj" )


# Making plots

gene <- plotCounts(dds,"MgS1_08g09760",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgS1_08g09760")

gene <- plotCounts(dds,"MgS1_05g10260",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgS1_05g10260")

gene <- plotCounts(dds,"MgS1_14g08070",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgS1_14g08070")

gene <- plotCounts(dds,"MgS1_14g08050",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgS1_14g08050")

gene <- plotCounts(dds,"MgS1_08g08930",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgS1_08g08930")


#Myb 2
gene <- plotCounts(dds,"MgS1_08g04810",intgroup=c("Genotype","Condition"),returnData = TRUE)
gene$genotype_condition <- paste(colData$Genotype,colData$Condition,sep="_")
ggplot(gene,aes(x=genotype_condition,y=count,color=genotype_condition)) + geom_boxplot(fill=NA) + geom_jitter() +  theme_bw() + ggtitle("MgS1_08g4810")


#### Karyoploter ####
library(karyoploteR)
library(pasilla)
library(apeglm)
library(ballgown)
library(GenomicRanges)
library(rtracklayer)

# Reading in GFF as a Grange object
S1_grange <- gffReadGR("./data/S1-v1.2_gene.gff")
names(S1_grange) <- S1_grange$ID
S1_grange_geno <- S1_grange
S1_grange_fieldgeno <- S1_grange

#Not filterd
mcols(S1_grange) <- res[names(S1_grange), c("log2FoldChange", "stat", "pvalue", "padj")]


# Differences in S1 vs S1
#genores
mcols(S1_grange_geno) <- genores[names(S1_grange_geno), c("log2FoldChange", "stat", "pvalue", "padj")]
S1_grange_geno

# Differences between genotypes and field sites
#field_siteres
mcols(S1_grange_fieldgeno) <- field_siteres[names(S1_grange_fieldgeno), c("log2FoldChange", "stat", "pvalue", "padj")]

#res_geno <- results(dds, alpha=0.05,name="Genotype_S1_vs_S1")
#res_lfc <- lfcShrink(dds, coef=2, res=res, type="apeglm")

mcols(S1_grange) <- res[names(S1_grange), c("log2FoldChange", "stat", "pvalue", "padj")]

# Making Kryotype plot
S1 <- read.table("./data/S1_chromosome.txt")
colnames(S1) <- c("chr", "start", "stop")
S1_cyto <- read.table("./data/S1-v1.2_gene.gff")
# S1_cytobands <- cbind(S1_cyto$V1, S1_cyto$V4, S1_cyto$V5, S1_cyto$V9)
# S1_cytobands[,4] <- gsub("\\.[1]*$", "", S1_cytobands[,4]) #Removing the .1
# S1_cytobands[,4] <- gsub("\\;.*", "", S1_cytobands[,4]) #Revoving the semi colon and should be everything after but thats not working...
# S1_cytobands[,4] <- gsub("Name=.*", "", S1_cytobands[,4]) #Removing Name= and everything after
# S1_cytobands[,4] <- gsub("ID=", "", S1_cytobands[,4]) #Removing "ID="
# colnames(S1_cytobands) <- c("chr", "start", "end", "name")
# S1_cytobands <- as.data.frame(S1_cytobands)
# write.csv(S1_cytobands, "S1_cytobands.csv")
S1_cytobands <- read.csv("./S1_cytobands.csv")
S1_cytobands <- as.data.frame(S1_cytobands[,2:6])

########################################## ALL #######################################################
#Adding the 10 most differentially expressed genes
ordered <- S1_grange[order(S1_grange$padj, na.last = TRUE),]
kp <- plotKaryotype(genome=S1,cex = .5, cytobands = S1_cytobands)
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal", cex = .5)

#Creating plot of all sig p values
filtered.S1_grange <- S1_grange[!is.na(S1_grange$padj)]
log.pval <- -log10(filtered.S1_grange$padj)
mcols(filtered.S1_grange)$log.pval <- log.pval
filtered.S1_grange

sign.genes <- filtered.S1_grange[filtered.S1_grange$padj < 0.05,]
kp <- plotKaryotype(genome=S1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, cex = .5)

# Range of y values
range(sign.genes$log.pval)

# Tame the floating points we will only plot those with padj less than 0.05
sign.genes <- filtered.S1_grange[filtered.S1_grange$padj < 0.05,]
kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))

# Plotting with log2foldchange because it has a smaller range and is easier to view.
fc.ymax <- ceiling(max((range(sign.genes$log2FoldChange)))) #Adjusting ymin and y max since distribution is above 0
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome=S1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)

# Adding y axis labl
kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Using the size of the point to represent significance.
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome=S1, cex =0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Adding names of top significant genes
top.genes <- ordered[1:20]

kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex =0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", cex= 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)

#Adjusting starting point size
points.top <- 0.8
kp <- plotKaryotype(genome=S1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

#Putting markers exactly where they are located
kp <- plotKaryotype(genome=S1, cex= 0.5)
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

kp <- plotKaryotype(genome=S1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

# Adding regions to the plot
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 0.5, cytobands = S1_cytobands)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(850429, 1032334, 13650670,5329939), end=c(7604769, 1246126, 17847181,7791197))
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))

#Adding gene desnity plots For chromosome 8
png("./figures/S1_Chromosome8_DESeqPlot_all.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr8")
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
kp <- kpPlotDensity(kp, data=S1_grange, window.size = 10e4, data.panel = 2)
dev.off()
# #Altering plotting parameters
# pp <- getDefaultPlotParams(plot.type = 2)
# pp$data2height <- 75
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome="S1", plot.type=2, plot.params = pp, cex = .5)
# kpAddMainTitle(kp, main = "S1 aligned reciprocal transplant experiment", cex = 1)
# #Data panel 1
# kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
# gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
# kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col="#777777")
# kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, label.dist = 0.008, label.color="#444444", line.color = "#777777")
# 
# #Data panel 2
# kp <- kpPlotDensity(kp, data=S1_grange, window.size = 10e4, data.panel = 2)
# 
#Adding gene desnity plots For chromosome 5
png("./figures/S1_Chromosome5_DESeqPlot_all.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr5")
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
kp <- kpPlotDensity(kp, data=S1_grange, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/S1_Chromosome14_DESeqPlot_all.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr14")
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
kp <- kpPlotDensity(kp, data=S1_grange, window.size = 10e4, data.panel = 2)
dev.off()

########################################## Genotype differences #######################################################
#Adding the 10 most differentially expressed genes
ordered <- S1_grange_geno[order(S1_grange_geno$padj, na.last = TRUE),]
kp <- plotKaryotype(genome=S1,cex = .5, cytobands = S1_cytobands)
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal", cex = .5)

#Creating plot of all sig p values
filtered.S1_grange_geno <- S1_grange_geno[!is.na(S1_grange_geno$padj)]
log.pval <- -log10(filtered.S1_grange_geno$padj)
mcols(filtered.S1_grange_geno)$log.pval <- log.pval
filtered.S1_grange_geno

sign.genes <- filtered.S1_grange_geno[filtered.S1_grange_geno$padj < 0.05,]
kp <- plotKaryotype(genome=S1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, cex = .5)

# Range of y values
range(sign.genes$log.pval)

# Tame the floating points we will only plot those with padj less than 0.05
sign.genes <- filtered.S1_grange_geno[filtered.S1_grange_geno$padj < 0.05,]
kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))

# Plotting with log2foldchange because it has a smaller range and is easier to view.
fc.ymax <- ceiling(max((range(sign.genes$log2FoldChange)))) #Adjusting ymin and y max since distribution is above 0
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome=S1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)

# Adding y axis labl
kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Using the size of the point to represent significance.
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome=S1, cex =0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Adding names of top significant genes
top.genes <- ordered[1:20]

kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex =0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", cex= 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)

#Adjusting starting point size
points.top <- 0.8
kp <- plotKaryotype(genome=S1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

#Putting markers exactly where they are located
kp <- plotKaryotype(genome=S1, cex= 0.5)
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

kp <- plotKaryotype(genome=S1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

# Adding regions to the plot
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 0.5, cytobands = S1_cytobands)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Adding gene desnity plots For chromosome 8
png("./figures/S1_Chromosome8_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr8")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.9, r1=1)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_geno, window.size = 10e4, data.panel = 2)
dev.off()
# #Altering plotting parameters
# pp <- getDefaultPlotParams(plot.type = 2)
# pp$data2height <- 75
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome="S1", plot.type=2, plot.params = pp, cex = .5)
# kpAddMainTitle(kp, main = "S1 aligned reciprocal transplant experiment", cex = 1)
# #Data panel 1
# kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
# gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
# kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col="#777777")
# kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, label.dist = 0.008, label.color="#444444", line.color = "#777777")
# 
# #Data panel 2
# kp <- kpPlotDensity(kp, data=S1_grange_geno, window.size = 10e4, data.panel = 2)
# 
#Adding gene desnity plots For chromosome 5
png("./figures/S1_Chromosome5_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr5")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.8, r1=0.9)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_geno, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/S1_Chromosome14_DESeqPlot_geno.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr14")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.8, r1=0.9)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_geno, window.size = 10e4, data.panel = 2)
dev.off()

png("./figures/S1_ChromosomeALL_DESeqPlot_geno.png", width = 2000, height = 1500)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1.2, cytobands = S1_cytobands)
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1.2)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1.2)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()
################################### Geno x Field site differences ######################################

#Adding the 10 most differentially expressed genes
ordered <- S1_grange_fieldgeno[order(S1_grange_fieldgeno$padj, na.last = TRUE),]
kp <- plotKaryotype(genome=S1,cex = .5, cytobands = S1_cytobands)
kp <- kpPlotMarkers(kp, ordered[1:10], labels = names(ordered[1:10]), text.orientation = "horizontal", cex = .5)

#Creating plot of all sig p values
filtered.S1_grange_fieldgeno <- S1_grange_fieldgeno[!is.na(S1_grange_fieldgeno$padj)]
log.pval <- -log10(filtered.S1_grange_fieldgeno$padj)
mcols(filtered.S1_grange_fieldgeno)$log.pval <- log.pval
filtered.S1_grange_fieldgeno

sign.genes <- filtered.S1_grange_fieldgeno[filtered.S1_grange_fieldgeno$padj < 0.05,]
kp <- plotKaryotype(genome=S1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, cex = .5)

# Range of y values
range(sign.genes$log.pval)

# Tame the floating points we will only plot those with padj less than 0.05
sign.genes <- filtered.S1_grange_fieldgeno[filtered.S1_grange_fieldgeno$padj < 0.05,]
kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log.pval, ymax=max(sign.genes$log.pval))

# Plotting with log2foldchange because it has a smaller range and is easier to view.
fc.ymax <- ceiling(max((range(sign.genes$log2FoldChange)))) #Adjusting ymin and y max since distribution is above 0
fc.ymin <- -fc.ymax

kp <- plotKaryotype(genome=S1, cex = .5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)

# Adding y axis labl
kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Using the size of the point to represent significance.
cex.val <- sqrt(sign.genes$log.pval)/3
kp <- plotKaryotype(genome=S1, cex =0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex = 0.5)

#Adding names of top significant genes
top.genes <- ordered[1:20]

kp <- plotKaryotype(genome=S1, cex = 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, cex =0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", cex= 0.5)
kpAddLabels(kp, labels = "log2 Fold Change", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, cex= 0.5)

#Adjusting starting point size
points.top <- 0.8
kp <- plotKaryotype(genome=S1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

#Putting markers exactly where they are located
kp <- plotKaryotype(genome=S1, cex= 0.5)
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

kp <- plotKaryotype(genome=S1, cex= 0.5)
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 0.5)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "horizontal", r0=points.top, cex= 0.5)

# Adding regions to the plot
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 0.5, cytobands = S1_cytobands)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Adding gene desnity plots For chromosome 8
png("./figures/S1_Chromosome8_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr8")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.9, r1=1)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()
# #Altering plotting parameters
# pp <- getDefaultPlotParams(plot.type = 2)
# pp$data2height <- 75
# pp$ideogramheight <- 10
# kp <- plotKaryotype(genome="S1", plot.type=2, plot.params = pp, cex = .5)
# kpAddMainTitle(kp, main = "S1 aligned reciprocal transplant experiment", cex = 1)
# #Data panel 1
# kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
# gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
# kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col="#777777")
# kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, label.dist = 0.008, label.color="#444444", line.color = "#777777")
# 
# #Data panel 2
# kp <- kpPlotDensity(kp, data=S1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
# 
#Adding gene desnity plots For chromosome 5
png("./figures/S1_Chromosome5_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr5")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.8, r1=0.9)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome 14
png("./figures/S1_Chromosome14_DESeqPlot_genobyfield.png", width = 1500, height = 700)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1, cytobands = S1_cytobands, chromosomes = "chr14")
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.8, r1=0.9)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()

#Adding gene desnity plots For chromosome all
png("./figures/S1_ChromosomeALL_DESeqPlot_genobyfield.png", width = 2000, height = 1500)
kp <- plotKaryotype(genome=S1, plot.type=2, cex= 1.2, cytobands = S1_cytobands)
#Data panel 1
kpPoints(kp, data=sign.genes, y=sign.genes$log2FoldChange, cex=cex.val, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, col=sign.col)
kpAxis(kp, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1.2)
kpAddLabels(kp, labels = "log2 FC", srt=90, pos=1, label.margin = 0.04, ymax=fc.ymax, ymin=fc.ymin, r1=points.top, cex= 1.2)
gene.mean <- start(top.genes) + (end(top.genes) - start(top.genes))/2
kpSegments(kp, chr=as.character(seqnames(top.genes)), x0=gene.mean, x1=gene.mean, y0=top.genes$log2FoldChange, y1=fc.ymax, ymax=fc.ymax, ymin=fc.ymin, r1=points.top)
kpPlotMarkers(kp, top.genes, labels = names(top.genes), text.orientation = "vertical", r0=points.top, cex= 1)
df <- data.frame(chr=c("chr8", "chr8", "chr5","chr14"), start=c(858736, 5998986, 14125309,5892556), end=c(6465310, 6260288, 18136144,8677481))
kpPlotRegions(kp, data=df, col=c("#1B9E77","#D95F02","#7570B3","#E7298A"), border=c("#1B9E77","#D95F02","#7570B3","#E7298A"), r0=0.7, r1=0.85)

#Data panel 2
kp <- kpPlotDensity(kp, data=S1_grange_fieldgeno, window.size = 10e4, data.panel = 2)
dev.off()