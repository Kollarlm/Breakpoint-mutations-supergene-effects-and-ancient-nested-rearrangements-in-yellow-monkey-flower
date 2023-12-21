## This script takes the output of annovar -> high impact filtering script
## and will find the number of genes within inversions that have snps that could affect expression differences.
## This script can but is not meant to calculate number of high impact snps within an inversion. It only takes
## unique genes so if multiple snps are within a gene or type of snp it will get collapsed. To calculate snps per 
## inversion you would have to remove the "unique" function.
## Author: Leslie M. Kollar
## Date: May 2, 2023

# packages
library(tidyr)
library(dplyr)

#Setwd

path <- "/Users/lesliekollar/Desktop/mimulus_genome_R/mimulus_genome/Mimulus_guttatus_nanopore_genome_v1.2"
setwd(path)

# Reading in genes within the inversion
chr_5_genes1 <- read.table("L1_genes_chr5_inversion.gff")
chr_5_genes2 <- read.table("S1_genes_chr5_inversion.gff")
chr_5_genes <- rbind(chr_5_genes1, chr_5_genes2)
chr_5_genes$V9 <- gsub(";Name=MgL1_05g[^\\d]+", "", chr_5_genes$V9)
chr_5_genes$V9 <- gsub(";Name=MgS1_05g[^\\d]+", "", chr_5_genes$V9)
colnames(chr_5_genes)[colnames(chr_5_genes)== "V9"] <- 'id'
chr_5_genes$id <- gsub("ID=", "", chr_5_genes$id)

chr_8_genes1 <- read.table("L1_genes_chr8_inversion.gff")
chr_8_genes2 <- read.table("S1_genes_chr8_inversion.gff")
chr_8_genes <- rbind(chr_8_genes1, chr_8_genes2)
chr_8_genes$V9 <- gsub(";Name=MgL1_08g[^\\d]+", "", chr_8_genes$V9)
chr_8_genes$V9 <- gsub(";Name=MgS1_08g[^\\d]+", "", chr_8_genes$V9)
colnames(chr_8_genes)[colnames(chr_8_genes)== "V9"] <- 'id'
chr_8_genes$id <- gsub("ID=", "", chr_8_genes$id)

chr_8small_genes1 <- read.table("L1_genes_chr8_inversion_small.gff")
chr_8small_genes2 <- read.table("S1_genes_chr8_inversion_small.gff")
chr_8small_genes <- rbind(chr_8small_genes1, chr_8small_genes2)
chr_8small_genes$V9 <- gsub(";Name=MgL1_08g[^\\d]+", "", chr_8small_genes$V9)
chr_8small_genes$V9 <- gsub(";Name=MgS1_08g[^\\d]+", "", chr_8small_genes$V9)
colnames(chr_8small_genes)[colnames(chr_8small_genes)== "V9"] <- 'id'
chr_8small_genes$id <- gsub("ID=", "", chr_8small_genes$id)

chr_14_genes1 <- read.table("L1_genes_chr14_inversion.gff")
chr_14_genes2 <- read.table("S1_genes_chr14_inversion.gff")
chr_14_genes <- rbind(chr_14_genes1, chr_14_genes2)
chr_14_genes$V9 <- gsub(";Name=MgL1_14g[^\\d]+", "", chr_14_genes$V9)
chr_14_genes$V9 <- gsub(";Name=MgS1_14g[^\\d]+", "", chr_14_genes$V9)
colnames(chr_14_genes)[colnames(chr_14_genes)== "V9"] <- 'id'
chr_14_genes$id <- gsub("ID=", "", chr_14_genes$id)

# Read in high impact snp annotation file

S1_snps <- read.csv("./data/high-impact-S1_nohets.csv", sep = ",", header = T)
colnames(S1_snps) <- c("SNP_type","Gene_impacted","Zygosity")
S1_snps$id <- gsub(":Mg.*", "", S1_snps$Gene_impacted)
S1_snps$id <- gsub("\\(Mg.*", "", S1_snps$id)
S1_snps_unique <- as.data.frame(unique(S1_snps$id))

L1_snps <- read.csv("./data/high-impact-L1_nohets.csv", sep = ",", header = T)
colnames(L1_snps) <- c("SNP_type","Gene_impacted","Zygosity")
L1_snps$id <- gsub(":Mg.*", "", L1_snps$Gene_impacted)
L1_snps$id <- gsub("\\(Mg.*", "", L1_snps$id)
L1_snps_unique <- as.data.frame(unique(L1_snps$id))

# SNP types
unique(S1_snps$SNP_type)
unique(L1_snps$SNP_type)

# "unknown"                 "nonsynonymous SNV"       "nonframeshift insertion" "stopgain"               
# "frameshift deletion"     "nonframeshift deletion"  "frameshift insertion"    "stoploss"               
# "intronic"                "UTR3"                    "UTR5"                    "splicing"  

# Filtering based on type of SNP and pulling out genes with 1 or more snps with in the inversion

##nonsynonymous SNV within genes
L1_genes_w_snps_nonsynSNV <- L1_snps %>%  ## Only snps that fall within known S1 genes (location because we aligned to S1 genome) which should be the case but just in case.
  filter(grepl("MgS1_", Gene_impacted)) %>%
  filter(SNP_type == "nonsynonymous SNV") ## Filter out snps that are a specific type
nrow(L1_genes_w_snps_nonsynSNV) # Number of nonsynonymous snps in total in L1

L1_gene_list_nonsynSNV <- as.data.frame(unique(L1_genes_w_snps_nonsynSNV$id)) # Even with in the type of snp there can be repeats (i.e. multiple snps). Use the unique function to count the number of genes rather than the number of snps
colnames(L1_gene_list_nonsynSNV)[1] <- 'id'
nrow(L1_gene_list_nonsynSNV)

L1_chr8small_inv_nonsynSNV <- L1_gene_list_nonsynSNV %>% # Pull out genes from master list of genes with snps that fall within inversion.
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_nonsynSNV)
L1_chr8_inv_nonsynSNV <- L1_gene_list_nonsynSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_nonsynSNV)
L1_chr5_inv_nonsynSNV <- L1_gene_list_nonsynSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_nonsynSNV)
L1_chr14_inv_nonsynSNV <- L1_gene_list_nonsynSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_nonsynSNV)

L1_chr8small_inv_nonsynSNV <- L1_genes_w_snps_nonsynSNV %>% # Pull out genes from master list of genes with snps that fall within inversion.
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_nonsynSNV)
L1_chr8_inv_nonsynSNV <- L1_genes_w_snps_nonsynSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_nonsynSNV)
L1_chr5_inv_nonsynSNV <- L1_genes_w_snps_nonsynSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_nonsynSNV)
L1_chr14_inv_nonsynSNV <- L1_genes_w_snps_nonsynSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_nonsynSNV)

# Unknown
L1_genes_w_snps_unknown <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "unknown")
L1_gene_list_unknown <- as.data.frame(unique(L1_genes_w_snps_unknown$id))
colnames(L1_gene_list_unknown)[1] <- 'id'
L1_chr8small_inv_unknown <- L1_gene_list_unknown %>% 
  filter(id %in% chr_8small_genes$id)
L1_chr8_inv_unknown <- L1_gene_list_unknown %>% 
  filter(id %in% chr_8_genes$id)
L1_chr5_inv_unknown <- L1_gene_list_unknown %>% 
  filter(id %in% chr_5_genes$id)
L1_chr14_inv_unknown <- L1_gene_list_unknown %>% 
  filter(id %in% chr_14_genes$id)

#Nonframeshift insertion
L1_genes_w_snps_nonframeshift_insertion <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "nonframeshift insertion")
L1_gene_list_nonframeshift_insertion <- as.data.frame(unique(L1_genes_w_snps_nonframeshift_insertion$id)) 
colnames(L1_gene_list_nonframeshift_insertion)[1] <- 'id'
L1_chr8small_inv_nonframeshift_insertion <- L1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_8small_genes$id)
L1_chr8_inv_nonframeshift_insertion <- L1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_8_genes$id)
L1_chr5_inv_nonframeshift_insertion <- L1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_5_genes$id)
L1_chr14_inv_nonframeshift_insertion <- L1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_14_genes$id)

# Stop gain 
L1_genes_w_snps_stopgain <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "stopgain")
nrow(L1_genes_w_snps_stopgain )

L1_gene_list_stopgain <- as.data.frame(unique(L1_genes_w_snps_stopgain$id)) 
colnames(L1_gene_list_stopgain)[1] <- 'id'
nrow(L1_gene_list_stopgain)

L1_chr8small_inv_stopgain <- L1_gene_list_stopgain %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_stopgain)
L1_chr8_inv_stopgain <- L1_gene_list_stopgain %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_stopgain)
L1_chr5_inv_stopgain <- L1_gene_list_stopgain %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_stopgain)
L1_chr14_inv_stopgain <- L1_gene_list_stopgain %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_stopgain)

L1_chr8small_inv_stopgain <- L1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_stopgain)
L1_chr8_inv_stopgain <- L1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_stopgain)
L1_chr5_inv_stopgain <- L1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_stopgain)
L1_chr14_inv_stopgain <- L1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_stopgain)

#Frameshift deletion
L1_genes_w_snps_frameshift_deletion <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "frameshift deletion")
nrow(L1_genes_w_snps_frameshift_deletion)

L1_gene_list_frameshift_deletion <- as.data.frame(unique(L1_genes_w_snps_frameshift_deletion$id)) 
colnames(L1_gene_list_frameshift_deletion)[1] <- 'id'
nrow(L1_gene_list_frameshift_deletion)

L1_chr8small_inv_frameshift_deletion <- L1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_frameshift_deletion)
L1_chr8_inv_frameshift_deletion <- L1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_frameshift_deletion)
L1_chr5_inv_frameshift_deletion <- L1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_frameshift_deletion)
L1_chr14_inv_frameshift_deletion <- L1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_frameshift_deletion)

L1_chr8small_inv_frameshift_deletion <- L1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_frameshift_deletion)
L1_chr8_inv_frameshift_deletion <- L1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_frameshift_deletion)
L1_chr5_inv_frameshift_deletion <- L1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_frameshift_deletion)
L1_chr14_inv_frameshift_deletion <- L1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_frameshift_deletion)

# Non frameshift deletion
L1_genes_w_snps_nonframeshift_deletion <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "nonframeshift deletion")
L1_gene_list_nonframeshift_deletion <- as.data.frame(unique(L1_genes_w_snps_nonframeshift_deletion$id)) 
colnames(L1_gene_list_nonframeshift_deletion)[1] <- 'id'
L1_chr8small_inv_nonframeshift_deletion <- L1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_8small_genes$id)
L1_chr8_inv_nonframeshift_deletion <- L1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_8_genes$id)
L1_chr5_inv_nonframeshift_deletion <- L1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_5_genes$id)
L1_chr14_inv_nonframeshift_deletion <- L1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_14_genes$id)

# Frameshift insertion
L1_genes_w_snps_frameshift_insertion <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "frameshift insertion")
nrow(L1_genes_w_snps_frameshift_insertion)

L1_gene_list_frameshift_insertion <- as.data.frame(unique(L1_genes_w_snps_frameshift_insertion$id)) 
colnames(L1_gene_list_frameshift_insertion)[1] <- 'id'
nrow(L1_gene_list_frameshift_insertion)

L1_chr8small_inv_frameshift_insertion <- L1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_frameshift_insertion)
L1_chr8_inv_frameshift_insertion <- L1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_frameshift_insertion)
L1_chr5_inv_frameshift_insertion <- L1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_frameshift_insertion)
L1_chr14_inv_frameshift_insertion <- L1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_frameshift_insertion)

L1_chr8small_inv_frameshift_insertion <- L1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_frameshift_insertion)
L1_chr8_inv_frameshift_insertion <- L1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_frameshift_insertion)
L1_chr5_inv_frameshift_insertion <- L1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_frameshift_insertion)
L1_chr14_inv_frameshift_insertion <- L1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_frameshift_insertion)

# Stop loss
L1_genes_w_snps_stoploss <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "stoploss")
nrow(L1_genes_w_snps_stoploss)

L1_gene_list_stoploss <- as.data.frame(unique(L1_genes_w_snps_stoploss$id))
colnames(L1_gene_list_stoploss)[1] <- 'id'
nrow(L1_gene_list_stoploss)

L1_chr8small_inv_stoploss <- L1_gene_list_stoploss %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_stoploss)
L1_chr8_inv_stoploss <- L1_gene_list_stoploss %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_stoploss)
L1_chr5_inv_stoploss <- L1_gene_list_stoploss %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_stoploss)
L1_chr14_inv_stoploss <- L1_gene_list_stoploss %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_stoploss)

L1_chr8small_inv_stoploss <- L1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_stoploss)
L1_chr8_inv_stoploss <- L1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_stoploss)
L1_chr5_inv_stoploss <- L1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_stoploss)
L1_chr14_inv_stoploss <- L1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_stoploss)

#Intronic
L1_genes_w_snps_intronic <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "intronic")
L1_gene_list_intronic <- as.data.frame(unique(L1_genes_w_snps_intronic$id))
colnames(L1_gene_list_intronic)[1] <- 'id'
L1_chr8small_inv_intronic <- L1_gene_list_intronic %>% 
  filter(id %in% chr_8small_genes$id)
L1_chr8_inv_intronic <- L1_gene_list_intronic %>% 
  filter(id %in% chr_8_genes$id)
L1_chr5_inv_intronic <- L1_gene_list_intronic %>% 
  filter(id %in% chr_5_genes$id)
L1_chr14_inv_intronic <- L1_gene_list_intronic %>% 
  filter(id %in% chr_14_genes$id)

# UTR3
L1_genes_w_snps_UTR3 <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "UTR3")
L1_gene_list_UTR3 <- as.data.frame(unique(L1_genes_w_snps_UTR3$id))
colnames(L1_gene_list_UTR3)[1] <- 'id'
L1_chr8small_inv_UTR3 <- L1_gene_list_UTR3 %>% 
  filter(id %in% chr_8small_genes$id)
L1_chr8_inv_UTR3 <- L1_gene_list_UTR3 %>% 
  filter(id %in% chr_8_genes$id)
L1_chr5_inv_UTR3 <- L1_gene_list_UTR3 %>% 
  filter(id %in% chr_5_genes$id)
L1_chr14_inv_UTR3 <- L1_gene_list_UTR3 %>% 
  filter(id %in% chr_14_genes$id)

#UTR5
L1_genes_w_snps_UTR5 <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "UTR5")
L1_gene_list_UTR5 <- as.data.frame(unique(L1_genes_w_snps_UTR5$id))
colnames(L1_gene_list_UTR5)[1] <- 'id'
L1_chr8small_inv_UTR5 <- L1_gene_list_UTR5 %>% 
  filter(id %in% chr_8small_genes$id)
L1_chr8_inv_UTR5 <- L1_gene_list_UTR5 %>% 
  filter(id %in% chr_8_genes$id)
L1_chr5_inv_UTR5 <- L1_gene_list_UTR5 %>% 
  filter(id %in% chr_5_genes$id)
L1_chr14_inv_UTR5 <- L1_gene_list_UTR5 %>% 
  filter(id %in% chr_14_genes$id)

#Splicing
L1_genes_w_snps_splicing <- L1_snps %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "splicing")
nrow(L1_genes_w_snps_splicing)

L1_gene_list_splicing <- as.data.frame(unique(L1_genes_w_snps_splicing$id))
colnames(L1_gene_list_splicing)[1] <- 'id'
nrow(L1_gene_list_splicing)

L1_chr8small_inv_splicing <- L1_gene_list_splicing %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_splicing)
L1_chr8_inv_splicing <- L1_gene_list_splicing %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_splicing)
L1_chr5_inv_splicing <- L1_gene_list_splicing %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_splicing)
L1_chr14_inv_splicing <- L1_gene_list_splicing %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_splicing)

L1_chr8small_inv_splicing <- L1_genes_w_snps_splicing %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_splicing)
L1_chr8_inv_splicing <- L1_genes_w_snps_splicing %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_splicing)
L1_chr5_inv_splicing <- L1_genes_w_snps_splicing %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_splicing)
L1_chr14_inv_splicing <- L1_genes_w_snps_splicing %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_splicing)

### S1 SNPs and Genes

#nonsynonymous SNV
S1_genes_w_snps_nonsynSNV <- S1_snps %>%  ## Only snps that fall within known S1 genes which should be the case but just in case.
  filter(grepl("MgL1_", Gene_impacted)) %>%
  filter(SNP_type == "nonsynonymous SNV") ## Filter out snps that are a specific type
nrow(S1_genes_w_snps_nonsynSNV)

S1_gene_list_nonsynSNV <- as.data.frame(unique(S1_genes_w_snps_nonsynSNV$id)) # Even with in the type of snp there can be repeats (i.e. multiple snps). Use the unique function to count the number of genes rather than the number of snps
colnames(S1_gene_list_nonsynSNV)[1] <- 'id'
nrow(S1_gene_list_nonsynSNV)

S1_chr8small_inv_nonsynSNV <- S1_gene_list_nonsynSNV %>% # Pull out genes from master list of genes with snps that fall within inversion.
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_nonsynSNV)
S1_chr8_inv_nonsynSNV <- S1_gene_list_nonsynSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_nonsynSNV)
S1_chr5_inv_nonsynSNV <- S1_gene_list_nonsynSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_nonsynSNV)
S1_chr14_inv_nonsynSNV <- S1_gene_list_nonsynSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_nonsynSNV)


S1_chr8small_inv_nonsynSNV <- S1_genes_w_snps_nonsynSNV %>% # Pull out genes from master list of genes with snps that fall within inversion.
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_nonsynSNV)
S1_chr8_inv_nonsynSNV <- S1_genes_w_snps_nonsynSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_nonsynSNV)
S1_chr5_inv_nonsynSNV <- S1_genes_w_snps_nonsynSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_nonsynSNV)
S1_chr14_inv_nonsynSNV <- S1_genes_w_snps_nonsynSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_nonsynSNV)

# Unknown
S1_genes_w_snps_unknown <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "unknown")
S1_gene_list_unknown <- as.data.frame(unique(S1_genes_w_snps_unknown$id))
colnames(S1_gene_list_unknown)[1] <- 'id'
S1_chr8small_inv_unknown <- S1_gene_list_unknown %>% 
  filter(id %in% chr_8small_genes$id)
S1_chr8_inv_unknown <- S1_gene_list_unknown %>% 
  filter(id %in% chr_8_genes$id)
S1_chr5_inv_unknown <- S1_gene_list_unknown %>% 
  filter(id %in% chr_5_genes$id)
S1_chr14_inv_unknown <- S1_gene_list_unknown %>% 
  filter(id %in% chr_14_genes$id)

#Nonframeshift insertion
S1_genes_w_snps_nonframeshift_insertion <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "nonframeshift insertion")
S1_gene_list_nonframeshift_insertion <- as.data.frame(unique(S1_genes_w_snps_nonframeshift_insertion$id)) 
colnames(S1_gene_list_nonframeshift_insertion)[1] <- 'id'
S1_chr8small_inv_nonframeshift_insertion <- S1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_8small_genes$id)
S1_chr8_inv_nonframeshift_insertion <- S1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_8_genes$id)
S1_chr5_inv_nonframeshift_insertion <- S1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_5_genes$id)
S1_chr14_inv_nonframeshift_insertion <- S1_gene_list_nonframeshift_insertion %>% 
  filter(id %in% chr_14_genes$id)

# Stop gain 
S1_genes_w_snps_stopgain <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "stopgain")
nrow(S1_genes_w_snps_stopgain)

S1_gene_list_stopgain <- as.data.frame(unique(S1_genes_w_snps_stopgain$id)) 
colnames(S1_gene_list_stopgain)[1] <- 'id'
nrow(S1_gene_list_stopgain)

S1_chr8small_inv_stopgain <- S1_gene_list_stopgain %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_stopgain)
S1_chr8_inv_stopgain <- S1_gene_list_stopgain %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_stopgain)
S1_chr5_inv_stopgain <- S1_gene_list_stopgain %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_stopgain)
S1_chr14_inv_stopgain <- S1_gene_list_stopgain %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_stopgain)

S1_chr8small_inv_stopgain <- S1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_stopgain)
S1_chr8_inv_stopgain <- S1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_stopgain)
S1_chr5_inv_stopgain <- S1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_stopgain)
S1_chr14_inv_stopgain <- S1_genes_w_snps_stopgain %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_stopgain)

#Frameshift deletion
S1_genes_w_snps_frameshift_deletion <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "frameshift deletion")
nrow(S1_genes_w_snps_frameshift_deletion)

S1_gene_list_frameshift_deletion <- as.data.frame(unique(S1_genes_w_snps_frameshift_deletion$id)) 
colnames(S1_gene_list_frameshift_deletion)[1] <- 'id'
nrow(S1_gene_list_frameshift_deletion)

S1_chr8small_inv_frameshift_deletion <- S1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_frameshift_deletion)
S1_chr8_inv_frameshift_deletion <- S1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_frameshift_deletion)
S1_chr5_inv_frameshift_deletion <- S1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_frameshift_deletion)
S1_chr14_inv_frameshift_deletion <- S1_gene_list_frameshift_deletion %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_frameshift_deletion)

S1_chr8small_inv_frameshift_deletion <- S1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_frameshift_deletion)
S1_chr8_inv_frameshift_deletion <- S1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_frameshift_deletion)
S1_chr5_inv_frameshift_deletion <- S1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_frameshift_deletion)
S1_chr14_inv_frameshift_deletion <- S1_genes_w_snps_frameshift_deletion %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_frameshift_deletion)

# Non frameshift deletion
S1_genes_w_snps_nonframeshift_deletion <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "nonframeshift deletion")
S1_gene_list_nonframeshift_deletion <- as.data.frame(unique(S1_genes_w_snps_nonframeshift_deletion$id)) 
colnames(S1_gene_list_nonframeshift_deletion)[1] <- 'id'
S1_chr8small_inv_nonframeshift_deletion <- S1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_8small_genes$id)
S1_chr8_inv_nonframeshift_deletion <- S1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_8_genes$id)
S1_chr5_inv_nonframeshift_deletion <- S1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_5_genes$id)
S1_chr14_inv_nonframeshift_deletion <- S1_gene_list_nonframeshift_deletion %>% 
  filter(id %in% chr_14_genes$id)

# Frameshift insertion
S1_genes_w_snps_frameshift_insertion <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "frameshift insertion")
nrow(S1_genes_w_snps_frameshift_insertion)

S1_gene_list_frameshift_insertion <- as.data.frame(unique(S1_genes_w_snps_frameshift_insertion$id)) 
colnames(S1_gene_list_frameshift_insertion)[1] <- 'id'
nrow(S1_gene_list_frameshift_insertion)

S1_chr8small_inv_frameshift_insertion <- S1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_frameshift_insertion)
S1_chr8_inv_frameshift_insertion <- S1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_frameshift_insertion)
S1_chr5_inv_frameshift_insertion <- S1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_frameshift_insertion)
S1_chr14_inv_frameshift_insertion <- S1_gene_list_frameshift_insertion %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_frameshift_insertion)

S1_chr8small_inv_frameshift_insertion <- S1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_frameshift_insertion)
S1_chr8_inv_frameshift_insertion <- S1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_frameshift_insertion)
S1_chr5_inv_frameshift_insertion <- S1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_frameshift_insertion)
S1_chr14_inv_frameshift_insertion <- S1_genes_w_snps_frameshift_insertion %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_frameshift_insertion)

# Stop loss
S1_genes_w_snps_stoploss <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "stoploss")
nrow(S1_genes_w_snps_stoploss)

S1_gene_list_stoploss <- as.data.frame(unique(S1_genes_w_snps_stoploss$id))
colnames(S1_gene_list_stoploss)[1] <- 'id'
nrow(S1_gene_list_stoploss)

S1_chr8small_inv_stoploss <- S1_gene_list_stoploss %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_stoploss)
S1_chr8_inv_stoploss <- S1_gene_list_stoploss %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_stoploss)
S1_chr5_inv_stoploss <- S1_gene_list_stoploss %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_stoploss)
S1_chr14_inv_stoploss <- S1_gene_list_stoploss %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_stoploss)

S1_chr8small_inv_stoploss <- S1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_stoploss)
S1_chr8_inv_stoploss <- S1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_stoploss)
S1_chr5_inv_stoploss <- S1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_stoploss)
S1_chr14_inv_stoploss <- S1_genes_w_snps_stoploss %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_stoploss)

#Intronic
S1_genes_w_snps_intronic <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "intronic")
S1_gene_list_intronic <- as.data.frame(unique(S1_genes_w_snps_intronic$id))
colnames(S1_gene_list_intronic)[1] <- 'id'
S1_chr8small_inv_intronic <- S1_gene_list_intronic %>% 
  filter(id %in% chr_8small_genes$id)
S1_chr8_inv_intronic <- S1_gene_list_intronic %>% 
  filter(id %in% chr_8_genes$id)
S1_chr5_inv_intronic <- S1_gene_list_intronic %>% 
  filter(id %in% chr_5_genes$id)
S1_chr14_inv_intronic <- S1_gene_list_intronic %>% 
  filter(id %in% chr_14_genes$id)

# UTR3
S1_genes_w_snps_UTR3 <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "UTR3")
S1_gene_list_UTR3 <- as.data.frame(unique(S1_genes_w_snps_UTR3$id))
colnames(S1_gene_list_UTR3)[1] <- 'id'
S1_chr8small_inv_UTR3 <- S1_gene_list_UTR3 %>% 
  filter(id %in% chr_8small_genes$id)
S1_chr8_inv_UTR3 <- S1_gene_list_UTR3 %>% 
  filter(id %in% chr_8_genes$id)
S1_chr5_inv_UTR3 <- S1_gene_list_UTR3 %>% 
  filter(id %in% chr_5_genes$id)
S1_chr14_inv_UTR3 <- S1_gene_list_UTR3 %>% 
  filter(id %in% chr_14_genes$id)

#UTR5
S1_genes_w_snps_UTR5 <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "UTR5")
S1_gene_list_UTR5 <- as.data.frame(unique(S1_genes_w_snps_UTR5$id))
colnames(S1_gene_list_UTR5)[1] <- 'id'
S1_chr8small_inv_UTR5 <- S1_gene_list_UTR5 %>% 
  filter(id %in% chr_8small_genes$id)
S1_chr8_inv_UTR5 <- S1_gene_list_UTR5 %>% 
  filter(id %in% chr_8_genes$id)
S1_chr5_inv_UTR5 <- S1_gene_list_UTR5 %>% 
  filter(id %in% chr_5_genes$id)
S1_chr14_inv_UTR5 <- S1_gene_list_UTR5 %>% 
  filter(id %in% chr_14_genes$id)

#Splicing
S1_genes_w_snps_splicing <- S1_snps %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "splicing")
nrow(S1_genes_w_snps_splicing)

S1_gene_list_splicing <- as.data.frame(unique(S1_genes_w_snps_splicing$id))
colnames(S1_gene_list_splicing)[1] <- 'id'
nrow(S1_gene_list_splicing)

S1_chr8small_inv_splicing <- S1_gene_list_splicing %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_splicing)
S1_chr8_inv_splicing <- S1_gene_list_splicing %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_splicing)
S1_chr5_inv_splicing <- S1_gene_list_splicing %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_splicing)
S1_chr14_inv_splicing <- S1_gene_list_splicing %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_splicing)

S1_chr8small_inv_splicing <- S1_genes_w_snps_splicing %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_splicing)
S1_chr8_inv_splicing <- S1_genes_w_snps_splicing %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_splicing)
S1_chr5_inv_splicing <- S1_genes_w_snps_splicing %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_splicing)
S1_chr14_inv_splicing <- S1_genes_w_snps_splicing %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_splicing)

# Nonsynonymous snps (only in annovar output genotype_exonic_funtion)
synv_annotated_snps_S1 <- read.table("data/high-impact-S1_synSNV_nothets.csv", sep = ",")
colnames(synv_annotated_snps_S1) <- c("SNP_type","Gene_impacted","Zygosity")
synv_annotated_snps_S1$id <- gsub(":Mg.*", "", synv_annotated_snps_S1$Gene_impacted)
synv_annotated_snps_S1$id <- gsub("\\(Mg.*", "", synv_annotated_snps_S1$id)
synv_annotated_snps_S1_unique <- as.data.frame(unique(synv_annotated_snps_S1$id))

S1_genes_w_snps_synSNV <- synv_annotated_snps_S1 %>%  ## Only snps that fall within known S1 genes
  filter(grepl("MgL1_", Gene_impacted)) %>% 
  filter(SNP_type == "synonymous SNV")
nrow(S1_genes_w_snps_synSNV)

S1_gene_list_synSNV <- as.data.frame(unique(S1_genes_w_snps_synSNV$id))
colnames(S1_gene_list_synSNV)[1] <- 'id'
nrow(S1_gene_list_synSNV)

S1_chr8small_inv_synSNV <- S1_gene_list_synSNV %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_synSNV)
S1_chr8_inv_synSNV <- S1_gene_list_synSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_synSNV)
S1_chr5_inv_synSNV <- S1_gene_list_synSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_synSNV)
S1_chr14_inv_synSNV <- S1_gene_list_synSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_synSNV)

S1_chr8small_inv_synSNV <- S1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_8small_genes$id)
nrow(S1_chr8small_inv_synSNV)
S1_chr8_inv_synSNV <- S1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(S1_chr8_inv_synSNV)
S1_chr5_inv_synSNV <- S1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(S1_chr5_inv_synSNV)
S1_chr14_inv_synSNV <- S1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(S1_chr14_inv_synSNV)

synv_annotated_snps_L1 <- read.table("data/high-impact-L1_synSNV_nothets.csv", sep = ",")
colnames(synv_annotated_snps_L1) <- c("SNP_type","Gene_impacted","Zygosity")
synv_annotated_snps_L1$id <- gsub(":Mg.*", "", synv_annotated_snps_L1$Gene_impacted)
synv_annotated_snps_L1$id <- gsub("\\(Mg.*", "", synv_annotated_snps_L1$id)
synv_annotated_snps_L1_unique <- as.data.frame(unique(synv_annotated_snps_L1$id))

L1_genes_w_snps_synSNV <- synv_annotated_snps_L1 %>%  ## Only snps that fall within known L1 genes
  filter(grepl("MgS1_", Gene_impacted)) %>% 
  filter(SNP_type == "synonymous SNV")
nrow(L1_genes_w_snps_synSNV)

L1_gene_list_synSNV <- as.data.frame(unique(L1_genes_w_snps_synSNV$id))
colnames(L1_gene_list_synSNV)[1] <- 'id'
nrow(L1_gene_list_synSNV)

L1_chr8small_inv_synSNV <- L1_gene_list_synSNV %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_synSNV)
L1_chr8_inv_synSNV <- L1_gene_list_synSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_synSNV)
L1_chr5_inv_synSNV <- L1_gene_list_synSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_synSNV)
L1_chr14_inv_synSNV <- L1_gene_list_synSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_synSNV)

L1_chr8small_inv_synSNV <- L1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_8small_genes$id)
nrow(L1_chr8small_inv_synSNV)
L1_chr8_inv_synSNV <- L1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_8_genes$id)
nrow(L1_chr8_inv_synSNV)
L1_chr5_inv_synSNV <- L1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_5_genes$id)
nrow(L1_chr5_inv_synSNV)
L1_chr14_inv_synSNV <- L1_genes_w_snps_synSNV %>% 
  filter(id %in% chr_14_genes$id)
nrow(L1_chr14_inv_synSNV)
