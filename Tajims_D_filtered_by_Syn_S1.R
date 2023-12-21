# Subsetting Tajimas D data for only syntenic genes
# Leslie M Kollar

## Loading in libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggtext)

## Loading in data files

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/S1_aligned"
setwd(path)

# Files for Tajimas D
files=list.files(pattern = ".tajD")
S1_aligned_dataset_tajD = do.call(rbind, lapply(files,fread))

# Removed NAs for tajimas D
S1_aligned_dataset_nona_tajD <- na.omit(S1_aligned_dataset_tajD)
colnames(S1_aligned_dataset_nona_tajD) <- c("gene", "scaffold", "gene_start", "gene_end","gene_len", "IA_gene_cov", "IA_S", "TajD_IA", "CP_gene_cov", "CP_S", "TajD_CP", "tot_gene_cov", "tot_S", "TajD_tot")

# Subsetting the data based on Syntenic genes

S1_syn_genes <- read.table("S1_syntelogs.txt")

S1_aligned_dataset_nona_tajD_synonly <- subset(S1_aligned_dataset_nona_tajD, gene %in% S1_syn_genes$V1)


######### barplots of Tajimas D ##########

S1_aligned_dataset_nona_tajD_synonly$scaffold <- gsub("chr","", as.character(S1_aligned_dataset_nona_tajD_synonly$scaffold))
S1_aligned_dataset_nona_tajD_synonly$scaffold <- as.numeric(S1_aligned_dataset_nona_tajD_synonly$scaffold)

data_cum_S1_D <- S1_aligned_dataset_nona_tajD_synonly %>% 
  group_by(scaffold) %>% 
  summarise(max_bp=max(gene_end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaffold, bp_add)

S1_aligned_dataset_nona_tajD_synonly <- S1_aligned_dataset_nona_tajD_synonly %>% 
  inner_join(data_cum_S1_D, by ="scaffold") %>% 
  mutate(bp_cum = gene_end + bp_add)

# Adding inverted vs noninverted to the chromosomes of interest.
S1_aligned_dataset_nona_tajD_synonly_5_NOTinv1 <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum < 66530930+14125309) 
S1_aligned_dataset_nona_tajD_synonly_5_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_5_NOTinv2 <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum > 66530930+18136144)
S1_aligned_dataset_nona_tajD_synonly_5_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_5_NotINV_fig <- rbind(S1_aligned_dataset_nona_tajD_synonly_5_NOTinv2, S1_aligned_dataset_nona_tajD_synonly_5_NOTinv1)

S1_aligned_dataset_nona_tajD_synonly_5_inv <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="5") %>% 
  filter(between(bp_cum, (66530930+14125309),(66530930+18136144) ))
S1_aligned_dataset_nona_tajD_synonly_5_inv$karyo <- "inv 5"

S1_aligned_dataset_nona_tajD_synonly_5 <- rbind(S1_aligned_dataset_nona_tajD_synonly_5_inv, S1_aligned_dataset_nona_tajD_synonly_5_NotINV_fig)

################################################################################

S1_aligned_dataset_nona_tajD_synonly_8_NOTinv1 <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118351995+858736)) 
S1_aligned_dataset_nona_tajD_synonly_8_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_8_NOTinv2 <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118351995+6465310))
S1_aligned_dataset_nona_tajD_synonly_8_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_8_NotINV_fig <- rbind(S1_aligned_dataset_nona_tajD_synonly_8_NOTinv1, S1_aligned_dataset_nona_tajD_synonly_8_NOTinv2)

S1_aligned_dataset_nona_tajD_synonly_8_inv <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118351995+858736),(118351995+6465310)))
S1_aligned_dataset_nona_tajD_synonly_8_inv$karyo <- "inv 8"

#removing small inversion from 8 inversion
S1_aligned_dataset_nona_tajD_synonly_8_inv_wosmalS1 <-  S1_aligned_dataset_nona_tajD_synonly_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118351995+5998986)) 

S1_aligned_dataset_nona_tajD_synonly_8_inv_wosmall2 <-  S1_aligned_dataset_nona_tajD_synonly_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118351995+6260288))

S1_aligned_dataset_nona_tajD_synonly_8 <- rbind(S1_aligned_dataset_nona_tajD_synonly_8_inv_wosmalS1, S1_aligned_dataset_nona_tajD_synonly_8_inv_wosmall2, S1_aligned_dataset_nona_tajD_synonly_8_NotINV_fig)

################################################################################

S1_aligned_dataset_nona_tajD_synonly_8small_inv <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118351995+5998986), (118351995+6260288)))
S1_aligned_dataset_nona_tajD_synonly_8small_inv$karyo <- "inv 8 SMALL"

############################################################################################################################

S1_aligned_dataset_nona_tajD_synonly_14_NOTinv1 <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum < (231994157+5892556)) 

S1_aligned_dataset_nona_tajD_synonly_14_NOTinv1$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_14_NOTinv2 <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum > (231994157+8677481))
S1_aligned_dataset_nona_tajD_synonly_14_NOTinv2$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_14_NotINV_fig <- rbind(S1_aligned_dataset_nona_tajD_synonly_14_NOTinv2, S1_aligned_dataset_nona_tajD_synonly_14_NOTinv1)

S1_aligned_dataset_nona_tajD_synonly_14_inv <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold=="14") %>% 
  filter(between(bp_cum, (231994157+5892556), (231994157+8677481)))
S1_aligned_dataset_nona_tajD_synonly_14_inv$karyo <- "inv 14"

S1_aligned_dataset_nona_tajD_synonly_14 <- rbind(S1_aligned_dataset_nona_tajD_synonly_14_inv, S1_aligned_dataset_nona_tajD_synonly_14_NotINV_fig)

#######################################################

# TajD dataset
### Inverted genes didnt pass?
tajd <- rbind(S1_aligned_dataset_nona_tajD_synonly_5, S1_aligned_dataset_nona_tajD_synonly_8,S1_aligned_dataset_nona_tajD_synonly_8small_inv,S1_aligned_dataset_nona_tajD_synonly_14)

S1_aligned_dataset_nona_tajD_synonly_noninv <-  S1_aligned_dataset_nona_tajD_synonly %>% 
  filter(scaffold %in% c("1","2","3","4","6","7","9","10","11","12","13"))

S1_aligned_dataset_nona_tajD_synonly_noninv$karyo <- "notinv"

S1_aligned_dataset_nona_tajD_synonly_final <- rbind(S1_aligned_dataset_nona_tajD_synonly_noninv,tajd)

data1 <- S1_aligned_dataset_nona_tajD_synonly_final[,c(2,8,17)]
colnames(data1) <- c("chromosome", "TajD", "Karyo")
data1$tajtype <- "IA"


data2 <- S1_aligned_dataset_nona_tajD_synonly_final[,c(2,11,17)]
colnames(data2) <- c("chromosome", "TajD", "Karyo")
data2$tajtype <- "CP"

data3 <- S1_aligned_dataset_nona_tajD_synonly_final[,c(2,14,17)]
colnames(data3) <- c("chromosome", "TajD", "Karyo")
data3$tajtype <- "total"

loath <- rbind(data1, data2, data3)

tajD_pot_S1 <- loath %>% 
  ggplot(aes(fill=tajtype, y= TajD, x= Karyo)) +
  geom_boxplot() +
  xlab("Genomic location")+
  ylab("Tajima's D")+
  scale_fill_manual(name = "Ecotype",labels = c("Coastal Perennial", "Inland Annual", "Total"),values = c("#D95F02","#1B9E77","#7570B3"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))


png("./figures/tajD_pot_S1_synonly.png", width = 15, height = 15, units = "cm",res = 300)
print(tajD_pot_S1)
dev.off()
