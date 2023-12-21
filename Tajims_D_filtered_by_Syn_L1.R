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

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned"
setwd(path)

# Files for Tajimas D
files=list.files(pattern = ".tajD")
L1_aligned_dataset_tajD = do.call(rbind, lapply(files,fread))

# Removed NAs for tajimas D
L1_aligned_dataset_nona_tajD <- na.omit(L1_aligned_dataset_tajD)
colnames(L1_aligned_dataset_nona_tajD) <- c("gene", "scaffold", "gene_start", "gene_end","gene_len", "IA_gene_cov", "IA_S", "TajD_IA", "CP_gene_cov", "CP_S", "TajD_CP", "tot_gene_cov", "tot_S", "TajD_tot")

# Subsetting the data based on Syntenic genes

L1_syn_genes <- read.table("L1_syntelogs.txt")

L1_aligned_dataset_nona_tajD_synonly <- subset(L1_aligned_dataset_nona_tajD, gene %in% L1_syn_genes$V1)


######### barplots of Tajimas D ##########

L1_aligned_dataset_nona_tajD$scaffold <- gsub("chr","", as.character(L1_aligned_dataset_nona_tajD$scaffold))
L1_aligned_dataset_nona_tajD$scaffold <- as.numeric(L1_aligned_dataset_nona_tajD$scaffold)

data_cum_L1_D <- L1_aligned_dataset_nona_tajD %>% 
  group_by(scaffold) %>% 
  summarise(max_bp=max(gene_end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaffold, bp_add)

L1_aligned_dataset_nona_tajD <- L1_aligned_dataset_nona_tajD %>% 
  inner_join(data_cum_L1_D, by ="scaffold") %>% 
  mutate(bp_cum = gene_end + bp_add)

# Adding inverted vs noninverted to the chromosomes of interest.
L1_aligned_dataset_nona_tajD_5_NOTinv1 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum < 67446946+13650670) 
L1_aligned_dataset_nona_tajD_5_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_5_NOTinv2 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum > 67446946+17847181)
L1_aligned_dataset_nona_tajD_5_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_5_NotINV_fig <- rbind(L1_aligned_dataset_nona_tajD_5_NOTinv2, L1_aligned_dataset_nona_tajD_5_NOTinv1)

L1_aligned_dataset_nona_tajD_5_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(between(bp_cum, (67446946+13650670),(67446946+17847181) ))
L1_aligned_dataset_nona_tajD_5_inv$karyo <- "inv 5"

L1_aligned_dataset_nona_tajD_5 <- rbind(L1_aligned_dataset_nona_tajD_5_inv, L1_aligned_dataset_nona_tajD_5_NotINV_fig)

################################################################################

L1_aligned_dataset_nona_tajD_8_NOTinv1 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118735039+850429 )) 
L1_aligned_dataset_nona_tajD_8_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_8_NOTinv2 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118735039+7604769))
L1_aligned_dataset_nona_tajD_8_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_8_NotINV_fig <- rbind(L1_aligned_dataset_nona_tajD_8_NOTinv1, L1_aligned_dataset_nona_tajD_8_NOTinv2)

L1_aligned_dataset_nona_tajD_8_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118735039+850429),(118735039+7604769)))
L1_aligned_dataset_nona_tajD_8_inv$karyo <- "inv 8"

#removing small inversion from 8 inversion
L1_aligned_dataset_nona_tajD_8_inv_wosmalL1 <-  L1_aligned_dataset_nona_tajD_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118735039+1032334)) 

L1_aligned_dataset_nona_tajD_8_inv_wosmall2 <-  L1_aligned_dataset_nona_tajD_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118735039+1246126))

L1_aligned_dataset_nona_tajD_8 <- rbind(L1_aligned_dataset_nona_tajD_8_inv_wosmalL1, L1_aligned_dataset_nona_tajD_8_inv_wosmall2, L1_aligned_dataset_nona_tajD_8_NotINV_fig)

################################################################################

L1_aligned_dataset_nona_tajD_8small_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118735039+1032334), (118735039+1246126)))
L1_aligned_dataset_nona_tajD_8small_inv$karyo <- "inv 8 SMALL"

############################################################################################################################

L1_aligned_dataset_nona_tajD_14_NOTinv1 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum < (232784694+5329939)) 

L1_aligned_dataset_nona_tajD_14_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_14_NOTinv2 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum > (232784694+7791197))
L1_aligned_dataset_nona_tajD_14_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_14_NotINV_fig <- rbind(L1_aligned_dataset_nona_tajD_14_NOTinv2, L1_aligned_dataset_nona_tajD_14_NOTinv1)

L1_aligned_dataset_nona_tajD_14_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(between(bp_cum, (232784694+5329939), (232784694+7791197)))
L1_aligned_dataset_nona_tajD_14_inv$karyo <- "inv 14"

L1_aligned_dataset_nona_tajD_14 <- rbind(L1_aligned_dataset_nona_tajD_14_inv, L1_aligned_dataset_nona_tajD_14_NotINV_fig)

#######################################################

# TajD dataset
### Inverted genes didnt pass?
tajd <- rbind(L1_aligned_dataset_nona_tajD_5, L1_aligned_dataset_nona_tajD_8,L1_aligned_dataset_nona_tajD_8small_inv,L1_aligned_dataset_nona_tajD_14)

L1_aligned_dataset_nona_tajD_noninv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold %in% c("1","2","3","4","6","7","9","10","11","12","13"))

L1_aligned_dataset_nona_tajD_noninv$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_final <- rbind(L1_aligned_dataset_nona_tajD_noninv,tajd)

data1 <- L1_aligned_dataset_nona_tajD_final[,c(2,8,17)]
colnames(data1) <- c("chromosome", "TajD", "Karyo")
data1$tajtype <- "IA"


data2 <- L1_aligned_dataset_nona_tajD_final[,c(2,11,17)]
colnames(data2) <- c("chromosome", "TajD", "Karyo")
data2$tajtype <- "CP"

data3 <- L1_aligned_dataset_nona_tajD_final[,c(2,14,17)]
colnames(data3) <- c("chromosome", "TajD", "Karyo")
data3$tajtype <- "total"

loath <- rbind(data1, data2, data3)

tajD_pot_L1 <- loath %>% 
  ggplot(aes(fill=tajtype, y= TajD, x= Karyo)) +
  geom_boxplot() +
  xlab("Genomic location")+
  ylab("Tajima's D")+
  scale_fill_manual(name = "Ecotype",labels = c("Coastal Perennial", "Inland Annual", "Total"),values = c("#D95F02","#1B9E77","#7570B3"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))


png("./figures/tajD_pot_L1_synonly.png", width = 15, height = 15, units = "cm",res = 300)
print(tajD_pot_L1)
dev.off()
# 