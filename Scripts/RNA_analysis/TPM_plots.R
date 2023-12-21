# Set working directory and load necessary packages
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(data.table)
library(ggpubr)

setwd("/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/")

####################### TPM ####################

######
# L1 #
######

L1_length <- read.table("./data/DEG_L1/L1_transcript_length.tsv")
colnames(L1_length) <- c("gene_model", "length")

abundance <- read.table("./data/DEG_L1/L1_star_2pass.tsv", header = T)
names(abundance)[names(abundance) == 'Genes'] <- 'gene_model'

L1_length_abundance <- L1_length %>% 
  inner_join(abundance, by ="gene_model", multiple = "all") # Adding L1 outlier data to dataset of syntenic only genes


L1_length_abundance_long <- data.frame(t(as.matrix(L1_length_abundance)))
colnames(L1_length_abundance_long) <- L1_length_abundance_long[1,]
L1_length_abundance_long <- L1_length_abundance_long[2:52,]

L1_length_abundance_standardized <- L1_length_abundance %>%
  mutate(
    across(c(3:52),
           .fns = ~./length))

L1_length_abundance_standardized <- L1_length_abundance_standardized %>% mutate_if(is.integer, as.numeric)
L1_length_abundance_standardized_onlynumeric <- L1_length_abundance_standardized[,3:52]

L1_TPM <- as.data.frame(apply(L1_length_abundance_standardized_onlynumeric,2,function(x){x/sum(x)*100000000}))

L1_TPM_final <- cbind(L1_length_abundance$gene_model, L1_TPM)
names(L1_TPM_final)[names(L1_TPM_final) == "L1_length_abundance$gene_model"] <- 'gene_model'

######
# S1 #
######

S1_length <- read.table("./data/DEG_S1/S1_transcript_length.tsv")
colnames(S1_length) <- c("gene_model", "length")

abundance <- read.table("./data/DEG_S1/S1_star_2pass.tsv", header = T)
names(abundance)[names(abundance) == 'Genes'] <- 'gene_model'

S1_length_abundance <- S1_length %>% 
  inner_join(abundance, by ="gene_model", multiple = "all") # Adding S1 outlier data to dataset of syntenic only genes

S1_length_abundance_long <- data.frame(t(as.matrix(S1_length_abundance)))
colnames(S1_length_abundance_long) <- S1_length_abundance_long[1,]
S1_length_abundance_long <- S1_length_abundance_long[2:52,]

S1_length_abundance_standardized <- S1_length_abundance %>%
  mutate(
    across(c(3:52),
           .fns = ~./length))

S1_length_abundance_standardized <- S1_length_abundance_standardized %>% mutate_if(is.integer, as.numeric)
S1_length_abundance_standardized_onlynumeric <- S1_length_abundance_standardized[,3:52]

S1_TPM <- as.data.frame(apply(S1_length_abundance_standardized_onlynumeric,2,function(x){x/sum(x)*100000000}))

S1_TPM_final <- cbind(S1_length_abundance$gene_model, S1_TPM)
names(S1_TPM_final)[names(S1_TPM_final) == "S1_length_abundance$gene_model"] <- 'gene_model'


#### Combining Data Sets ####

L1_S1_TPM <- rbind(S1_TPM_final, L1_TPM_final)

### Plotting DEG at the inversion breakpionts ###
L1_S1_TPM_long <- data.frame(t(as.matrix(L1_S1_TPM)))
colnames(L1_S1_TPM_long) <- L1_S1_TPM_long[1,]
L1_S1_TPM_long <- L1_S1_TPM_long[2:52,]
variable <- data.frame(rownames(L1_S1_TPM_long))
variable$rownames.L1_S1_TPM_long. <- gsub("*_[0-9]+", " ", variable$rownames.L1_S1_TPM_long.)
colnames(variable) <- "condition"
L1_S1_TPM_long_final <- cbind(variable,L1_S1_TPM_long )
#L1_S1_TPM_long_final$condition <- as.factor(L1_S1_TPM_long_final$condition)
L1_S1_TPM_long_final <- L1_S1_TPM_long_final %>% select(-condition) %>% mutate_if(is.character, as.numeric)
L1_S1_TPM_long_final_2 <- cbind(variable,L1_S1_TPM_long_final )
L1_S1_TPM_long_final_3 <- na.omit(L1_S1_TPM_long_final_2)


L1_S1_TPM_long_final_3$condition <- gsub("L1_Coastal", "LMC-L1 Coastal", L1_S1_TPM_long_final_3$condition)
L1_S1_TPM_long_final_3$condition <- gsub("L1_Inland", "LMC-L1 Inland", L1_S1_TPM_long_final_3$condition)
L1_S1_TPM_long_final_3$condition <- gsub("S1_Coastal", "SWB-S1 Coastal", L1_S1_TPM_long_final_3$condition)
L1_S1_TPM_long_final_3$condition <- gsub("S1_Inland", "SWB-S1 Inland", L1_S1_TPM_long_final_3$condition)

## DEG at genes at the inversion breakpoint ##

MgL1_05g16120 <- ggplot(L1_S1_TPM_long_final_3,aes(y= MgL1_05g16120, x = condition, fill=condition))+
geom_boxplot()+
labs(x="Treatment", y = "MgL1_05g16120 expression (TPM)")+
theme_minimal()+
scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())

MgS1_05g10260 <- ggplot(L1_S1_TPM_long_final_3,aes(y= MgS1_05g10260, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = "MgS1_05g10260 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())




MgL1_08g01750 <- ggplot(L1_S1_TPM_long_final_3,aes(y=  MgL1_08g01750, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgL1_08g01750 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())

MgS1_08g08930 <-  ggplot(L1_S1_TPM_long_final_3,aes(y= MgS1_08g08930, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgS1_08g08930 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())

MgL1_08g01750<-  ggplot(L1_S1_TPM_long_final_3,aes(y= MgL1_08g01750, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgL1_08g01750 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())

MgL1_14g10910 <-  ggplot(L1_S1_TPM_long_final_3,aes(y= MgL1_14g10910, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgL1_14g10910 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())

MgS1_14g08070 <-  ggplot(L1_S1_TPM_long_final_3,aes(y= MgS1_14g08070, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgS1_14g08070 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())

MgS1_14g08050 <-  ggplot(L1_S1_TPM_long_final_3,aes(y= MgS1_14g08050, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgS1_14g08050 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.title.x=element_blank())


MgS1_14g08070 <-  ggplot(L1_S1_TPM_long_final_3,aes(y= MgS1_14g08070, x = condition, fill=condition))+
  geom_boxplot()+
  labs(x="Treatment", y = " MgS1_14g08070 expression (TPM)")+
  theme_minimal()+
  scale_fill_manual( values = c("#D95F02","#1B9E77","#7570B3","#E7298A"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"),axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x=element_blank())


bps_DEG <- ggarrange(MgS1_05g10260,MgL1_05g16120,  MgS1_08g08930, MgL1_08g01750,  MgL1_14g10910, MgS1_14g08070, MgS1_14g08050, common.legend = T )
 