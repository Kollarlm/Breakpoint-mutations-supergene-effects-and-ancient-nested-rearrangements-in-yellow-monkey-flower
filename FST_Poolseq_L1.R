## # Plotting output for F python scripts from Gould et al 2017 
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

# Files for pi
files=list.files(pattern = ".fst")
L1_aligned_dataset_fst = do.call(rbind, lapply(files,fread))

# Removing windowns with 2 SD difference from mean

L1_SD <- L1_aligned_dataset_fst %>% 
  summarise(win_depth_sd=sd(win_depth, na.rm=TRUE)) # standard deviation = 227.9529
L1_filt <- (2*L1_SD)

L1_aligned_dataset_fst_filtered <- L1_aligned_dataset_fst %>% # removes data greater than 2 SD from mean
  filter(abs(scale(win_depth)) < 2)

#Taking the top 1% of fst

obs_L1 <-  nrow(L1_aligned_dataset_fst_filtered)

L1_aligned_dataset_fst_filtered <- L1_aligned_dataset_fst_filtered %>% arrange(desc(win_Fst))
L1_fst_top1 <- L1_aligned_dataset_fst_filtered %>% 
  filter(row_number() < obs_L1 * 0.01)

#Changing chromosome so they are in order

L1_aligned_dataset_fst_filtered$scaff <- gsub("chr","", as.character(L1_aligned_dataset_fst_filtered$scaff))
L1_aligned_dataset_fst_filtered$scaff <- as.numeric(L1_aligned_dataset_fst_filtered$scaff)

data_cum_L1_fst <- L1_aligned_dataset_fst_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win.end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

L1_aligned_dataset_fst_filtered <- L1_aligned_dataset_fst_filtered %>% 
  inner_join(data_cum_L1_fst, by ="scaff") %>% 
  mutate(bp_cum = win.end + bp_add)

# Scatter plot for FST

L1_fst_plot <- L1_aligned_dataset_fst_filtered %>% 
  ggplot(aes(x=bp_cum, y= win_Fst)) +
  geom_point() +
  geom_hline(yintercept =0.01738773, color="black", linetype = 2) +
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "#7570B3", 
           fill = "#7570B3", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(232813426+5329939), xmax=(232813426+7791197), ymin=-Inf, ymax=Inf,
           color = "#E7298A", 
           fill = "#E7298A", alpha = 0.5, linetype = 2)+
  ylab(paste("FST"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))
  
png("./figures/L1_fst_plot.png", width = 15, height = 12, units = "cm",res = 300)
  print(L1_fst_plot)
  dev.off()
