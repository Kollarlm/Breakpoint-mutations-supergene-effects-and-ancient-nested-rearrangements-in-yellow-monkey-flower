# Plotting output pf python scripts from Gould et al 2017 
# Leslie M Kollar

## Loading in libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(ggtext)
library(knitr)

## Loading in data files

path <- "/Users/lesliekollar/Desktop/mimulus-gould2018-reanalysis/data/poolseq_output/L1_aligned"
setwd(path)

# Files for pi
# We use the fst files because windows were properly filtered.
files=list.files(pattern = ".fst")
L1_aligned_dataset_pi = do.call(rbind, lapply(files,fread))


# Files for gstat
files=list.files(pattern = ".gstat")
L1_aligned_dataset_gstat = do.call(rbind, lapply(files,fread))

# Files for Tajimas D
files=list.files(pattern = ".tajD")
L1_aligned_dataset_tajD = do.call(rbind, lapply(files,fread))

# Removed NAs for tajimas D
L1_aligned_dataset_nona_tajD <- na.omit(L1_aligned_dataset_tajD)
colnames(L1_aligned_dataset_nona_tajD) <- c("gene", "scaffold", "gene_start", "gene_end","gene_len", "IA_gene_cov", "IA_S", "TajD_IA", "CP_gene_cov", "CP_S", "TajD_CP", "tot_gene_cov", "tot_S", "TajD_tot")


## Removed windows with average depth greater than 2 standard deviations
L1_SD <- L1_aligned_dataset_pi %>% 
  summarise(win_depth_sd=sd(win_depth, na.rm=TRUE)) # standard deviation = 128.0278
L1_filt <- (2*L1_SD)

L1_mean <- L1_aligned_dataset_pi %>% 
  summarise(win_depth_sd=mean(win_depth, na.rm=TRUE)) # mean = 246.5186

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi %>% # removes data greater than 2 SD from mean
  filter(abs(scale(win_depth)) < 2)

# gstat
L1_SD_g <- L1_aligned_dataset_gstat %>% 
  summarise(win_depth_sd=sd(avg.win.depth, na.rm=TRUE)) # standard deviation = 179.8778
L1_filt_g <- (2*L1_SD_g)

L1_mean_g <- L1_aligned_dataset_gstat %>% 
  summarise(win_depth_sd=mean(avg.win.depth, na.rm=TRUE)) # mean = 

L1_aligned_dataset_gstat_filtered <- L1_aligned_dataset_gstat %>% # removes data greater than 2 SD from mean
  filter(abs(scale(avg.win.depth)) < 2)

## Adding in column for the ratio of pi

L1_aligned_dataset_pi_filtered$piRatio <- (L1_aligned_dataset_pi_filtered$piIA/L1_aligned_dataset_pi_filtered$piCP)

## Take top and bottom 1% of piRatio
obs_L1 <-  nrow(L1_aligned_dataset_pi_filtered)

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi_filtered %>% arrange(desc(piRatio))
L1_pi_top1 <- L1_aligned_dataset_pi_filtered %>% 
  filter(row_number() < obs_L1 * 0.01)


L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi_filtered %>% arrange((piRatio))
L1_pi_bottom1 <- L1_aligned_dataset_pi_filtered %>% 
  filter(row_number() < obs_L1 * 0.01)

## Take top 1% of gstat

obs_L1_gstat <-  nrow(L1_aligned_dataset_gstat_filtered)

L1_aligned_dataset_gstat_filtered <- L1_aligned_dataset_gstat_filtered %>% arrange(desc(G_stat))
L1_pi_top1_gstat <- L1_aligned_dataset_gstat_filtered %>% 
  filter(row_number() < obs_L1_gstat * 0.01)

## Preparing the data to plot manhattan plot in ggplot

# Remvoing Chr so that the chromsomes are in order

L1_aligned_dataset_pi_filtered$scaff <- gsub("chr","", as.character(L1_aligned_dataset_pi_filtered$scaff))
L1_aligned_dataset_pi_filtered$scaff <- as.numeric(L1_aligned_dataset_pi_filtered$scaff)

L1_aligned_dataset_gstat_filtered$scaff <- gsub("chr","", as.character(L1_aligned_dataset_gstat_filtered$scaff))
L1_aligned_dataset_gstat_filtered$scaff <- as.numeric(L1_aligned_dataset_gstat_filtered$scaff)

data_cum_L1_pi <- L1_aligned_dataset_pi_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win.end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

L1_aligned_dataset_pi_filtered <- L1_aligned_dataset_pi_filtered %>% 
  inner_join(data_cum_L1_pi, by ="scaff") %>% 
  mutate(bp_cum = win.end + bp_add)
  
data_cum_L1_gstat <- L1_aligned_dataset_gstat_filtered %>% 
  group_by(scaff) %>% 
  summarise(max_bp=max(win_end)) %>% 
  mutate(bp_add=lag(cumsum(max_bp), default = 0)) %>% 
  select(scaff, bp_add)

L1_aligned_dataset_gstat_filtered <- L1_aligned_dataset_gstat_filtered %>% 
  inner_join(data_cum_L1_gstat, by ="scaff") %>% 
  mutate(bp_cum=win_end + bp_add)

## Scatter plot of piRatio
## "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
pi_data <- subset(L1_aligned_dataset_pi_filtered, win.start==3586980)
  

L1_pi_plot <- ggplot(L1_aligned_dataset_pi_filtered, aes(x=bp_cum, y= log10(piRatio))) +
  geom_point() +
  geom_hline(yintercept =log10(2.567310), color="black", linetype=2) +
  geom_hline(yintercept = log10(0.5397863), color="black", linetype=2) +
  annotate(geom = "rect", xmin=(118943830+850429), xmax=(118943830+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype =2)+
  annotate(geom = "rect", xmin=(118943830+1032334), xmax=(118943830+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype=2)+
  annotate(geom = "rect", xmin=(67485786+13650670), xmax=(67485786+17847181), ymin=-Inf, ymax=Inf,
           color = "#7570B3", 
           fill = "#7570B3", alpha = 0.5, linetype=2)+
  annotate(geom = "rect", xmin=(232897274+5329939), xmax=(232897274+7791197), ymin=-Inf, ymax=Inf,
           color = "#E7298A", 
           fill = "#E7298A", alpha = 0.5, linetype=2)+
 # annotate(geom = "pointrange", x=124949634, y=log10(0.6037378), ymin = log10(0.6037378), ymax = log10(0.6037378),  color= "red",size=.5)+
#  annotate(geom = "pointrange", x=122530809, y=log10(0.5685423), ymin = log10(0.5685423), ymax = log10(0.5685423),  color= "blue",size=.5)+
 # annotate(geom = "pointrange", x=122538343, y=log10(0.6485823), ymin = log10(0.6485823), ymax = log10(0.6485823),  color= "blue",size=.5)+
  #annotate(geom = "pointrange", x=122610093, y=log10(0.8409326), ymin = log10(0.8409326), ymax = log10(0.8409326),  color= "blue",size=.5)+
  ylab(paste("log10(","\u03c0","ratio)"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

L1_dxy_plot <- ggplot(L1_aligned_dataset_pi_filtered, aes(x=bp_cum, y= Dxy_win)) +
  geom_point() +
  annotate(geom = "rect", xmin=(118943830+850429), xmax=(118943830+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype =2)+
  annotate(geom = "rect", xmin=(118943830+1032334), xmax=(118943830+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype=2)+
  annotate(geom = "rect", xmin=(67485786+13650670), xmax=(67485786+17847181), ymin=-Inf, ymax=Inf,
           color = "#7570B3", 
           fill = "#7570B3", alpha = 0.5, linetype=2)+
  annotate(geom = "rect", xmin=(232897274+5329939), xmax=(232897274+7791197), ymin=-Inf, ymax=Inf,
           color = "#E7298A", 
           fill = "#E7298A", alpha = 0.5, linetype=2)+
ylab(paste("Dxy"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))


# Scatter plot of read depth
L1_depth_plot <- ggplot(L1_aligned_dataset_gstat_filtered, aes(x=bp_cum, y = avg.win.depth)) +
  stat_smooth( method = "loess", span = 0.1, se = FALSE, color = "black")+
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "#7570B3", 
           fill = "#7570B3", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(233048607+5329939), xmax=(233048607+7791197), ymin=-Inf, ymax=Inf,
           color = "#E7298A", 
           fill = "#E7298A", alpha = 0.5, linetype = 2)+
  ylab(paste("Read depth"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))
  

## Scatter plot of gstat
L1_gstat_plot <- ggplot(L1_aligned_dataset_gstat_filtered, aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="black", linetype = 2) +
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "#7570B3", 
           fill = "#7570B3", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(233048607+5329939), xmax=(233048607+7791197), ymin=-Inf, ymax=Inf,
           color = "#E7298A", 
           fill = "#E7298A", alpha = 0.5, linetype = 2)+
  #annotate(geom = "pointrange", x=122678833, y=10.580692, ymin = 10.580692, ymax = 10.580692,  color= "red",size=.5)+
  #annotate(geom = "pointrange", x=122573736, y=19.076968, ymin = 19.076968, ymax = 19.076968,  color= "blue",size=.5)+
  #annotate(geom = "pointrange", x=122678833, y=10.580692, ymin = 10.580692, ymax = 10.580692,  color= "blue",size=.5)+   ylab(paste("G statistic"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

## Arrange plots on one graph (Figure2)

L1_alinged_Fig2 <- ggarrange(L1_gstat_plot, L1_pi_plot,L1_depth_plot,nrow = 3, ncol = 1, labels = c("D", "E", "F" ))

png("./figures/L1_aligned_gouldFig2.png", width = 30, height = 20, units = "cm",res = 300)
print(L1_alinged_Fig2)
dev.off()

## Figure 3 from Gould et al 2017 (G stat values across chromosmes 5, 8, and 14)

L1_gstat_plot_8 <- L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="black", linetype = 2) +
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype =2)+
 # annotate(geom = "pointrange", x=122678833, y=10.580692, ymin = 10.580692, ymax = 10.580692,  color= "red",size=.5)+
#  annotate(geom = "pointrange", x=122573736, y=19.076968, ymin = 19.076968, ymax = 19.076968,  color= "blue",size=.5)+
  ylab(paste("G statistic"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

L1_dxy_plot_8 <- L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  ggplot(aes(x=bp_cum, y = Dxy_win)) +
  geom_point() +
  annotate(geom = "rect", xmin=(118986430+850429), xmax=(118986430+7604769), ymin=-Inf, ymax=Inf,
           color = "#D95F02", 
           fill = "#D95F02", alpha = 0.5, linetype = 2)+
  annotate(geom = "rect", xmin=(118986430+1032334), xmax=(118986430+1246126), ymin=-Inf, ymax=Inf,
           color = "#1B9E77", 
           fill = "#1B9E77", alpha = 0.5, linetype =2)+
   ylab(paste("Dxy on chromosome 8"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

L1_gstat_plot_5 <- L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="black", linetype = 2) +
  annotate(geom = "rect", xmin=(67527869+13650670), xmax=(67527869+17847181), ymin=-Inf, ymax=Inf,
           color = "#7570B3", 
           fill = "#7570B3", alpha = 0.5, linetype = 2)+
  ylab(paste("G statistic"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

L1_gstat_plot_14 <- L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  ggplot(aes(x=bp_cum, y = G_stat)) +
  geom_point() +
  geom_hline(yintercept =(28.63410), color="black", linetype = 2) +
  annotate(geom = "rect", xmin=(233048607+5329939), xmax=(233048607+7791197), ymin=-Inf, ymax=Inf,
           color = "#E7298A", 
           fill = "#E7298A", alpha = 0.5, linetype = 2)+
  ylab(paste("G statistic"))+
  xlab("")+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

## Arrange plots on one graph (Figure3)

L1_alinged_Fig3 <- ggarrange(L1_gstat_plot_5, L1_gstat_plot_8, L1_gstat_plot_14,nrow = 3, ncol = 1, labels = c("D", "E", "F") )

png("./figures/L1_aligned_gouldFig3.png", width = 30, height = 20, units = "cm",res = 300)
print(L1_alinged_Fig3)
dev.off()

## Distribution of G Statistic window values inside and outside the chromosome inversion (figure 4)
mean(L1_aligned_dataset_gstat_filtered$G_stat)

# Add column for within inversion or not within
L1_aligned_dataset_gstat_filtered_5_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum < 67527869+13650670) 
  L1_aligned_dataset_gstat_filtered_5_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_5_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum > 67527869+17847181)
L1_aligned_dataset_gstat_filtered_5_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_5_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_5_NOTinv2, L1_aligned_dataset_gstat_filtered_5_NOTinv1)

mean(L1_aligned_dataset_gstat_filtered_5_NotINV_fig$G_stat) #Colinear regions of chromosome 5

L1_aligned_dataset_gstat_filtered_5_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="5") %>% 
  filter(between(bp_cum, (67527869+13650670),(67527869+17847181) ))
L1_aligned_dataset_gstat_filtered_5_inv$karyo <- "inv"

mean(L1_aligned_dataset_gstat_filtered_5_inv$G_stat) # chromosome 5 inversion

L1_aligned_dataset_gstat_filtered_genome_no5 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff!="5") 
L1_aligned_dataset_gstat_filtered_genome_no5$karyo <- "notinv"


L1_aligned_dataset_gstat_filtered_5_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_genome_no5,L1_aligned_dataset_gstat_filtered_5_inv, L1_aligned_dataset_gstat_filtered_5_NotINV_fig)
noninv_5 <- L1_aligned_dataset_gstat_filtered_5_figure4 %>% filter(karyo =="notinv") 


Gstat_chr5_density <- L1_aligned_dataset_gstat_filtered_5_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
 # ylab("Density")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear region"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

################################################################################

L1_aligned_dataset_gstat_filtered_8_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+850429)) 

L1_aligned_dataset_gstat_filtered_8_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+7604769))
L1_aligned_dataset_gstat_filtered_8_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_8_NOTinv1, L1_aligned_dataset_gstat_filtered_8_NOTinv2)

mean(L1_aligned_dataset_gstat_filtered_8_NotINV_fig$G_stat) # Colinear regions of chromosome 8

L1_aligned_dataset_gstat_filtered_8_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118986430+850429),(118986430+7604769)))
L1_aligned_dataset_gstat_filtered_8_inv$karyo <- "inv"

mean(L1_aligned_dataset_gstat_filtered_8_inv$G_stat) # Chromosome 8 inversion

L1_aligned_dataset_gstat_filtered_8_inv_wosmall1 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+1032334)) 

L1_aligned_dataset_gstat_filtered_8_inv_wosmall2 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+1246126))

mean(L1_aligned_dataset_gstat_filtered_8_inv_wosmall2$G_stat)

L1_aligned_dataset_gstat_filtered_genome_no8 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff!="8") 
L1_aligned_dataset_gstat_filtered_genome_no8$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_genome_no8,L1_aligned_dataset_gstat_filtered_8_inv_wosmall1, L1_aligned_dataset_gstat_filtered_8_inv_wosmall2, L1_aligned_dataset_gstat_filtered_8_NotINV_fig)

noninv_8 <- L1_aligned_dataset_gstat_filtered_8_figure4 %>% filter(karyo =="notinv") 
mean(noninv_8$G_stat)

Gstat_chr8_density <- L1_aligned_dataset_gstat_filtered_8_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
 # ylab("Density")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear region"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

################################################################################

L1_aligned_dataset_gstat_filtered_8small_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+1032334)) 

L1_aligned_dataset_gstat_filtered_8small_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8small_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+1246126))

L1_aligned_dataset_gstat_filtered_8small_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_8small_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_8_NOTinv1, L1_aligned_dataset_gstat_filtered_8_NOTinv2)

L1_aligned_dataset_gstat_filtered_8small_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118986430+1032334), (118986430+1246126)))
L1_aligned_dataset_gstat_filtered_8small_inv$karyo <- "inv"

mean(L1_aligned_dataset_gstat_filtered_8small_inv$G_stat) #Chromosome 8 small inversion

L1_aligned_dataset_gstat_filtered_8small_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_genome_no8, L1_aligned_dataset_gstat_filtered_8small_inv, L1_aligned_dataset_gstat_filtered_8_NotINV_fig)

noninv_8small <- L1_aligned_dataset_gstat_filtered_8small_figure4 %>% filter(karyo =="notinv") 
mean(noninv_8small$G_stat)

Gstat_chr8small_density <- L1_aligned_dataset_gstat_filtered_8small_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("G-statistic")+
 # ylab("Density")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear region"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))
############################################################################################################################

L1_aligned_dataset_gstat_filtered_14_NOTinv1 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum < (233048607+5329939)) 

L1_aligned_dataset_gstat_filtered_14_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_14_NOTinv2 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum > (233048607+7791197))
L1_aligned_dataset_gstat_filtered_14_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_gstat_filtered_genome_no14 <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff!="14") 
L1_aligned_dataset_gstat_filtered_genome_no14$karyo <- "notinv"


L1_aligned_dataset_gstat_filtered_14_NotINV_fig <- rbind(L1_aligned_dataset_gstat_filtered_14_NOTinv2, L1_aligned_dataset_gstat_filtered_14_NOTinv1)

mean(L1_aligned_dataset_gstat_filtered_14_NotINV_fig$G_stat)

L1_aligned_dataset_gstat_filtered_14_inv <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="14") %>% 
  filter(between(bp_cum, (233048607+5329939), (233048607+7791197)))
L1_aligned_dataset_gstat_filtered_14_inv$karyo <- "inv"

mean(L1_aligned_dataset_gstat_filtered_14_inv$G_stat)

L1_aligned_dataset_gstat_filtered_14_figure4 <- rbind(L1_aligned_dataset_gstat_filtered_genome_no14,L1_aligned_dataset_gstat_filtered_14_inv, L1_aligned_dataset_gstat_filtered_14_NotINV_fig)

noninv_14 <- L1_aligned_dataset_gstat_filtered_14_figure4 %>% filter(karyo =="notinv") 
mean(noninv_14$G_stat)

Gstat_chr14_density <- L1_aligned_dataset_gstat_filtered_14_figure4 %>% 
  ggplot(aes(x=G_stat, fill=karyo)) +
  geom_density(alpha=.3)+
 # ylab("Density")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear region"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

figure4_density_gstat <- ggarrange(Gstat_chr5_density, Gstat_chr8_density, Gstat_chr8small_density, Gstat_chr14_density, nrow = 4, ncol = 1)


png("./figures/figure4_density_gstat.png", width = 20, height = 30, units = "cm",res = 300)
print(figure4_density_gstat)
dev.off()

#######################################################
## Checking for significance of elevated G stat in non inverted, inverted, and small inv in chr8
L1_aligned_dataset_gstat_filtered_8_inv_wosmall1 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118986430+1032334)) 

L1_aligned_dataset_gstat_filtered_8_inv_wosmall2 <-  L1_aligned_dataset_gstat_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118986430+1246126))

L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined <- rbind(L1_aligned_dataset_gstat_filtered_8_inv_wosmall2, L1_aligned_dataset_gstat_filtered_8_inv_wosmall2)

L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined$karyo <- "invLARGE"

mean(L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined$G_stat)

L1_aligned_dataset_gstat_filtered_8small_inv_anova <-  L1_aligned_dataset_gstat_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118986430+1032334), (118986430+1246126)))
L1_aligned_dataset_gstat_filtered_8small_inv_anova$karyo <- "invSMALL"

L1_aligned_dataset_gstat_filtered_8ALL_inv_anova <- rbind(L1_aligned_dataset_gstat_filtered_genome_no8, L1_aligned_dataset_gstat_filtered_8small_inv_anova,L1_aligned_dataset_gstat_filtered_8_inv_wosmall_combined)

L1_aligned_dataset_gstat_filtered_8small_anova <- rbind(L1_aligned_dataset_gstat_filtered_genome_no8,L1_aligned_dataset_gstat_filtered_8ALL_inv_anova, L1_aligned_dataset_gstat_filtered_8_NotINV_fig)

library(car)
summary(aov(G_stat ~ karyo,L1_aligned_dataset_gstat_filtered_8small_anova))
boxplot.8.gstat <- boxplot(G_stat ~ karyo,L1_aligned_dataset_gstat_filtered_8small_anova, names= c("Large \n Chromosome 8 \n inversion",
                  "\n Small \n Chromosome 8 \n Inversion", "Not  \n inverted \n regions" ),
                  col = c("#66A61E", "#E6AB02", "#7570B3"), xaxt ="n", ylab = "G statistic", xlab = "")
axis(side = 1, at = 1:3, labels =  c("Large \n Chromosome 8 \n inversion",
                                "\n Small \n Chromosome 8 \n Inversion", "Not \n inverted  \n regions"), tcl = 0)

png("./figures/boxplot.8.gstat.png", width = 15, height = 15, units = "cm",res = 300)
boxplot.8.gstat
dev.off()


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
  filter(bp_cum < 67487358+13650670) 
L1_aligned_dataset_nona_tajD_5_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_5_NOTinv2 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(bp_cum > 67487358+17847181)
L1_aligned_dataset_nona_tajD_5_NOTinv2$karyo <- "notinv"


L1_aligned_dataset_nona_tajD_5_NotINV_fig <- rbind(L1_aligned_dataset_nona_tajD_5_NOTinv2, L1_aligned_dataset_nona_tajD_5_NOTinv1)

L1_aligned_dataset_nona_tajD_5_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="5") %>% 
  filter(between(bp_cum, (67487358+13650670),(67487358+17847181) ))
L1_aligned_dataset_nona_tajD_5_inv$karyo <- "inv 5"

mean(L1_aligned_dataset_nona_tajD_5_inv$TajD_IA)
mean(L1_aligned_dataset_nona_tajD_5_inv$TajD_CP)

L1_aligned_dataset_nona_tajD_5 <- rbind(L1_aligned_dataset_nona_tajD_5_inv, L1_aligned_dataset_nona_tajD_5_NotINV_fig)

################################################################################

L1_aligned_dataset_nona_tajD_8_NOTinv1 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118900276+850429 )) 
L1_aligned_dataset_nona_tajD_8_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_8_NOTinv2 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118900276+7604769))
L1_aligned_dataset_nona_tajD_8_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_8_NotINV_fig <- rbind( L1_aligned_dataset_nona_tajD_8_NOTinv1, L1_aligned_dataset_nona_tajD_8_NOTinv2)

L1_aligned_dataset_nona_tajD_8_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118900276+850429),(118900276+7604769)))
L1_aligned_dataset_nona_tajD_8_inv$karyo <- "inv 8"

mean(L1_aligned_dataset_nona_tajD_8_inv$TajD_IA)
mean(L1_aligned_dataset_nona_tajD_8_inv$TajD_CP)

#removing small inversion from 8 inversion
L1_aligned_dataset_nona_tajD_8_inv_wosmalL1 <-  L1_aligned_dataset_nona_tajD_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum < (118900276+1032334)) 

L1_aligned_dataset_nona_tajD_8_inv_wosmall2 <-  L1_aligned_dataset_nona_tajD_8_inv %>% 
  filter(scaffold=="8") %>% 
  filter(bp_cum > (118900276+1246126))

L1_aligned_dataset_nona_tajD_8 <- rbind(L1_aligned_dataset_nona_tajD_8_inv_wosmalL1, L1_aligned_dataset_nona_tajD_8_inv_wosmall2, L1_aligned_dataset_nona_tajD_8_NotINV_fig)

L1_aligned_dataset_nona_tajD_8_wosmall_combined <- rbind(L1_aligned_dataset_nona_tajD_8_inv_wosmalL1, L1_aligned_dataset_nona_tajD_8_inv_wosmall2)

mean(L1_aligned_dataset_nona_tajD_8_wosmall_combined$TajD_IA)
mean(L1_aligned_dataset_nona_tajD_8_wosmall_combined$TajD_CP)

################################################################################

L1_aligned_dataset_nona_tajD_8small_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="8") %>% 
  filter(between(bp_cum, (118900276+1032334), (118900276+1246126)))

L1_aligned_dataset_nona_tajD_8small_inv$karyo <- "inv 8 SMALL"

mean(L1_aligned_dataset_nona_tajD_8small_inv$TajD_IA)
mean(L1_aligned_dataset_nona_tajD_8small_inv$TajD_CP)

############################################################################################################################

L1_aligned_dataset_nona_tajD_14_NOTinv1 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum < (232978286+5329939)) 

L1_aligned_dataset_nona_tajD_14_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_nona_tajD_14_NOTinv2 <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(bp_cum > (232978286+7791197))
L1_aligned_dataset_nona_tajD_14_NOTinv2$karyo <- "notinv"


L1_aligned_dataset_nona_tajD_14_NotINV_fig <- rbind(L1_aligned_dataset_nona_tajD_14_NOTinv2, L1_aligned_dataset_nona_tajD_14_NOTinv1)

L1_aligned_dataset_nona_tajD_14_inv <-  L1_aligned_dataset_nona_tajD %>% 
  filter(scaffold=="14") %>% 
  filter(between(bp_cum, (232978286+5329939), (232978286+7791197)))
L1_aligned_dataset_nona_tajD_14_inv$karyo <- "inv 14"

mean(L1_aligned_dataset_nona_tajD_14_inv$TajD_IA)
mean(L1_aligned_dataset_nona_tajD_14_inv$TajD_CP)

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


png("./figures/tajD_pot_L1.png", width =15, height = 15, units = "cm",res = 300)
print(tajD_pot_L1)
dev.off()
# 
# loath %>% 
#   group_by(Karyo, tajtype) %>% 
#   summarise(m=mean(TajD))

write.csv(L1_pi_top1_gstat, file = "./data_files/L1_pi_top1_gstat.csv")
write.csv(L1_pi_top1, file="./data_files/L1_pi_top1.csv")
write.csv(L1_pi_bottom1, file="./data_files/L1_pi_bottom1.csv")

## Distribution of Pi Statistic window values inside and outside the chromosome inversion
# Add column for within inversion or not within

mean(L1_aligned_dataset_pi_filtered$piIA)
mean(L1_aligned_dataset_pi_filtered$piCP)

L1_aligned_dataset_pi_filtered_5_NOTinv1 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum < 67485786+13650670) 
L1_aligned_dataset_pi_filtered_5_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_5_NOTinv2 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="5") %>% 
  filter(bp_cum > 67485786+17847181)
L1_aligned_dataset_pi_filtered_5_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_5_NotINV_fig <- rbind(L1_aligned_dataset_pi_filtered_5_NOTinv2, L1_aligned_dataset_pi_filtered_5_NOTinv1)

mean(L1_aligned_dataset_pi_filtered_5_NotINV_fig$piRatio)
mean(L1_aligned_dataset_pi_filtered_5_NotINV_fig$piCP)
mean(L1_aligned_dataset_pi_filtered_5_NotINV_fig$piIA)

L1_aligned_dataset_pi_filtered_5_inv <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="5") %>% 
  filter(between(bp_cum, (67485786+13650670),(67485786+17847181) ))
L1_aligned_dataset_pi_filtered_5_inv$karyo <- "inv"

mean(L1_aligned_dataset_pi_filtered_5_inv$piRatio)
mean(L1_aligned_dataset_pi_filtered_5_inv$piCP)
mean(L1_aligned_dataset_pi_filtered_5_inv$piIA)

L1_aligned_dataset_pi_filtered_no5 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff!="5")
L1_aligned_dataset_pi_filtered_no5$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_5_figure4 <- rbind(L1_aligned_dataset_pi_filtered_no5,L1_aligned_dataset_pi_filtered_5_inv, L1_aligned_dataset_pi_filtered_5_NotINV_fig)

noninv_5pi <- L1_aligned_dataset_pi_filtered_5_figure4 %>% filter(karyo =="notinv") 
mean(noninv_5pi$piCP)
mean(noninv_5pi$piIA)
mean(noninv_5pi$piRatio)

pi_chr5_density <- L1_aligned_dataset_pi_filtered_5_figure4 %>% 
  ggplot(aes(x=piRatio, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("\u03c0 Ratio")+
  ggtitle("inv_chr5A")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(plot.title = element_text(color = "black", size=14, face="bold", hjust = 0.5),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

################################################################################

L1_aligned_dataset_pi_filtered_8_NOTinv1 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118943830+850429)) 

L1_aligned_dataset_pi_filtered_8_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_8_NOTinv2 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118943830+7604769))
L1_aligned_dataset_pi_filtered_8_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_8_NotINV_fig <- rbind(L1_aligned_dataset_pi_filtered_8_NOTinv1, L1_aligned_dataset_pi_filtered_8_NOTinv2)

mean(L1_aligned_dataset_pi_filtered_8_NotINV_fig$piRatio)
mean(L1_aligned_dataset_pi_filtered_8_NotINV_fig$piCP)
mean(L1_aligned_dataset_pi_filtered_8_NotINV_fig$piIA)

mean(L1_aligned_dataset_pi_filtered_8_NotINV_fig$Dxy_win)


L1_aligned_dataset_pi_filtered_8_inv <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118943830+850429),(118943830+7604769)))
L1_aligned_dataset_pi_filtered_8_inv$karyo <- "inv"

mean(L1_aligned_dataset_pi_filtered_8_inv$piRatio)
mean(L1_aligned_dataset_pi_filtered_8_inv$piCP)
mean(L1_aligned_dataset_pi_filtered_8_inv$piIA)

mean(L1_aligned_dataset_pi_filtered_8_inv$Dxy_win)

L1_aligned_dataset_pi_filtered_8_inv_wosmall1 <-  L1_aligned_dataset_pi_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118943830+1032334)) 

L1_aligned_dataset_pi_filtered_8_inv_wosmall2 <-  L1_aligned_dataset_pi_filtered_8_inv %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118943830+1246126))

L1_aligned_dataset_pi_filtered_8_combined <- rbind(L1_aligned_dataset_pi_filtered_8_inv_wosmall1, L1_aligned_dataset_pi_filtered_8_inv_wosmall2)

mean(L1_aligned_dataset_pi_filtered_8_combined$piRatio)
mean(L1_aligned_dataset_pi_filtered_8_combined$piIA)
mean(L1_aligned_dataset_pi_filtered_8_combined$piCP)

mean(L1_aligned_dataset_pi_filtered_8_combined$Dxy_win)

L1_aligned_dataset_pi_filtered_no8 <-  L1_aligned_dataset_pi_filtered_8_inv %>% 
  filter(scaff!="8")
L1_aligned_dataset_pi_filtered_no8$karyo <- "notinv"


L1_aligned_dataset_pi_filtered_8_figure4 <- rbind(L1_aligned_dataset_pi_filtered_no8,L1_aligned_dataset_pi_filtered_8_inv_wosmall1, L1_aligned_dataset_pi_filtered_8_inv_wosmall2, L1_aligned_dataset_pi_filtered_8_NotINV_fig)


L1_aligned_dataset_pi_filtered_8_wosmall_combined <- rbind(L1_aligned_dataset_pi_filtered_8_inv_wosmall1, L1_aligned_dataset_pi_filtered_8_inv_wosmall2)
mean(L1_aligned_dataset_pi_filtered_8_wosmall_combined$piRatio)

noninv_8pi <- L1_aligned_dataset_pi_filtered_8_figure4 %>% filter(karyo =="notinv") 
mean(noninv_8pi$piCP)
mean(noninv_8pi$piIA)
mean(noninv_8pi$piRatio)

pi_chr8_density <- L1_aligned_dataset_pi_filtered_8_figure4 %>% 
  ggplot(aes(x=piRatio, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("\u03c0 Ratio")+
  ggtitle("inv_chr8A")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(plot.title = element_text(color = "black", size=14, face="bold", hjust = 0.5),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

################################################################################

L1_aligned_dataset_pi_filtered_8small_NOTinv1 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum < (118943830+1032334)) 

L1_aligned_dataset_pi_filtered_8small_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_8small_NOTinv2 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  filter(bp_cum > (118943830+1246126))
L1_aligned_dataset_pi_filtered_8small_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_8small_NotINV_fig <- rbind(L1_aligned_dataset_pi_filtered_8_NOTinv1, L1_aligned_dataset_pi_filtered_8_NOTinv2)

L1_aligned_dataset_pi_filtered_8small_inv <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="8") %>% 
  filter(between(bp_cum, (118943830+1032334), (118943830+1246126)))
L1_aligned_dataset_pi_filtered_8small_inv$karyo <- "inv"

mean(L1_aligned_dataset_pi_filtered_8small_inv$piRatio)
mean(L1_aligned_dataset_pi_filtered_8small_inv$piCP)
mean(L1_aligned_dataset_pi_filtered_8small_inv$piIA)

mean(L1_aligned_dataset_pi_filtered_8small_inv$Dxy_win)

L1_aligned_dataset_pi_filtered_8small_figure4 <- rbind(L1_aligned_dataset_pi_filtered_no8,L1_aligned_dataset_pi_filtered_8small_inv, L1_aligned_dataset_pi_filtered_8small_NotINV_fig)

pi_chr8small_density <- L1_aligned_dataset_pi_filtered_8small_figure4 %>% 
  ggplot(aes(x=piRatio, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("\u03c0 Ratio")+
  ggtitle("inv_chr8B")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(plot.title = element_text(color = "black", size=14, face="bold", hjust = 0.5),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))
############################################################################################################################

L1_aligned_dataset_pi_filtered_14_NOTinv1 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum < (233048607+5329939)) 

L1_aligned_dataset_pi_filtered_14_NOTinv1$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_14_NOTinv2 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="14") %>% 
  filter(bp_cum > (233048607+7791197))
L1_aligned_dataset_pi_filtered_14_NOTinv2$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_14_NotINV_fig <- rbind(L1_aligned_dataset_pi_filtered_14_NOTinv2, L1_aligned_dataset_pi_filtered_14_NOTinv1)

mean(L1_aligned_dataset_pi_filtered_14_NotINV_fig$piRatio)
mean(L1_aligned_dataset_pi_filtered_14_NotINV_fig$piCP)
mean(L1_aligned_dataset_pi_filtered_14_NotINV_fig$piIA)

L1_aligned_dataset_pi_filtered_14_inv <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff=="14") %>% 
  filter(between(bp_cum, (233048607+5329939), (233048607+7791197)))
L1_aligned_dataset_pi_filtered_14_inv$karyo <- "inv"

mean(L1_aligned_dataset_pi_filtered_14_inv$piRatio)
mean(L1_aligned_dataset_pi_filtered_14_inv$piCP)
mean(L1_aligned_dataset_pi_filtered_14_inv$piIA)
 
L1_aligned_dataset_pi_filtered_no14 <-  L1_aligned_dataset_pi_filtered %>% 
  filter(scaff!="14")
L1_aligned_dataset_pi_filtered_no14$karyo <- "notinv"

L1_aligned_dataset_pi_filtered_14_figure4 <- rbind(L1_aligned_dataset_pi_filtered_no14,L1_aligned_dataset_pi_filtered_14_inv, L1_aligned_dataset_pi_filtered_14_NotINV_fig)

noninv_14pi <- L1_aligned_dataset_pi_filtered_14_figure4 %>% filter(karyo =="notinv") 
mean(noninv_14pi$piCP)
mean(noninv_14pi$piIA)
mean(noninv_14pi$piRatio)

pi_chr14_density <- L1_aligned_dataset_pi_filtered_14_figure4 %>% 
  ggplot(aes(x=piRatio, fill=karyo)) +
  geom_density(alpha=.3)+
  xlab("\u03c0 Ratio")+
  ggtitle("inv_chr14A")+
  scale_fill_manual(name = "Orientation",labels = c("inversion", "colinear"),values = c("#D95F02","#1B9E77"))+
  theme_bw()+ theme(plot.title = element_text(color = "black", size=14, face="bold", hjust = 0.5),
                    panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), 
                    axis.line = element_line(colour = "black"))

figure4_density_pi <- ggarrange(pi_chr5_density, pi_chr8_density, pi_chr8small_density, pi_chr14_density, nrow = 4, ncol = 1) 


png("./figures/figure4_density_pi.png", width = 20, height = 30, units = "cm",res = 300)
print(figure4_density_pi)
dev.off()

## Inversion specific plots

chr5_inv_g_pi_plot <- ggarrange(pi_chr5_density,Gstat_chr5_density, nrow = 2, legend = "none")


chr8_inv_g_pi_plot <- ggarrange(pi_chr8_density,Gstat_chr8_density, nrow = 2, common.legend = F, legend = "none")


chr8small_inv_g_pi_plot <- ggarrange(pi_chr8small_density,Gstat_chr8small_density, nrow = 2, common.legend = F, legend = "none")


chr14_inv_g_pi_plot <- ggarrange(pi_chr14_density,Gstat_chr14_density, nrow = 2, common.legend = F, legend = "none")


inversion_stat_figure <- ggarrange(chr5_inv_g_pi_plot, chr8_inv_g_pi_plot, chr8small_inv_g_pi_plot, 
  chr14_inv_g_pi_plot, #labels = c("inv_chr5A", "inv_chr8A","inv_chr8B", "inv_chr14"), 
  nrow = 1, common.legend = T) %>% 
gridExtra::grid.arrange(get_legend(pi_chr5_density), heights = unit(c(80, 5), "mm"))

png("./figures/inversion_stat_figure.png", width = 28, height = 15, units = "cm",res = 300)
print(inversion_stat_figure)
dev.off()
