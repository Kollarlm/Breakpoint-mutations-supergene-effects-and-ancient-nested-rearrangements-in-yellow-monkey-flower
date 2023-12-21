#Set some variables
tandem_cutoff=5 #Must be within this number of genes of a gene with same arrayID to be classified as "tandem"

#Load libraries
library(GENESPACE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(venneuler)

runwd <- file.path("//Users/lesliekollar/Desktop/mimulus_genome_R/mimulus_genome/Mimulus_guttatus_nanopore_genome_v1.2")
setwd(runwd)

#Read in the pangenome database
pgdb <- fread("results/S1_pangenomeDB.txt.gz",na.strings = c("", "NA"))
pgff <- fread("results/gffWithOgs.txt.gz",na.strings = c("", "NA"))

#Create two new tables for each of the genomes, removing genes not in the pangenome
L1 <- subset(pgdb, genome=="L1" & !is.na(pgChr))
S1 <- subset(pgdb, genome=="S1" & !is.na(pgChr))

#Create a separate list of unique genes to each genome, not in pangenome
L1u <- subset(pgdb, genome=="L1" & is.na(pgChr))
S1u <- subset(pgdb, genome=="S1" & is.na(pgChr))

#Combine the files based on the pangenome ID
LS <- full_join(L1,S1,by=("pgID"))

#Calculate putative CNV & PAV from pgID
cnv <- dcast(subset(pgdb, !is.na(pgChr)), pgChr + pgOrd + pgID ~ genome, value.var="id", fun.aggregate=length)

#If isDirectSyn & isArrayRep is TRUE for both paired genes in both genomes, it is a syntenic pair
#Subset syntenic genes and write out lists of these genes
synt <- subset(LS, isArrayRep.x & isDirectSyn.x & isArrayRep.y & isDirectSyn.y)
write.table(data.frame(L1=synt$id.x,S1=synt$id.y),file="results/syntelogs.tsv",
            sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(unique(synt$id.x),file="results/L1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(unique(synt$id.y),file="results/S1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#Get High-confidence 1x1 syntelogs and write out lists of these genes
synt1x1 <- subset(synt, pgID %in% subset(cnv, L1==1 & S1==1)$pgID)
write.table(data.frame(L1=synt1x1$id.x,S1=synt1x1$id.y),file="results/syntelogs1x1.tsv",
            sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(synt1x1$id.x,file="results/L1_1x1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(synt1x1$id.y,file="results/S1_1x1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#Get sets of non-syntenic genes and write out lists of these genes
nonsynt <- setdiff(LS, synt)
L1nonsynt <- subset(nonsynt, !(id.x %in% synt$id.x))
S1nonsynt <- subset(nonsynt, !(id.y %in% synt$id.y))
write.table(c(unique(L1nonsynt$id.x),L1u$id),file="results/L1_nonsyntenic.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(c(unique(S1nonsynt$id.y),S1u$id),file="results/S1_nonsyntenic.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

#Get genes in arraypairs arrays
#Add the arrayID to pgdb
pge <- setorder(merge(pgdb,pgff[,c(2,11)],by="id"), pgChr, pgOrd, na.last = T)
#Count count arrays with more than one gene in them
arraycount <- as.data.frame(table(unique(pge[,c(1,15)])$arrayID))
#subset out those genes with an array count of more than 1
arraypairs <- setorder(subset(pge, arrayID %in% subset(arraycount, Freq > 1)$Var1), arrayID)
#Create output dataframe
columns=c("genome","arrayID","id","isTandem","geneCount","arrayGenes")
tandem <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(tandem) <- columns
#Now for a complicated loop to identify tandem arrays
for(i in unique(arraypairs$arrayID)){
  #some arraypairs arrays may belong to multiple orthogroups, collapse these so each gene is represented once
  x <- setorder(unique(subset(arraypairs, arrayID==i)[,c(1,8,9,10,15)]), ord)
  genome <- unique(x$genome)
  #Check to make sure all of these are on the same chromosome!
  if(uniqueN(x$chr)==1){
    #Most array pairs are simple arraypairs duplicates and faster to process
    #If the number of genes is 2 and the distance is <= our cutoff classify as arraypairs duplicates
    if(nrow(x)==2 & abs(x[1]$ord - x[2]$ord) <= tandem_cutoff){
      for(gene in x$id){
        tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=TRUE,geneCount=2,
                                          arrayGenes=paste(x$id,collapse=", ")))
      } 
      #If the number of genes is 2 and distance > cutoff, classify as a dispersed
    } else if(nrow(x)==2 & abs(x[1]$ord - x[2]$ord) > tandem_cutoff) {
      tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=FALSE,geneCount=NA,
                                        arrayGenes=NA))
      #If the arraypairs have more than one gene, it requires more work
    } else {
      #Create an empty dataframe to calculate pairwise distances between all genes
      x2 <- data.frame(matrix(ncol=nrow(x),nrow=nrow(x)))
      #Set row & column names for the table
      row.names(x2) <- colnames(x2) <- x$id
      #Calculate absolute distances between all genes
      for(gene in 1:nrow(x)){
        x2[,gene] <- abs(x$ord - x[gene]$ord)
      }
      for(gene in x$id){
        #tandem if a gene's closest array pair is greater than the cutoff
        #If so, classify as FALSE
        if(min(x2[x2[,gene] != 0,gene]) > tandem_cutoff){
          tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=FALSE,geneCount=NA,
                                            arrayGenes=NA))
          #Otherwise, the gene is part of an arraypairs array
        } else {
          #Get a list of all those genes within the cutoff distance of the arraypairs array
          x3 <- row.names(x2[x2[gene,] <= tandem_cutoff,])
          #Now we are going to loop back over those genes in that initial list and search
          #for other genes which may be within range of the arraypairs array, but not the initial gene
          #these will get added to the list
          newgenes <- 1 #set this to initialize the for loop
          #Keep the loop going until no new genes are added
          while(newgenes != 0){
            #Create x4 list same as x3
            x4 <- x3
            #Loop over genes in x3
            for(gene2 in x3){
              #Extend out the array and add the new genes to x4
              x4 <- unique(c(x4,row.names(x2[x2[gene2,] <= tandem_cutoff,])))
            }
            #Set newgenes to the difference between x4 & x3
            #The loop will terminate when no new genes are added
            newgenes <- length(x4) - length(x3)
            #Set x3 to now be equal to x4
            x3 <- x4
          }
          #Now we wll report it
          tandem <- rbind(tandem,data.frame(genome=genome,arrayID=i,id=gene,isTandem=TRUE,
                                            geneCount=length(x3),arrayGenes=paste(x3,collapse=", ")))
        }
      }
    }	
    #If the genes are not on the same chromosome, they can't be a arraypairs array. Throw a warning.
  } else {
    print(paste("Warning! Array ID",unique(x$arrayID),"genes on different chromosomes. Skipping",sep=" "))
  }
}

#Cleanup
rm(x,x2,x3,x4,gene,gene2,i,columns)
#Output tables of the tandem arrays
write.table(subset(tandem,genome=="L1" & isTandem),file="results/L1_tandem_arrays.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(subset(tandem,genome=="S1" & isTandem),file="results/S1_tandem_arrays.tsv",sep="\t",quote=FALSE,row.names=FALSE)

# Adding gene id to CNV
dictionary <- pgdb[,1:3] 
dictionary <- cbind(dictionary,pgdb$id) 
dictionary$id <- dictionary$V2 # renaming column


dictionary2 <- dictionary %>% 
  group_by(pgChr, pgOrd, pgID ) %>% #Gather by the columns that are repetative
  summarise_at(vars(-group_cols()), ~toString(.[!is.na(.)], collapse="_")) # Collapse the id column into a single cell seperated by ","

dictionary2$V2 <- NULL #removing column

cnv_withID <- cnv %>% left_join(dictionary2,by=c("pgChr", "pgOrd" ,"pgID"))

### Sub-setting the  dataframe 

# Genes absent in one ecotype but not the other... specific to L1 or S1
# This method keeps duplicated genes in a single column
# absent_s1_present_L1 <- cnv_withID[sapply(lapply(cnv$S1, `%in%`, "0"), any),]
# write.table(absent_s1_present_L1, file="supplemental_tables/absent_S1_present_L1.txt",
#             quote=FALSE,col.names=FALSE,row.names=FALSE)
# 
# absent_L1_present_S1 <- cnv_withID[sapply(lapply(cnv$L1, `%in%`, "0"), any),]
# write.table(absent_L1_present_S1, file="supplemental_tables/absent_L1_present_S1.txt",
#             quote=FALSE,col.names=FALSE,row.names=FALSE)

# Copy number of genes in S1, and L1
cnv_withID$cn_diff <- cnv_withID$L1 - cnv_withID$S1 
# write.table(cnv_withID, file="supplemental_tables/L1_S1_copy_numbers.txt",
#             quote=FALSE,col.names=FALSE,row.names=FALSE)

# Same number of copies and not a CNV
same_cn <- cnv_withID %>% 
  subset(cn_diff == 0)

# Different number of copies and thus CNV.
# This includes PAVs as well
diff_cn <- cnv_withID %>%
  subset(cn_diff != 0)
write.table(diff_cn, file = "supplemental_tables/L1_S1_different_copy_numbers.txt",
            quote=FALSE, col.names = FALSE, row.names = FALSE)

# Tandem arrays
tandem_arrays <- as.data.frame(table((tandem$arrayID)))
#colnames(tandem_arrays) <- c("Number of genes", "Frequency")
write.table(tandem_arrays, file = "supplemental_tables/tandem_arrays_S1L1.txt",
            quote=FALSE, col.names = FALSE, row.names = FALSE)

L1_tandem <- tandem %>% 
  filter(genome == "L1") %>% # How many tandem duplicates are in L1?
  select(arrayID) %>% 
  unique()
  
S1_tandem <- tandem %>% 
  filter(genome == "S1") %>% # How many tandem duplicates are in L1?
  select(arrayID) %>% 
  unique()

#Plotting tandem arrays
ggplot(data=tandem_arrays, aes( x=tandem_arrays$Freq)) +
  geom_histogram()+
  #geom_text(aes(label=tandem_arrays$Freq), vjust=1.6, color="black", size=3) +
  labs(y= "Frequency", x= "Size of Tandem Array")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    

## Gene Content Summary
chr <- c("chr1", "chr2"," chr3", "chr4", "chr5", "chr6",  "chr7", "chr8", "chr9",  "chr10", 
           "chr11", "chr12", "chr13", "chr14")

cnv_withID_4inversions <- cnv %>% left_join(dictionary,by=c("pgChr", "pgOrd" ,"pgID"), multiple="all") #File for filtering out inversion genes. Not collapsed

cnv_chr <- cnv_withID_4inversions%>% filter(
  pgChr %in% chr)

cnv_chr$id <- gsub("\\.[1]*$", "", cnv_chr$id)

# genes_shared <- cnv_chr %>%   filter(L1 != 0) %>% # Use syntenic genes
#           filter(S1 != 0)

genes_L1_specific <- cnv_chr %>% 
  filter(S1 == 0)  


genes_S1_specific <- cnv_chr %>% 
  filter(L1 == 0)  

#Filtering by chr8
chr_8_genes1 <- read.table("L1_genes_chr8_inversion.gff")
chr_8_genes2 <- read.table("S1_genes_chr8_inversion.gff")
chr_8_genes <- rbind(chr_8_genes1, chr_8_genes2)
chr_8_genes$V9 <- gsub(";Name=MgL1_08g[^\\d]+", "", chr_8_genes$V9)
chr_8_genes$V9 <- gsub(";Name=MgS1_08g[^\\d]+", "", chr_8_genes$V9)
colnames(chr_8_genes)[colnames(chr_8_genes)== "V9"] <- 'id'
chr_8_genes$id <- gsub("ID=", "", chr_8_genes$id)

#Only in the chr 8 inversion
cnv_chr8 <- cnv_chr %>%  
  filter(pgChr == "chr8") %>% # It doesnt matter if you use pgChr.x or .y because they are the same
  filter(id %in% chr_8_genes$id) %>% 
  group_by(pgChr, pgOrd, pgID,L1,S1 ) %>% #Gather by the columns that are repetitive
  summarise_at(vars(-group_cols()), ~toString(.[!is.na(.)], collapse="_")) # Collapse the id column into a single cell seperated by ","

# Both L1 and S1 have to have copies
cnv_chr8_shared <- cnv_chr8 %>% 
  filter(L1 != 0) %>% 
  filter(S1 != 0)

# L1 doesnt have
cnv_chr8_notshared_L1 <- cnv_chr8 %>% 
  filter(L1 == 0) 


# S1 doesnt have 
cnv_chr8_notshared_S1 <- cnv_chr8 %>% 
  filter(S1 == 0) 

# CNVs
cnv_chr8$cn_diff <- cnv_chr8$L1 - cnv_chr8$S1 
cnv_chr8_DIFF <- cnv_chr8 %>% 
  filter(cn_diff != 0)

# PAVs
PAV_cnv_chr8_A <- cnv_chr8 %>% 
  filter(S1 == 0) %>% 
  filter(L1 != 0)

PAV_cnv_chr8_B <- cnv_chr8 %>% 
  filter(S1 != 0) %>% 
  filter(L1 == 0)

PAV_cnv_chr8 <- rbind(PAV_cnv_chr8_B, PAV_cnv_chr8_A)

#Filtering by chr8small
chr_8small_genes1 <- read.table("L1_genes_chr8_inversion_small.gff")
chr_8small_genes2 <- read.table("S1_genes_chr8_inversion_small.gff")
chr_8small_genes <- rbind(chr_8small_genes1, chr_8small_genes2)
chr_8small_genes$V9 <- gsub(";Name=MgL1_08g[^\\d]+", "", chr_8small_genes$V9)
chr_8small_genes$V9 <- gsub(";Name=MgS1_08g[^\\d]+", "", chr_8small_genes$V9)
colnames(chr_8small_genes)[colnames(chr_8small_genes)== "V9"] <- 'id'
chr_8small_genes$id <- gsub("ID=", "", chr_8small_genes$id)

cnv_chr8small <- cnv_chr %>%  
  filter(pgChr == "chr8") %>% # It doesnt matter if you use pgChr.x or .y because they are the same
  filter(id %in% chr_8small_genes$id) %>% 
  group_by(pgChr, pgOrd, pgID,L1,S1 ) %>% #Gather by the columns that are repetative
  summarise_at(vars(-group_cols()), ~toString(.[!is.na(.)], collapse="_")) # Collapse the id column into a single cell seperated by ","

# Both L1 and S1 have to have copies
cnv_chr8small_shared <- cnv_chr8small %>% 
  filter(L1 != 0) %>% 
  filter(S1 != 0)

# L1 doesnt have
cnv_chr8small_notshared_L1 <- cnv_chr8small %>% 
  filter(L1 == 0) 


# S1 doesnt have 
cnv_chr8small_notshared_S1 <- cnv_chr8small %>% 
  filter(S1 == 0)  

# CNVs
cnv_chr8small$cn_diff <- cnv_chr8small$L1 - cnv_chr8small$S1 
cnv_chr8small_DIFF <- cnv_chr8small %>% 
  filter(cn_diff != 0)

# PAVs
PAV_cnv_chr8small_A <- cnv_chr8small %>% 
  filter(S1 == 0) %>% 
  filter(L1 != 0)

PAV_cnv_chr8small_B <- cnv_chr8small %>% 
  filter(S1 != 0) %>% 
  filter(L1 == 0)

PAV_cnv_chr8small <- rbind(PAV_cnv_chr8small_B, PAV_cnv_chr8small_A)

#Filtering by chr5
chr_5_genes1 <- read.table("L1_genes_chr5_inversion.gff")
chr_5_genes2 <- read.table("S1_genes_chr5_inversion.gff")
chr_5_genes <- rbind(chr_5_genes1, chr_5_genes2)
chr_5_genes$V9 <- gsub(";Name=MgL1_05g[^\\d]+", "", chr_5_genes$V9)
chr_5_genes$V9 <- gsub(";Name=MgS1_05g[^\\d]+", "", chr_5_genes$V9)
colnames(chr_5_genes)[colnames(chr_5_genes)== "V9"] <- 'id'
chr_5_genes$id <- gsub("ID=", "", chr_5_genes$id)

cnv_chr5 <- cnv_chr %>%  
  filter(pgChr == "chr5") %>% # It doesnt matter if you use pgChr.x or .y because they are the same
  filter(id %in% chr_5_genes$id) %>% 
  group_by(pgChr, pgOrd, pgID,L1,S1 ) %>% #Gather by the columns that are repetative
  summarise_at(vars(-group_cols()), ~toString(.[!is.na(.)], collapse="_")) # Collapse the id column into a single cell seperated by ","

# Both L1 and S1 have to have copies
cnv_chr5_shared <- cnv_chr5 %>% 
  filter(L1 != 0) %>% 
  filter(S1 != 0)

# L1 doesnt have
cnv_chr5_notshared_L1 <- cnv_chr5 %>% 
  filter(L1 == 0) 


# S1 doesnt have 
cnv_chr5_notshared_S1 <- cnv_chr5 %>% 
  filter(S1 == 0)  

# CNVs
cnv_chr5$cn_diff <- cnv_chr5$L1 - cnv_chr5$S1 
cnv_chr5_DIFF <- cnv_chr5 %>% 
  filter(cn_diff != 0)

# PAVs
PAV_cnv_chr5_A <- cnv_chr5 %>% 
  filter(S1 == 0) %>% 
  filter(L1 != 0)

PAV_cnv_chr5_B <- cnv_chr5 %>% 
  filter(S1 != 0) %>% 
  filter(L1 == 0)

PAV_cnv_chr5 <- rbind(PAV_cnv_chr5_B, PAV_cnv_chr5_A)

#Filtering by chr14
chr_14_genes1 <- read.table("L1_genes_chr14_inversion.gff")
chr_14_genes2 <- read.table("S1_genes_chr14_inversion.gff")
chr_14_genes <- rbind(chr_14_genes1, chr_14_genes2)
chr_14_genes$V9 <- gsub(";Name=MgL1_14g[^\\d]+", "", chr_14_genes$V9)
chr_14_genes$V9 <- gsub(";Name=MgS1_14g[^\\d]+", "", chr_14_genes$V9)
colnames(chr_14_genes)[colnames(chr_14_genes)== "V9"] <- 'id'
chr_14_genes$id <- gsub("ID=", "", chr_14_genes$id)

cnv_chr14 <- cnv_chr %>%  
  filter(pgChr == "chr14") %>% # It doesnt matter if you use pgChr.x or .y because they are the same
  filter(id %in% chr_14_genes$id) %>% 
  group_by(pgChr, pgOrd, pgID,L1,S1 ) %>% #Gather by the columns that are repetative
  summarise_at(vars(-group_cols()), ~toString(.[!is.na(.)], collapse="_")) # Collapse the id column into a single cell seperated by ","

# Both L1 and S1 have to have copies
cnv_chr14_shared <- cnv_chr14 %>% 
  filter(L1 != 0) %>% 
  filter(S1 != 0)

# L1 doesnt have
cnv_chr14_notshared_L1 <- cnv_chr14 %>% 
  filter(L1 == 0) 


# S1 doesnt have 
cnv_chr14_notshared_S1 <- cnv_chr14 %>% 
  filter(S1 == 0)  

# CNVs
cnv_chr14$cn_diff <- cnv_chr14$L1 - cnv_chr14$S1 
cnv_chr14_DIFF <- cnv_chr14 %>% 
  filter(cn_diff != 0)

# PAVs
PAV_cnv_chr14_A <- cnv_chr14 %>% 
  filter(S1 == 0) %>% 
  filter(L1 != 0)

PAV_cnv_chr14_B <- cnv_chr14 %>% 
  filter(S1 != 0) %>% 
  filter(L1 == 0)

PAV_cnv_chr14 <- rbind(PAV_cnv_chr14_B, PAV_cnv_chr14_A)


