#####################
##### Hi-C Maps #####
#####################

# Set a working directory
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

# Load libraries
library(HiContacts)
library(dplyr)
library(ggplot2)
library(purrr)
library(GenomicRanges)
library(patchwork)

#####################
###### Samples ######
#####################

# Prepare metadata file
metadata <- as.data.frame(matrix(nrow = 8, ncol = 3))
colnames(metadata) <- c("strain", "h.p.i", "replicate")
rownames(metadata) <- list.files("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
metadata$strain <- c("KO", "KO", "WT", "WT", "KO", "KO", "WT", "WT")
metadata$h.p.i <- c(40, 40, 40, 40, 16, 16, 16, 16)
metadata$replicate <- c(2, 1, 2, 1, 2, 1, 2, 1)
metadata$ID <- paste0(metadata$strain, "_", metadata$h.p.i, "h_rep", metadata$replicate)

samples <- rownames(metadata)

########################
##### Compartments #####
########################

# setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
# 
# # Import Hi-C contact matrices
# hics <- list()
# for(i in 1:nrow(metadata)){
#   hic <- import(paste0(rownames(metadata)[i],"/inter_30.mcool"), resolution = 100000)
#   hics <- append(hics, hic)
# }
# names(hics) <- metadata$ID
# 
# for(i in 1:nrow(metadata)){
#   sample <- rownames(metadata)[i]
#   setwd(paste0("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results/",sample))
#   
#   # Chromosome list
#   chr <- c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3")
#   
#   compartmetns <- data.frame()
#   for(j in 1:length(chr)){
#     compar <- read.table(paste0(chr[j],"_eigen.txt"), sep = "\t")
#     compar$start <- as.numeric(seq(1, nrow(compar)*10000 - 1, 10000))
#     compar$end <- as.numeric(seq(10000, nrow(compar)*10000 + 1, 10000))
#     compar$chr <- chr[j]
#     colnames(compar)[1] <- "E1"
#     compar$compartment <- "A"
#     compar[compar$E1 < 0,]$compartment <- "B"
#     compartmetns <- rbind(compartmetns, compar)
#   }
#   
#   compartmetns <- makeGRangesFromDataFrame(compartmetns, keep.extra.columns = T) 
#   
#   setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
#   topologicalFeatures(hics[[metadata$ID[i]]], "compartments") <- compartmetns
#   metadata(hics[[metadata$ID[i]]])$eigens <- compartmetns 
#   
#   # Plot
#   chr.list <- toupper(c("Pf3D7_01_v3","Pf3D7_02_v3", 
#                         "Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3"))
#   
#   for(k in 1:length(chr.list)){
#     setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
#     corr.map <- autocorrelate(hics[[metadata$ID[i]]][chr.list[k]])
#     ##  
#     
#     p1 <- plotMatrix(corr.map, use.scores = 'autocorrelated', scale = 'linear', limits = c(-1, 1), caption = FALSE)
#     eigen <- coverage(metadata(hics[[metadata$ID[i]]])$eigens, weight = 'E1')[[k]]
#     eigen_df <- tibble(pos = cumsum(runLength(eigen)), eigen = runValue(eigen))
#     p2 <- ggplot(eigen_df, aes(x = pos, y = eigen)) + 
#       geom_area() + 
#       theme_void() + 
#       coord_cartesian(expand = FALSE) + 
#       labs(x = "Genomic position", y = "Eigenvector value")
#     p <- wrap_plots(p1, p2, ncol = 1, heights = c(10, 1))
#     
#     setwd(paste0("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/compartments/",metadata$ID[i]))
#     png(paste0(chr.list[k], "_compartments.png"), width = 200, height = 200, units = "mm", res = 200)
#     p
#     dev.off()
#   }
# }