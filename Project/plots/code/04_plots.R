#####################
##### Hi-C Maps #####
#####################

# Set a working directory
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

# Load libraries
library(HiContacts)
library(dplyr)
library(ggplot2)

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

# Chromosome list
chr.list <- toupper(c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3"))

# Import Hi-C contact matrices
hics <- list()
for(i in 1:nrow(metadata)){
  hic <- import(paste0(rownames(metadata)[i],"/inter_30.mcool"), resolution = 10000)
  hics <- append(hics, hic)
}
names(hics) <- metadata$ID

# Plot 40 h.p.i.
for(i in 1:length(chr.list)){
  chr <- chr.list[i]
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
  merged_replicates <- list(
    KO_40h = merge(hics[["KO_40h_rep2"]][chr], hics[["KO_40h_rep1"]][chr]),
    WT_40h = merge(hics[["WT_40h_rep2"]][chr], hics[["WT_40h_rep1"]][chr]))
  maps <- imap(merged_replicates, ~ plotMatrix(.x, use.scores = 'balanced', caption = FALSE, limits = c(0, 4)) + ggtitle(.y))
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/40")
  p <- cowplot::plot_grid(plotlist = maps, nrow = 1)
  ggsave(plot = p, filename = paste0(chr, "_KO_40_WT_40.png"))
}
  

# Plot 16 h.p.i.
for(i in 1:length(chr.list)){
  chr <- chr.list[i]
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
  merged_replicates <- list(
    KO_16h = merge(hics[["KO_16h_rep2"]][chr], hics[["KO_16h_rep1"]][chr]),
    WT_16h = merge(hics[["WT_16h_rep2"]][chr], hics[["WT_16h_rep1"]][chr]))
  maps <- imap(merged_replicates, ~ plotMatrix(.x, use.scores = 'balanced', caption = FALSE, limits = c(0, 4)) + ggtitle(.y))
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/16")
  p <- cowplot::plot_grid(plotlist = maps, nrow = 1)
  ggsave(plot = p, filename = paste0(chr, "_KO_16_WT_16.png"))
}


# Plot log scale KO_40 vs WT_40 for differential interactions

for(i in 1:length(chr.list)){
  chr <- chr.list[i]
  
  # get differential results
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/multiHiCcompare/Results/KO_40_vs_WT_40")
  load(paste0(chr,".RData"))
  hic_table <- results
 
  colnames(hic_table)[1:3] <- c("seqnames1", "start1", "start2")
  # prepare data to plot
  gis <- hic_table |> 
    mutate(
      seqnames2 = seqnames1, 
      end1 = start1 + 9999, 
      end2 = start2 + 9999
    ) |> 
    filter(abs(logFC) >= 0.58) |>
    df2gi() 
  
  # extract hic matrix
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
  merged_replicates <- list(
    KO_40h = merge(hics[["KO_40h_rep2"]][chr], hics[["KO_40h_rep1"]][chr]),
    WT_40h = merge(hics[["WT_40h_rep2"]][chr], hics[["WT_40h_rep1"]][chr]))
  
  # save plots
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/differential/KO_40_vs_WT_40")
  p <- cowplot::plot_grid(
    plotMatrix(merged_replicates[['KO_40h']], use.scores = 'balanced', caption = FALSE, limits = c(0, 4)) + ggtitle("KO_40h"),
    plotMatrix(merged_replicates[['WT_40h']], use.scores = 'balanced', caption = FALSE, limits = c(0, 4)) + ggtitle("WT_40h"),
    plotMatrix(gis, use.scores = 'logFC', scale = 'linear', limits = c(-1.5, 1.5), cmap = bgrColors()) + ggtitle("differential interaction"), align = "hv", axis = 'tblr', nrow = 1)
  ggsave(plot = p, filename = paste0(chr, "_KO_40_vs_WT_40.png"), width = 8, height = 6, dpi = 300, units = "in")
  
}


# Plot log scale KO_16 vs WT_16 for differential interactions

for(i in 1:length(chr.list)){
  chr <- chr.list[i]
  
  # get differential results
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/multiHiCcompare/Results/KO_16_vs_WT_16")
  load(paste0(chr,".RData"))
  hic_table <- results
  
  colnames(hic_table)[1:3] <- c("seqnames1", "start1", "start2")
  # prepare data to plot
  gis <- hic_table |> 
    mutate(
      seqnames2 = seqnames1, 
      end1 = start1 + 9999, 
      end2 = start2 + 9999
    ) |> 
    filter(abs(logFC) >= 0.58) |>
    df2gi() 
 
  # extract hic matrix
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
  merged_replicates <- list(
    KO_16h = merge(hics[["KO_16h_rep2"]][chr], hics[["KO_16h_rep1"]][chr]),
    WT_16h = merge(hics[["WT_16h_rep2"]][chr], hics[["WT_16h_rep1"]][chr]))
  
  # save plots
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/differential/KO_16_vs_WT_16")
  p <- cowplot::plot_grid(
    plotMatrix(merged_replicates[['KO_16h']], use.scores = 'balanced', caption = FALSE, limits = c(0, 4)) + ggtitle("KO_16h"),
    plotMatrix(merged_replicates[['WT_16h']], use.scores = 'balanced', caption = FALSE, limits = c(0, 4)) + ggtitle("WT_16h"),
    plotMatrix(gis, use.scores = 'logFC', scale = 'linear', limits = c(-1.5, 1.5), cmap = bgrColors()) + ggtitle("differential interaction"), align = "hv", axis = 'tblr', nrow = 1)
  ggsave(plot = p, filename = paste0(chr, "_KO_16_vs_WT_16.png"), width = 8, height = 6, dpi = 300, units = "in")
  
}