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
library(InteractionSet)

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

########################
######### Loops ########
########################

# Get significant interactions info

for(i in 1:nrow(metadata)){
  sample <- rownames(metadata)[i]
  setwd(paste0("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/FitHiC/Results/",sample))

  fithic <- read.table(paste0(sample,".",metadata$ID[i],".spline_pass2.significances.txt.gz"), sep = "\t", header = T)
  fithic$start1 <- fithic$fragmentMid1 - 4999
  fithic$end1 <- fithic$fragmentMid1 + 5000
  fithic$start2 <- fithic$fragmentMid2 - 4999
  fithic$end2 <- fithic$fragmentMid2 + 5000

  cutoff <- quantile(fithic$contactCount, 0.90)
  fithic <- fithic[fithic$contactCount > cutoff,]
  
  # filter significant interactions fdr < 0.001
  fithic <- fithic[fithic$q_value < 0.001,]

  index.1 <- c(1,8,9)
  index.2 <- c(3,10,11)
  
  region.1 <- fithic[index.1]
  region.2 <- fithic[index.2]
  
  colnames(region.1) <- c("chr", "start", "end")
  colnames(region.2) <- c("chr", "start", "end")
  
  region.1 <- makeGRangesFromDataFrame(region.1)
  region.2 <- makeGRangesFromDataFrame(region.2)
  
  # Get interactions  
  gi <- GInteractions(region.1, region.2)
  gi$counts <- fithic$contactCount
  gi$qvalue <- fithic$q_value
  gi$pvalue <- fithic$p_value
  
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
  topologicalFeatures(hics[[metadata$ID[i]]], "loops") <- gi

}

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

# Merge replicates for ploting
merged_replicates <- list(
  KO_40h = merge(hics[["KO_40h_rep2"]]["PF3D7_11_V3"], hics[["KO_40h_rep1"]]["PF3D7_11_V3"]),
  WT_40h = merge(hics[["WT_40h_rep2"]]["PF3D7_11_V3"], hics[["WT_40h_rep1"]]["PF3D7_11_V3"]),
  KO_16h = merge(hics[["KO_16h_rep2"]]["PF3D7_11_V3"], hics[["KO_16h_rep1"]]["PF3D7_11_V3"]),
  WT_16h = merge(hics[["WT_16h_rep2"]]["PF3D7_11_V3"], hics[["WT_16h_rep1"]]["PF3D7_11_V3"]))

# find overlaping loops between replicates 

rep2 <- topologicalFeatures(hics[["KO_40h_rep2"]], "loops")
rep1 <- topologicalFeatures(hics[["KO_40h_rep1"]], "loops")

loops_KO_40 <- subsetByOverlaps(rep1, rep2, type = "start", use.region="same")

rep2 <- topologicalFeatures(hics[["KO_16h_rep2"]], "loops")
rep1 <- topologicalFeatures(hics[["KO_16h_rep1"]], "loops")

loops_KO_16 <- subsetByOverlaps(rep1, rep2, type = "start", use.region="same")

rep2 <- topologicalFeatures(hics[["WT_40h_rep2"]], "loops")
rep1 <- topologicalFeatures(hics[["WT_40h_rep1"]], "loops")

loops_WT_40 <- subsetByOverlaps(rep1, rep2, type = "start", use.region="same")

rep2 <- topologicalFeatures(hics[["WT_16h_rep2"]], "loops")
rep1 <- topologicalFeatures(hics[["WT_16h_rep1"]], "loops")

loops_WT_16 <- subsetByOverlaps(rep1, rep2, type = "start", use.region="same")


# Zoom in to the region of interest
# save plots
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
toPlot_KO <- zoom(merged_replicates[["KO_40h"]], 10000) |> refocus('PF3D7_11_V3:100000-500000')
toPlot_WT <- zoom(merged_replicates[["WT_40h"]], 10000) |> refocus('PF3D7_11_V3:100000-500000')


p <- cowplot::plot_grid(
    plotMatrix(toPlot_KO, use.scores = 'balanced', loops = loops_KO_40, limits = c(0, 5), caption = FALSE) + ggtitle("KO_40h"),
    plotMatrix(toPlot_WT, use.scores = 'balanced', loops = loops_WT_40, limits = c(0, 5), caption = FALSE) + ggtitle("WT_40h"))

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/loops")
ggsave(plot = p, filename = "chr11_KO_40_vs_WT_40.png", width = 12, height = 8, dpi = 300, units = "in")


# save plots
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
toPlot_KO <- zoom(merged_replicates[["KO_16h"]], 10000) |> refocus('PF3D7_11_V3:100000-500000')
toPlot_WT <- zoom(merged_replicates[["WT_16h"]], 10000) |> refocus('PF3D7_11_V3:100000-500000')


p <- cowplot::plot_grid(
  plotMatrix(toPlot_KO, use.scores = 'balanced', loops = loops_KO_16, limits = c(0, 5), caption = FALSE) + ggtitle("KO_16h"),
  plotMatrix(toPlot_WT, use.scores = 'balanced', loops = loops_WT_16, limits = c(0, 5), caption = FALSE) + ggtitle("WT_16h"))

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/loops")
ggsave(plot = p, filename = "chr11_KO_16_vs_WT_16.png", width = 12, height = 8, dpi = 300, units = "in")








setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

toPlot <- zoom(hics[["KO_16h_rep2"]], 10000) |> 
  refocus('PF3D7_11_V3:100000-500000') 

p <- plotMatrix(
  toPlot,
  loops = topologicalFeatures(hics[["KO_16h_rep2"]], "loops"),
  limits = c(0, 5),
  caption = FALSE
)

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/loops")

ggsave(plot = p, filename = "chr11_KO_16_rep2.png", width = 6, height = 6, dpi = 300, units = "in")



setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

toPlot <- zoom(hics[["KO_16h_rep1"]], 10000) |> 
  refocus('PF3D7_11_V3:100000-500000') 

p <- plotMatrix(
  toPlot,
  loops = topologicalFeatures(hics[["KO_16h_rep1"]], "loops"),
  limits = c(0, 5),
  caption = FALSE
)

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/plots/loops")

ggsave(plot = p, filename = "chr11_KO_16_rep1.png", width = 6, height = 6, dpi = 300, units = "in")
