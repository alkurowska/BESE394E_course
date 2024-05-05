#####################
## Quality Metrics ##
#####################

# Set a working directory
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

# Load libraries
library(ggplot2)
library(hicrep)
library(HiContacts)
library(purrr)
library(HiCExperiment)
library(reshape2)

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

#########################
##### Mapping Stats #####
#########################

# Plot Sequenced Reads Stats: Normal Paired, Chimeric Paired, Chimeric Ambiguous, Unmapped
# Normal Paired and Chimeric Paired constitute the Valid Pairs 
stats <- data.frame()
for(i in 1:length(samples)){
  qc <- read.table(paste0(samples[i],"/qc.txt"), header = F, sep = "\t")
  toPlot <- as.data.frame(qc[2:5,1:2])
  toPlot$id <- metadata[i,]$ID
  stats <- rbind(stats, toPlot)
}
colnames(stats) <- c("reads", "value", "sample")

# plot
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")

png("Sequenced_reads_stats.png", width = 200, height = 100, units = "mm", res = 300)
ggplot(stats, aes(fill=reads, y=sample, x=value)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("") + 
  ylab("")
dev.off()

#------------------------
# Plot Valid Pairs Stats: Unique Reads, PCR Duplicates, Optical Duplicates
# Duplicates are Removed
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

stats <- data.frame()
for(i in 1:length(samples)){
  qc <- read.table(paste0(samples[i],"/qc.txt"), header = F, sep = "\t")
  toPlot <- as.data.frame(qc[8:10,1:2])
  toPlot$id <- metadata[i,]$ID
  stats <- rbind(stats, toPlot)
}
colnames(stats) <- c("reads", "value", "sample")

# plot
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")

png("Valid_pairs_stats.png", width = 200, height = 100, units = "mm", res = 300)
ggplot(stats, aes(fill=reads, y=sample, x=value)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("") + 
  ylab("")
dev.off()

#--------------------
# Unique Reads Stats: Intra-fragment Reads (unligated), Below MAPQ Threshold (< 30), Hi-C Contacts
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

stats <- data.frame()
for(i in 1:length(samples)){
  qc <- read.table(paste0(samples[i],"/qc.txt"), header = F, sep = "\t")
  toPlot <- as.data.frame(qc[12:14,1:2])
  toPlot$id <- metadata[i,]$ID
  stats <- rbind(stats, toPlot)
}
colnames(stats) <- c("reads", "value", "sample")

# plot
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")

png("Unique_pairs_stats.png", width = 200, height = 100, units = "mm", res = 300)
ggplot(stats, aes(fill=reads, y=sample, x=value)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("") + 
  ylab("")
dev.off()

#--------------------
# Hi-C contact Stats: Inter-chromosomal, Intra-chromosomal
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

stats <- data.frame()
for(i in 1:length(samples)){
  qc <- read.table(paste0(samples[i],"/qc.txt"), header = F, sep = "\t")
  toPlot <- as.data.frame(qc[16:17,1:2])
  toPlot$id <- metadata[i,]$ID
  stats <- rbind(stats, toPlot)
}
colnames(stats) <- c("reads", "value", "sample")

# plot
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")

png("Hi-C_stats.png", width = 200, height = 100, units = "mm", res = 300)
ggplot(stats, aes(fill=reads, y=sample, x=value)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("") + 
  ylab("")
dev.off()

#-------------------------
# Intra-chromosomal Stats: Short Range (<20Kb), Long Range (>20Kb)
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

stats <- data.frame()
for(i in 1:length(samples)){
  qc <- read.table(paste0(samples[i],"/qc.txt"), header = F, sep = "\t")
  toPlot <- as.data.frame(qc[18:19,1:2])
  toPlot$id <- metadata[i,]$ID
  stats <- rbind(stats, toPlot)
}
colnames(stats) <- c("reads", "value", "sample")

# plot
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")

png("Intra_chrom_stats.png", width = 200, height = 100, units = "mm", res = 300)
ggplot(stats, aes(fill=reads, y=sample, x=value)) + 
  geom_bar(position="fill", stat="identity") +
  xlab("") + 
  ylab("")
dev.off()


###########################
##### Reproducibility #####
###########################

# Set a working directory 
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results/")

# Chromosome list
chr.list <- toupper(c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3"))

# Import Hi-C contact matrices
hics <- list()
for(i in 1:nrow(metadata)){
  hic <- import(paste0(rownames(metadata)[i],"/inter.mcool"))
  hics <- append(hics, hic)
}
names(hics) <- metadata$ID
  
# Compute correlation between experiments
# stratum-adjusted correlations between Hi-C datasets. 
# “Stratum” refers to the distance from the main diagonal: with increase distance from the main diagonal, interactions of the DNA polymer are bound to decrease. hicrep computes a “per-stratum” correlation score and computes a weighted average correlation for entire chromosomes.

for(k in 1:length(chr.list)){
  chr <- chr.list[k]
  sc <- data.frame(matrix(NA, ncol=8, nrow=8))
  colnames(sc) <- metadata$ID
  rownames(sc) <- metadata$ID
  for(i in 1:length(hics)){
    for(j in 1:length(hics)){
      name1 <- names(hics)[i] #select sample 1
      name2 <- names(hics)[j] #select sample 2
      setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results/")
      # Depending on matrix dimensions, compare them in different order. Workaround of https://github.com/TaoYang-dev/hicrep/issues/70
      mat1 <- hics[[name1]][chr] |> as.matrix(sparse = TRUE, use.scores = 'count') #extract matrix for sample 1 for chr k
      mat2 <- hics[[name2]][chr] |> as.matrix(sparse = TRUE, use.scores = 'count') #extract matrix for sample 2 for chr k
      if (nrow(mat1) > nrow(mat2)) {
        scc <- get.scc(mat2, mat1, resol = 10000, h = 2, lbr = 20000, ubr = 200000)$scc #calculate correlation
      } else {
        scc <- get.scc(mat1, mat2, resol = 10000, h = 2, lbr = 20000, ubr = 200000)$scc
      }
      sc[name1,name2] <- scc
    }
  }
  new_results <- setNames(melt(as.matrix(sc)), c('sample1', 'sample2', chr))
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")
  write.table(new_results, sep = "\t", paste0(chr,"_scc.txt"))
}

# Find mean correlation value from all of the chromosomes
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/QC/Results")
scc.final <- data.frame(matrix(NA, nrow = nrow(metadata)*nrow(metadata), ncol = length(chr.list)))
colnames(scc.final) <- chr.list
for(i in 1:length(chr.list)){
  chr <- chr.list[i]
  scc.res <- read.table(paste0(chr.list[i],"_scc.txt"), sep = "\t")
  scc.final[,chr] <- scc.res[,chr]
}
scc.res$sample2==data$sample2
data <- as.data.frame(rowSums(scc.final)/length(chr.list))
colnames(data) <- "scc"
data$sample1 <- metadata$ID
data$sample2 <- c(rep(metadata$ID[1], 8), rep(metadata$ID[2], 8), rep(metadata$ID[3], 8), rep(metadata$ID[4], 8), 
                  rep(metadata$ID[5], 8), rep(metadata$ID[6], 8), rep(metadata$ID[7], 8), rep(metadata$ID[8], 8))

png("Stratum_correlation.png", width = 100, height = 100, units = "mm", res = 300)
ggplot(data, aes(x = sample1, y = sample2, fill = scc)) + 
  geom_tile() + 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  theme_bw() + 
  coord_fixed(ratio = 1) + 
  xlab("") + ylab("") +
  scale_fill_gradientn(colours = bgrColors())
dev.off()

