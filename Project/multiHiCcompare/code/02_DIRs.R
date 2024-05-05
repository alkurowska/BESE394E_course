########################################
## Differentially Interacting Regions ##
########################################

# Load libraries 
library(multiHiCcompare)
library(strawr)
library(edgeR)

# Set a working directory
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

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
metadata$group <- paste0(metadata$strain, "_", metadata$h.p.i)

# list of sampels to analyze
samples <- rownames(metadata)

# Extract raw martices for every chromosome per sample

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

chr.list <- toupper(c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3"))

inter.list <- list() #create a list of all interactions per chr per sample
for(i in 1:length(samples)){
  sample.id <- samples[i] #define the sample id
  for(j in 1: length(chr.list)){
    chr <- as.character(chr.list[j]) #defined the chromosome for intra-chromosomal interactions 
    hic.file <- strawr::straw("NONE", paste0(sample.id,"/inter_30.hic"), chr, chr, "BP", 10000) #extract raw matrix, at 10k resolution 
    chr.name <- as.data.frame(rep(j, nrow(hic.file))) #add column with a chr name
    hic.file <- cbind(chr.name, hic.file)
    colnames(hic.file) <- c("chr", "region1", "region2", "IF") ## matrix ready to be input into multiHiCcompare
    newname <- paste0(sample.id,"_",chr)
    assign(newname, hic.file)
    hic.file <- list(hic.file)
    names(hic.file) <- newname
    
    inter.list <- append(inter.list, hic.file) # append to the list

  }
}


# ALL SAMPLES
# SRR19611534     KO_40
# SRR19611535     KO_40
# SRR19611536     WT_40
# SRR19611537     WT_40
# SRR19611538     KO_16
# SRR19611539     KO_16
# SRR19611540     WT_16
# SRR19611541     WT_16

#########################
##### KO_40 vs WT_40 ####
#########################
# SRR19611534     KO_40
# SRR19611535     KO_40
# SRR19611536     WT_40
# SRR19611537     WT_40

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/multiHiCcompare/Results/KO_40_vs_WT_40")

# all interactions per chr per samples

for(i in 1:length(chr.list)){
  # make groups & covariate input
  groups <- factor(metadata[metadata$group%in%c("KO_40", "WT_40"),]$strain)
  covariates <- data.frame(h.p.i = factor(metadata[metadata$group%in%c("KO_40", "WT_40"),]$h.p.i))
  # create hicexp object
  samples.id <- rownames(metadata[metadata$group%in%c("KO_40", "WT_40"),])
  data <- paste0(samples.id,"_",chr.list[i])
  hicexp <- make_hicexp(data_list = inter.list[data],
                         groups = groups, #covariates = covariates, #zero.p = 0, A.min = 5, 
                         filter = FALSE)
  
  # make MD plot before norm
  png(paste0(chr.list[i],".png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Library normalization FASTLO
  # Joint normalization algorithm for all of the samples
  hicexp <- fastlo(hicexp, verbose = FALSE, parallel = FALSE)
  
  # make MD plot
  png(paste0(chr.list[i],"_fastlo.png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Differentially Interacting Regions
  # GLM methods
  
  # Create design matrix
  design.matrix <- model.matrix(~ 0 + groups)
  colnames(design.matrix) <- gsub("groups", "", colnames(design.matrix))
  
  ## KO-WT in 40 h.p.i.
  # Specify contrasts
  contr <- makeContrasts(KO_40_vs_WT_40 = KO-WT, levels=colnames(design.matrix))
  
  hicexp <- hic_glm(hicexp, design = design.matrix, contrast = contr,
                                   method = "QLFTest", p.method = "fdr", parallel = FALSE)
  # save results
  results <- results(hicexp)
  save(results, file = paste0(chr.list[i],".RData"))
  
  # Get top DIRs
  td <- topDirs(hicexp, logfc_cutoff = 0.58, logcpm_cutoff = 0.58, pval_aggregate = "max", p.adj_cutoff = 0.05, return_df = 'pairedbed')
  write.table(td, paste0(chr.list[i],"_DIRs.txt"), sep = "\t")
  
  # Plot a composite MD plot of the results of a comparison where the significant differences are highlighted.
  png(paste0(chr.list[i],"_composite.png"), width = 100, height = 100, units = "mm", res = 300)
  MD_composite(hicexp)
  dev.off()
}



#########################
##### KO_16 vs WT_16 ####
#########################
# SRR19611538     KO_16
# SRR19611539     KO_16
# SRR19611540     WT_16
# SRR19611541     WT_16

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/multiHiCcompare/Results/KO_16_vs_WT_16")

# all interactions per chr per samples

for(i in 1:length(chr.list)){
  # make groups & covariate input
  groups <- factor(metadata[metadata$group%in%c("KO_16", "WT_16"),]$strain)
  covariates <- data.frame(h.p.i = factor(metadata[metadata$group%in%c("KO_16", "WT_16"),]$h.p.i))
  # create hicexp object
  samples.id <- rownames(metadata[metadata$group%in%c("KO_16", "WT_16"),])
  data <- paste0(samples.id,"_",chr.list[i])
  hicexp <- make_hicexp(data_list = inter.list[data],
                        groups = groups, #covariates = covariates, #zero.p = 0, A.min = 5, 
                        filter = FALSE)
  
  # make MD plot before norm
  png(paste0(chr.list[i],".png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Library normalization FASTLO
  # Joint normalization algorithm for all of the samples
  hicexp <- fastlo(hicexp, verbose = FALSE, parallel = FALSE)
  
  # make MD plot
  png(paste0(chr.list[i],"_fastlo.png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Differentially Interacting Regions
  # GLM methods
  
  # Create design matrix
  design.matrix <- model.matrix(~ 0 + groups)
  colnames(design.matrix) <- gsub("groups", "", colnames(design.matrix))
  
  ## KO-WT in 16 h.p.i.
  # Specify contrasts
  contr <- makeContrasts(KO_16_vs_WT_16 = KO-WT, levels=colnames(design.matrix))
  
  hicexp <- hic_glm(hicexp, design = design.matrix, contrast = contr,
                    method = "QLFTest", p.method = "fdr", parallel = FALSE)
  
  # save results
  results <- results(hicexp)
  save(results, file = paste0(chr.list[i],".RData"))
  
  # Get top DIRs
  td <- topDirs(hicexp, logfc_cutoff = 0.58, logcpm_cutoff = 0.58, pval_aggregate = "max", p.adj_cutoff = 0.05, return_df = 'pairedbed')
  write.table(td, paste0(chr.list[i],"_DIRs.txt"), sep = "\t")
  
  # Plot a composite MD plot of the results of a comparison where the significant differences are highlighted.
  png(paste0(chr.list[i],"_composite.png"), width = 100, height = 100, units = "mm", res = 300)
  MD_composite(hicexp)
  dev.off()
}


#########################
##### KO_40 vs KO_16 ####
#########################
# SRR19611534     KO_40
# SRR19611535     KO_40
# SRR19611538     KO_16
# SRR19611539     KO_16

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/multiHiCcompare/Results/KO_40_vs_KO_16")

# all interactions per chr per samples

for(i in 1:length(chr.list)){
  # make groups & covariate input
  groups <- factor(metadata[metadata$group%in%c("KO_40", "KO_16"),]$group)
  # create hicexp object
  samples.id <- rownames(metadata[metadata$group%in%c("KO_40", "KO_16"),])
  data <- paste0(samples.id,"_",chr.list[i])
  hicexp <- make_hicexp(data_list = inter.list[data],
                        groups = groups, #zero.p = 0, A.min = 5, 
                        filter = FALSE)
  
  # make MD plot before norm
  png(paste0(chr.list[i],".png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Library normalization FASTLO
  # Joint normalization algorithm for all of the samples
  hicexp <- fastlo(hicexp, verbose = FALSE, parallel = FALSE)
  
  # make MD plot
  png(paste0(chr.list[i],"_fastlo.png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Differentially Interacting Regions
  # GLM methods
  
  # Create design matrix
  design.matrix <- model.matrix(~ 0 + groups)
  colnames(design.matrix) <- gsub("groups", "", colnames(design.matrix))
  
  ## 40 h.p.i. vs 16 h.p.i. in KO
  # Specify contrasts
  contr <- makeContrasts(KO_40_vs_KO_16 = KO_40-KO_16, levels=colnames(design.matrix))
  
  hicexp <- hic_glm(hicexp, design = design.matrix, contrast = contr,
                    method = "QLFTest", p.method = "fdr", parallel = FALSE)
  
  # save results
  results <- results(hicexp)
  save(results, file = paste0(chr.list[i],".RData"))
  
  # Get top DIRs
  td <- topDirs(hicexp, logfc_cutoff = 0.58, logcpm_cutoff = 0.58, pval_aggregate = "max", p.adj_cutoff = 0.05, return_df = 'pairedbed')
  write.table(td, paste0(chr.list[i],"_DIRs.txt"), sep = "\t")
  
  # Plot a composite MD plot of the results of a comparison where the significant differences are highlighted.
  png(paste0(chr.list[i],"_composite.png"), width = 100, height = 100, units = "mm", res = 300)
  MD_composite(hicexp)
  dev.off()
}



#########################
##### WT_40 vs WT_16 ####
#########################
# SRR19611536     WT_40
# SRR19611537     WT_40
# SRR19611540     WT_16
# SRR19611541     WT_16

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/multiHiCcompare/Results/WT_40_vs_WT_16")

# all interactions per chr per samples

for(i in 1:length(chr.list)){
  # make groups & covariate input
  groups <- factor(metadata[metadata$group%in%c("WT_40", "WT_16"),]$group)
  # create hicexp object
  samples.id <- rownames(metadata[metadata$group%in%c("WT_40", "WT_16"),])
  data <- paste0(samples.id,"_",chr.list[i])
  hicexp <- make_hicexp(data_list = inter.list[data],
                        groups = groups, #zero.p = 0, A.min = 5, 
                        filter = FALSE)
  
  # make MD plot before norm
  png(paste0(chr.list[i],".png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Library normalization FASTLO
  # Joint normalization algorithm for all of the samples
  hicexp <- fastlo(hicexp, verbose = FALSE, parallel = FALSE)
  
  # make MD plot
  png(paste0(chr.list[i],"_fastlo.png"), width = 300, height = 200, units = "mm", res = 300)
  MD_hicexp(hicexp, prow = 2, pcol = 3, plot.loess = T)
  dev.off()
  
  # Differentially Interacting Regions
  # GLM methods
  
  # Create design matrix
  design.matrix <- model.matrix(~ 0 + groups)
  colnames(design.matrix) <- gsub("groups", "", colnames(design.matrix))
  
  ## 40 h.p.i. vs 16 h.p.i. in WT
  # Specify contrasts
  contr <- makeContrasts(WT_40_vs_WT_16 = WT_40-WT_16, levels=colnames(design.matrix))
  
  hicexp <- hic_glm(hicexp, design = design.matrix, contrast = contr,
                    method = "QLFTest", p.method = "fdr", parallel = FALSE)
  
  # save results
  results <- results(hicexp)
  save(results, file = paste0(chr.list[i],".RData"))
  
  # Get top DIRs
  td <- topDirs(hicexp, logfc_cutoff = 0.58, logcpm_cutoff = 0.58, pval_aggregate = "max", p.adj_cutoff = 0.05, return_df = 'pairedbed')
  write.table(td, paste0(chr.list[i],"_DIRs.txt"), sep = "\t")
  
  # Plot a composite MD plot of the results of a comparison where the significant differences are highlighted.
  png(paste0(chr.list[i],"_composite.png"), width = 100, height = 100, units = "mm", res = 300)
  MD_composite(hicexp)
  dev.off()
}

