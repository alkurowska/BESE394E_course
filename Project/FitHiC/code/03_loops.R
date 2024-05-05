##############################
## Singificant Interactions ##
##############################

# Set a working directory
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")

# Load libraries
library(strawr)
library(FitHiC)

# ALL SAMPLES
# SRR19611534     KO_40
# SRR19611535     KO_40
# SRR19611536     WT_40
# SRR19611537     WT_40
# SRR19611538     KO_16
# SRR19611539     KO_16
# SRR19611540     WT_16
# SRR19611541     WT_16


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


############################
###### Prepare inputs ######
############################

######## INTERFILE #########
# This file stores the information about interactions between fragment pairs. 

chr.list <- toupper(c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3","Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3"))
samples.id <- rownames(metadata)

for(i in 1:length(samples.id)){
  intersfile <- data.frame()
  for(j in 1:length(chr.list)){
    chr <- as.character(chr.list[j])
    # Extract interaction file
    setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/hic-results")
    hic.file <- strawr::straw("NONE", paste0(samples.id[i],"/inter_30.hic"), chr, chr, "BP", 10000) # 10k resolution, intra-chromosomal
    
    # INTERSFILE
    inter.file <- matrix(0, ncol = 5, nrow = nrow(hic.file))
    inter.file <- as.data.frame(inter.file)
    colnames(inter.file) <- c("Chromosome1.Name",	"Mid.Point.1",	"Chromosome2.Name",	"Mid.Point.2",	"Hit.Count")
    inter.file$Chromosome1.Name <- chr
    inter.file$Mid.Point.1 <- hic.file$x + 5000
    inter.file$Chromosome2.Name <- chr
    inter.file$Mid.Point.2 <- hic.file$y + 5000
    inter.file$Hit.Count <- hic.file$counts
    
    intersfile <- rbind(intersfile,inter.file)
  }
  
  
  colnames(intersfile) <- c("Chromosome1.Name",	"Mid.Point.1",	"Chromosome2.Name",	"Mid.Point.2",	"Hit.Count")
  newname <- paste0("intersfile.",samples.id[i])
  assign(newname,intersfile)
  setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/FitHiC/Input")
  write.table(intersfile, paste0(newname, ".bed"), sep = "\t", row.names = F, col.names = F)
}

######## FRAGSFILE #########
# This file stores the information about midpoints (or start indices) of the fragments 

setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/FitHiC/Input")

files <- list(intersfile.SRR19611534, intersfile.SRR19611535, intersfile.SRR19611536, intersfile.SRR19611537,
              intersfile.SRR19611538, intersfile.SRR19611539, intersfile.SRR19611540, intersfile.SRR19611541)

names(files) <- c("intersfile.SRR19611534", "intersfile.SRR19611535", "intersfile.SRR19611536", "intersfile.SRR19611537",
                  "intersfile.SRR19611538", "intersfile.SRR19611539", "intersfile.SRR19611540", "intersfile.SRR19611541")

for(i in 1:length(files)){
  # Combine all fragments 
  file <- files[[i]]
  region1 <- paste0(file$Chromosome1.Name,".", file$Mid.Point.1)
  region2 <- paste0(file$Chromosome2.Name,".", file$Mid.Point.2)
  
  all.regions <- unique(c(region1, region2))
  regions <- data.frame(do.call("rbind", strsplit(as.character(all.regions), ".", fixed = TRUE)))
  
  fragsfile <- as.data.frame(matrix(0, ncol = 5, nrow = nrow(regions)))
  colnames(fragsfile) <- c("Chromosome.Name",	"Column.2", "Mid.Point",	"Hit.Count", "Column.5")
  fragsfile$Chromosome.Name <- regions$X1
  fragsfile$Mid.Point <- regions$X2
  fragsfile$Hit.Count <- 1
  
  newname <- gsub("inters", "frags", names(files)[i])
  assign(newname,fragsfile)
  write.table(fragsfile, paste0(newname, ".bed"), sep = "\t", row.names = F, col.names = F)
}


############################
######## RUN FitHiC ########
############################

# Set a working directory
setwd("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/FitHiC/Input")

for( i in 1:length(samples.id)){
  fragsfile <- paste0("fragsfile.",samples.id[i],".bed.gz")
  intersfile <- paste0("intersfile.",samples.id[i],".bed.gz")
  outdir <- paste0("/Users/kurowsaa/OneDrive/Documents/KAUST/BESE394E_homework/BESE394E_course/FINAL/FitHiC/Results/",samples.id[i])
  FitHiC(fragsfile, intersfile, outdir, distUpThres = 5000000, distLowThres = 50000,
         libname=paste0(samples.id[i],".",metadata$ID[i]), noOfBins=100,
         visual=TRUE)
}





