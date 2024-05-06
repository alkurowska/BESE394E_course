#check if var genes are within our differential interactions

###########################################################
###             run in bash
### awk 'FNR>1' *DIRs.txt > WT_40_vs_WT_16DIRs.txt
### awk 'FNR>1' *DIRs.txt > KO_40_vs_KO_16DIRs.txt
### awk 'FNR>1' *DIRs.txt > KO_40_vs_WT_40DIRs.txt
### awk 'FNR>1' *DIRs.txt > KO_16_vs_WT_16DIRs.txt
###########################################################

setwd("../multiHiCcompare/")

#libraries
library(dplyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(scales)
library(topGO)
library(AnnotationHub)
library(biomaRt)
library(scales)
set.seed(123456)

chr.df <- data.frame("seqnames"=c("PF3D7_01_V3","PF3D7_02_V3","PF3D7_03_V3","PF3D7_04_V3","PF3D7_05_V3","PF3D7_06_V3","PF3D7_07_V3","PF3D7_08_V3","PF3D7_09_V3","PF3D7_10_V3","PF3D7_11_V3","PF3D7_12_V3","PF3D7_13_V3","PF3D7_14_V3"),
                      "chr"= paste("chr",1:14, sep = ""))

#functions
annotateIRs2<-function(queryGRanges, outFileName){
  #find overlaps with annotation - not only proteins5
  olaps <- findOverlaps(queryGRanges, PlasmoDB_Annot, type = "any", select = "all") #all features
  olaps_prot <- findOverlaps(queryGRanges, PlasmoDB_Annot_Prot, type = "any", select = "all") #only protein coding genes
  
  #detect features and plot
  features <- PlasmoDB_Annot@elementMetadata$feature[queryHits(olaps)]
  pdf(paste("Results/Annotations/",outFileName,".pdf",sep = ""))
  print({
    barplot(table(features), cex.names = 0.5, las=2, main = outFileName)
  })
  dev.off()
  
  #generate GRanges with annotation for each queryHit
  genes_olap <- olaps_prot %>% as.data.frame %>% group_by(queryHits) %>%
    mutate(genes = PlasmoDB_Annot@elementMetadata$attribute[subjectHits]) %>%
    dplyr::select(queryHits, genes)
  
  GRanges_annotatedRegions <- queryGRanges[queryHits(olaps_prot)] 
  mcols(GRanges_annotatedRegions)$Description <- genes_olap$genes
  
  #save to file
  write.table(GRanges_annotatedRegions, file=paste("Results/Annotations/",outFileName, "annotation.GRanges", sep = ""))
  
  return(GRanges_annotatedRegions)
}

ORA2 <- function(geneList, nameResults){
  GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO, geneSel=selection)
  
  results <- runTest(GOdata, algorithm="weight01", statistic="ks")
  goEnrichment <- GenTable(GOdata, KS=results, orderBy="KS", topNodes=length(usedGO(GOdata)))
  goEnrichment$KS <- as.numeric(goEnrichment$KS)
  goEnrichment <- goEnrichment[goEnrichment$KS<0.9,]
  goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
  goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
  goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
  goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))
  
  ## To also get Gene IDs
  get.genes <- plyr::ldply(genesInTerm(GOdata,goEnrichment$GO.ID), rbind)
  tidy.genes <- get.genes %>% tidyr::unite(col = "col", 2:ncol(get.genes), sep = ",", na.rm = T)
  goEnrichment$associated_genes <- tidy.genes$col
  
  goEnrichment$Term<-pRoloc::goIdToTerm(goEnrichment$GO.ID, names = TRUE, keepNA = TRUE)
  goEnrichment$Term<-factor(goEnrichment$Term, levels = c(goEnrichment$Term[order(goEnrichment$KS, decreasing = TRUE)]))
  
  pdf(paste("Results/ORA/",nameResults,".pdf",sep = ""), height = 10)
  print({
    dot(goEnrichment[1:20,])
  })
  dev.off()
  
  a <- goEnrichment
  
  writexl::write_xlsx(list("Up"=a ),paste("Results/ORA/",nameResults,".xlsx",sep = ""))
}

dot <- function(res){
  ggplot(res,
         aes(x = Term, y = -log10(KS), size = -log10(KS), fill = -log10(KS))) +
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'royalblue', high = 'red4') +
    xlab('') + ylab('Enrichment score') +
    labs(
      title = 'GO Biological processes',
      subtitle = paste('Top',length(res$Term),'terms ordered by Kolmogorov-Smirnov p-value'),
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 11, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 10, face = 'bold', vjust = 1),
      axis.text.x = element_text(angle = 0, size = 10, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 10, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 10, face = 'bold'),
      axis.title.x = element_text(size = 10, face = 'bold'),
      axis.title.y = element_text(size = 10, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      #Legend
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 11, face = "bold"), # Text size
      title = element_text(size = 11, face = "bold")) +
    coord_flip()+scale_x_discrete(labels = label_wrap(40))
}

get_genes2 <- function(GRangesDF){
  genes <- c(GRangesDF$Description[grepl("gene_id=PF3D7_", GRangesDF$Description)])
  qvals <- c(GRangesDF$p.adj[grepl("gene_id=PF3D7_", GRangesDF$Description)])
  genes <- str_split(genes,";")
  genes <- lapply(genes, '[[',3)
  genes <- unlist(genes, recursive=FALSE)
  genes <- substring(genes, 9, 21)
  idx <- !duplicated(genes)
  genes_qvals <- qvals[idx]
  names(genes_qvals) <- genes[idx]
  return(genes_qvals)
}

selection <- function(allScore){ 
  return(allScore < 0.05)}

###########################################################
KO_40_vs_WT_40 <- read.table("KO_40_vs_WT_40/KO_40_vs_WT_40DIRs.txt")[2:12]
colnames(KO_40_vs_WT_40) <- c("chr1","start1","end1","chr2","start2","end2","D","logFC","logCPM","p.value","p.adj")

#makeGRanges
KO_40_vs_WT_40_GR <- makeGRangesFromDataFrame(data.frame(
  "Chr"= chr.df$seqnames[match(c(KO_40_vs_WT_40$chr1, KO_40_vs_WT_40$chr2),chr.df$chr)],
  "Start"= c(KO_40_vs_WT_40$start1, KO_40_vs_WT_40$start2),
  "End"= c(KO_40_vs_WT_40$end1, KO_40_vs_WT_40$end2), 
  "logFC" = rep(KO_40_vs_WT_40$logFC, 2),
  "p.adj" = rep(KO_40_vs_WT_40$p.adj, 2)),
  keep.extra.columns = TRUE)

#annotation
KO_40_vs_WT_40_GR_Ann <- annotateIRs2(KO_40_vs_WT_40_GR, "KO_40_vs_WT_40")

#ORA
KO_40_vs_WT_40_genes <- get_genes2(KO_40_vs_WT_40_GR_Ann) #get genes ith associated p.adj value
#ORA2(geneList = KO_40_vs_WT_40_genes, nameResults = "KO_40_vs_WT_40")
#no significant terms found


###########################################################
KO_16_vs_WT_16 <- read.table("KO_16_vs_WT_16/KO_16_vs_WT_16DIRs.txt")[2:12]
colnames(KO_16_vs_WT_16) <- c("chr1","start1","end1","chr2","start2","end2","D","logFC","logCPM","p.value","p.adj")

#makeGRanges
KO_16_vs_WT_16_GR <- makeGRangesFromDataFrame(data.frame(
  "Chr"= chr.df$seqnames[match(c(KO_16_vs_WT_16$chr1, KO_16_vs_WT_16$chr2),chr.df$chr)],
  "Start"= c(KO_16_vs_WT_16$start1, KO_16_vs_WT_16$start2),
  "End"= c(KO_16_vs_WT_16$end1, KO_16_vs_WT_16$end2), 
  "logFC" = rep(KO_16_vs_WT_16$logFC, 2),
  "p.adj" = rep(KO_16_vs_WT_16$p.adj, 2)),
  keep.extra.columns = TRUE)

#annotation
KO_16_vs_WT_16_GR_Ann <- annotateIRs2(KO_16_vs_WT_16_GR, "KO_16_vs_WT_16")

#ORA
KO_16_vs_WT_16_genes <- get_genes2(KO_16_vs_WT_16_GR_Ann) #get genes ith associated p.adj value
#ORA2(geneList = KO_16_vs_WT_16_genes, nameResults = "KO_16_vs_WT_16")
#no significat terms were found (KS < 0.05)

###########################################################
#check if var genes are within DIRs
var_genes <- gff[grepl(pattern = "VAR",x = gff$attribute),]
tmp <- var_genes$attribute[grepl("ID=PF3D7_", var_genes$attribute)]
tmp <- str_split(tmp,";")
tmp <- lapply(tmp, '[[',1)
tmp <- unlist(tmp, recursive=FALSE)
tmp <- substring(tmp, 4, 16)
var_genes$GeneID <- tmp
var_genes_GR <- makeGRangesFromDataFrame(var_genes, keep.extra.columns = TRUE)

#overlap with KO_16_vs_WT_16
var_KO_16_vs_WT_16 <- subsetByOverlaps(var_genes_GR, KO_16_vs_WT_16_GR)
write.table(var_KO_16_vs_WT_16, file=paste("Results/VAR_DIRs_overlap_KO_16_vs_WT_16.GRanges", sep = ""))
#at 16hrs 14 DIRs overlap with VAR genes

#overlap with KO_40_vs_WT_40
var_KO_40_vs_WT_40 <- subsetByOverlaps(var_genes_GR, KO_40_vs_WT_40_GR)
write.table(var_KO_40_vs_WT_40, file=paste("Results/VAR_DIRs_overlap_KO_40_vs_WT_40.GRanges", sep = ""))
#at 40 hrs no DIRs overlap with a VAR gene
