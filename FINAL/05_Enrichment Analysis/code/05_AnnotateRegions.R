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

setwd("../FitHiC/")

#relevant chromosomes
chr.list <- c("PF3D7_01_V3","PF3D7_02_V3","PF3D7_03_V3","PF3D7_04_V3","PF3D7_05_V3","PF3D7_06_V3","PF3D7_07_V3","PF3D7_08_V3","PF3D7_09_V3","PF3D7_10_V3","PF3D7_11_V3","PF3D7_12_V3","PF3D7_13_V3","PF3D7_14_V3")

#Remove api and mit chromosome and replace chromosome name: start from PF into PF and ending from v3 into V3                             
gff <- read.delim("../PFalciparum_3D7_Annotation/PlasmoDB/PlasmoDB-67_PFalciparum3D7.gff", header=F, comment.char="#", col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"))
gff$seqname <- str_replace(gff$seqname , "Pf", "PF")
gff$seqname <- str_replace(gff$seqname , "v3", "V3")
gff_flitered <- gff[gff$seqname %in% chr.list,]


#GRanges object with annotation on P.falciparum
PlasmoDB_Annot <- makeGRangesFromDataFrame(gff_flitered,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=chr.list,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"))

gff_flitered_prot <- gff_flitered[gff_flitered$feature == "protein_coding_gene",]

PlasmoDB_Annot_Prot <- makeGRangesFromDataFrame(gff_flitered_prot,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=chr.list,
                                           seqnames.field=c("seqnames", "seqname",
                                                            "chromosome", "chrom",
                                                            "chr", "chromosome_name",
                                                            "seqid"))

#####
#list of FitHiC results
files <- list.files(path = "Data/",pattern = ".gz")

#filter by p.adj < 0.001, counts, and merge replicates
for (i in c(1,3,5,7)) {
  fithic1 <- read.table(gzfile(files[i]), header = 1)
  fithic2 <- read.table(gzfile(files[i+1]), header = 1)
  file_name2 <- str_split(files[i+1], "\\.", simplify = TRUE)
  #filter by contact counts
  cutoff1 <- quantile(fithic1$contactCount, 0.90)
  fithic1 <- fithic1[fithic1$contactCount > cutoff1,]
  cutoff2 <- quantile(fithic2$contactCount, 0.90)
  fithic2 <- fithic2[fithic2$contactCount > cutoff2,]
  #filter by p.adj
  fithic1 <- filter(fithic1, q_value < 0.001)
  fithic2 <- filter(fithic2, q_value < 0.001)
  #keep only contacts in both replicates
  fithic_filtered <- merge(fithic1,fithic2, by = c("fragmentMid1" = "fragmentMid1", "fragmentMid2" = "fragmentMid2", "chr1"="chr1", "chr2"="chr2"))
  #write to file
  write.table(fithic_filtered, file=paste("Data/Filtered/",file_name2[2],"_rep2_filtered.txt", sep = ""))
}

rm(list = c("fithic1", "fithic2","file_name2", "i", "cutoff1", "cutoff2", "fithic_filtered"))

#####
#annotate each anchor
KO_40hrs <- read.table("Data/Filtered/KO_40h_rep1_rep2_filtered.txt")

KO_40_anchor1 <- makeGRangesFromDataFrame(data.frame(
                         "Chr"= KO_40hrs$chr1,
                         "Start"= KO_40hrs$fragmentMid1-4999,
                         "End"= KO_40hrs$fragmentMid1+5000, 
                         "q.value" = KO_40hrs$q_value.x), 
                         keep.extra.columns=TRUE)

KO_40_anchor2 <- makeGRangesFromDataFrame(data.frame(
                                                     "Chr"= KO_40hrs$chr2,
                                                     "Start"= KO_40hrs$fragmentMid2-4999,
                                                     "End"= KO_40hrs$fragmentMid2+5000, 
                                                     "q.value" = KO_40hrs$q_value.x),
                                      keep.extra.columns=TRUE)

WT_40hrs <- read.table("Data/Filtered/WT_40h_rep1_rep2_filtered.txt")

WT_40_anchor1 <- makeGRangesFromDataFrame(data.frame("ID"= paste("Region", 1:dim(WT_40hrs)[1]),
                                                     "Chr"= WT_40hrs$chr1,
                                                     "Start"= WT_40hrs$fragmentMid1-4999,
                                                     "End"= WT_40hrs$fragmentMid1+5000, 
                                                     "q.value" = WT_40hrs$q_value.x),
                                          keep.extra.columns=TRUE)


WT_40_anchor2 <- makeGRangesFromDataFrame(data.frame("ID"= paste("Region", 1:dim(WT_40hrs)[1]),
                                                     "Chr"= WT_40hrs$chr1,
                                                     "Start"= WT_40hrs$fragmentMid2-4999,
                                                     "End"= WT_40hrs$fragmentMid2+5000, 
                                                     "q.value" = WT_40hrs$q_value.x),
                                          keep.extra.columns=TRUE)

KO_16hrs <- read.table("Data/Filtered/KO_16h_rep1_rep2_filtered.txt")

KO_16_anchor1 <- makeGRangesFromDataFrame(data.frame("ID"= paste("Region", 1:dim(KO_16hrs)[1]),
                                                     "Chr"= KO_16hrs$chr1,
                                                     "Start"= KO_16hrs$fragmentMid1-4999,
                                                     "End"= KO_16hrs$fragmentMid1+5000, 
                                                     "q.value" = KO_16hrs$q_value.x),
                                          keep.extra.columns=TRUE)

KO_16_anchor2 <- makeGRangesFromDataFrame(data.frame("ID"= paste("Region", 1:dim(KO_16hrs)[1]),
                                                     "Chr"= KO_16hrs$chr1,
                                                     "Start"= KO_16hrs$fragmentMid2-4999,
                                                     "End"= KO_16hrs$fragmentMid2+5000, 
                                                     "q.value" = KO_16hrs$q_value.x),
                                          keep.extra.columns=TRUE)

WT_16hrs <- read.table("Data/Filtered/WT_16h_rep1_rep2_filtered.txt")

WT_16_anchor1 <- makeGRangesFromDataFrame(data.frame("ID"= paste("Region", 1:dim(WT_16hrs)[1]),
                                                     "Chr"= WT_16hrs$chr1,
                                                     "Start"= WT_16hrs$fragmentMid1-4999,
                                                     "End"= WT_16hrs$fragmentMid1+5000, 
                                                     "q.value" = WT_16hrs$q_value.x),
                                          keep.extra.columns=TRUE)

WT_16_anchor2 <- makeGRangesFromDataFrame(data.frame("ID"= paste("Region", 1:dim(WT_16hrs)[1]),
                                                     "Chr"= WT_16hrs$chr1,
                                                     "Start"= WT_16hrs$fragmentMid2-4999,
                                                     "End"= WT_16hrs$fragmentMid2+5000, 
                                                     "q.value" = WT_16hrs$q_value.x),
                                          keep.extra.columns=TRUE)
#####
annotateIRs<-function(queryGRanges, outFileName){
  #find overlaps with annotation - not only proteins5
  olaps <- findOverlaps(queryGRanges, PlasmoDB_Annot, type = "any", select = "all") #all features
  olaps_prot <- findOverlaps(queryGRanges, PlasmoDB_Annot_Prot, type = "any", select = "all") #only protein coding genes
  
  #detect features and plot
  features <- PlasmoDB_Annot@elementMetadata$feature[queryHits(olaps)]
  pdf(paste("Data/Filtered/Annotations/",outFileName,".pdf",sep = ""))
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
  write.table(GRanges_annotatedRegions, file=paste("Data/Filtered/Annotations/",outFileName, "annotation.GRanges", sep = ""))
  
  return(GRanges_annotatedRegions)
}

KO_40hrs_Anchor1_Ann <- annotateIRs(KO_40_anchor1, "KO_40hrs_Anchor1")
KO_40hrs_Anchor2_Ann <- annotateIRs(KO_40_anchor2, "KO_40hrs_Anchor2")

WT_40hrs_Anchor1_Ann <- annotateIRs(WT_40_anchor1, "WT_40hrs_Anchor1")
WT_40hrs_Anchor2_Ann <- annotateIRs(WT_40_anchor2, "WT_40hrs_Anchor2")

KO_16hrs_Anchor1_Ann <- annotateIRs(KO_16_anchor1, "KO_16hrs_Anchor1")
KO_16hrs_Anchor2_Ann <- annotateIRs(KO_16_anchor1, "KO_16hrs_Anchor2")

WT_16hrs_Anchor1_Ann <- annotateIRs(WT_16_anchor1, "WT_16hrs_Anchor1")
WT_16hrs_Anchor2_Ann <- annotateIRs(WT_16_anchor2, "WT_16hrs_Anchor2")


#get unique genes for each condition
get_genes2 <- function(GRangesDF){
  genes <- c(GRangesDF$Description[grepl("gene_id=PF3D7_", GRangesDF$Description)])
  qvals <- c(GRangesDF$p.val[grepl("gene_id=PF3D7_", GRangesDF$Description)])
  genes <- str_split(genes,";")
  genes <- lapply(genes, '[[',3)
  genes <- unlist(genes, recursive=FALSE)
  genes <- substring(genes, 9, 21)
  idx <- !duplicated(genes)
  genes_qvals <- qvals[idx]
  names(genes_qvals) <- genes[idx]
  return(genes_qvals)
}

#####
#ORA

## Plotting function
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

## background genes. If you don't have this use all the genes
ia5.bk <- PlasmoDB_Annot_Prot$attribute[grepl("ID=PF3D7_", PlasmoDB_Annot_Prot$attribute)]  
ia5.bk <- str_split(ia5.bk,";")
ia5.bk <- lapply(ia5.bk, '[[',1)
ia5.bk <- unlist(ia5.bk, recursive=FALSE)
ia5.bk <- unique(ia5.bk)
ia5.bk <- substring(ia5.bk, 4, 16)
ia5.bk <- unique(ia5.bk)

## Get the GO terms from Biomart
all <- biomaRt::listDatasets(biomaRt::useMart( biomart="protists_mart", host="https://protists.ensembl.org"))
db <- useMart(biomart="protists_mart", host="https://protists.ensembl.org",dataset = "pfalciparum_eg_gene")

## Getch GO IDs
go_ids= getBM(attributes=c('go_id', 'ensembl_gene_id', 'namespace_1003'), filters='ensembl_gene_id', values=ia5.bk, mart=db)
temp <- listAttributes(db)
gene_2_GO=unstack(go_ids[,c(1,2)])

## remove genes without GO term
go_ids <- go_ids[!go_ids$go_id=="",]
gene_2_GO=unstack(go_ids[,c(1,2)])

selection <- function(allScore){ return(allScore < 0.05)}

# make named factor showing which genes are of interest
KO_40_genes <- get_genes(KO_40hrs_Anchor1_Ann, KO_40hrs_Anchor2_Ann)
KO_16_genes <- get_genes(KO_16hrs_Anchor1_Ann, KO_16hrs_Anchor2_Ann)
WT_40_genes <- get_genes(WT_40hrs_Anchor1_Ann, WT_40hrs_Anchor2_Ann)
WT_16_genes <- get_genes(WT_16hrs_Anchor1_Ann, WT_16hrs_Anchor2_Ann)

## Use only significantly perturbed genes
ORA <- function(geneList, nameResults){
GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO, geneSel=selection)

results <- runTest(GOdata, algorithm="weight01", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results, orderBy="KS", topNodes=length(usedGO(GOdata)))
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
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

pdf(paste("ORA/",nameResults,".pdf",sep = ""), height = 10)
print({
  dot(goEnrichment[1:dim(goEnrichment)[1],])
})
dev.off()

a <- goEnrichment

writexl::write_xlsx(list("Up"=a ),paste("ORA/GO_Enrichment",nameResults,".xlsx",sep = ""))
}

ORA(geneList = KO_40_genes, nameResults = "KO_40hrs")
ORA(geneList = WT_40_genes, nameResults = "WT_40hrs")
ORA(geneList = KO_16_genes, nameResults = "KO_16hrs")
ORA(geneList = WT_16_genes, nameResults = "WT_16hrs")

#####
#start with AP2 gene only - chr 11, 328877 - 336810 - find all the interactions that have this gene, annotate genes from regions interacting with AP2 and ORA
#check if var genes are within our significant interactions
#See what happens near PfAP2-P gene
#Genomic location = Chr11 (PF3D7_11_V3) 328877 - 336810
#the binning consists of 5000 bases windows
#the fragment that contains the gene is PF3D7_11_V3, 320 000-330 000, Midpoint is 325 000

annotateIRs2<-function(queryGRanges, outFileName){
  olaps_prot <- findOverlaps(queryGRanges, PlasmoDB_Annot_Prot, type = "any", select = "all") #only protein coding genes
  
  #generate GRanges with annotation for each queryHit
  genes_olap <- olaps_prot %>% as.data.frame %>% group_by(queryHits) %>%
    mutate(genes = PlasmoDB_Annot@elementMetadata$attribute[subjectHits]) %>%
    dplyr::select(queryHits, genes)
  
  GRanges_annotatedRegions <- queryGRanges[queryHits(olaps_prot)] 
  mcols(GRanges_annotatedRegions)$Description <- genes_olap$genes
  
  #save to file
  write.table(GRanges_annotatedRegions, file=paste("PfAP2_P_ORA/",outFileName, ".GRanges", sep = ""))
  
  return(GRanges_annotatedRegions)
}

#find overlaping IRs and PfAP2-P gene by midFragment==325000 and chr==11
PfAP2.P_WT_40hrs <- WT_40hrs %>% filter(chr1 == "PF3D7_11_V3" | chr2 == "PF3D7_11_V3") %>% 
  filter(fragmentMid1 == "325000" | fragmentMid2 == "325000")
PfAP2.P_KO_40hrs <- KO_40hrs %>% filter(chr1 == "PF3D7_11_V3" | chr2 == "PF3D7_11_V3") %>% 
  filter(fragmentMid1 == "325000" | fragmentMid2 == "325000")
PfAP2.P_WT_16hrs <- WT_16hrs %>% filter(chr1 == "PF3D7_11_V3" | chr2 == "PF3D7_11_V3") %>% 
  filter(fragmentMid1 == "325000" | fragmentMid2 == "325000")
PfAP2.P_KO_16hrs <- KO_16hrs %>% filter(chr1 == "PF3D7_11_V3" | chr2 == "PF3D7_11_V3") %>% 
  filter(fragmentMid1 == "325000" | fragmentMid2 == "325000")

#get IRs midFragments
locus_40WT <- c(PfAP2.P_WT_40hrs$fragmentMid1[PfAP2.P_WT_40hrs$fragmentMid1 != "325000"], PfAP2.P_WT_40hrs$fragmentMid2[PfAP2.P_WT_40hrs$fragmentMid2 != "325000"])
locus_40KO <- c(PfAP2.P_KO_40hrs$fragmentMid1[PfAP2.P_KO_40hrs$fragmentMid1 != "325000"], PfAP2.P_KO_40hrs$fragmentMid2[PfAP2.P_KO_40hrs$fragmentMid2 != "325000"])
locus_16WT <- c(PfAP2.P_WT_16hrs$fragmentMid1[PfAP2.P_WT_16hrs$fragmentMid1 != "325000"], PfAP2.P_WT_16hrs$fragmentMid2[PfAP2.P_WT_16hrs$fragmentMid2 != "325000"])
locus_16KO <- c(PfAP2.P_KO_16hrs$fragmentMid1[PfAP2.P_KO_16hrs$fragmentMid1 != "325000"], PfAP2.P_KO_16hrs$fragmentMid2[PfAP2.P_KO_16hrs$fragmentMid2 != "325000"])

#make GRanges for IRs
locus_40WT_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = rep("PF3D7_11_V3", length(locus_40WT)),
                                                           "start" = locus_40WT - 4999,
                                                           "end"= locus_40WT + 5000,
                                                     "p.val"=PfAP2.P_WT_40hrs$p_value.x), keep.extra.columns = TRUE)

locus_40KO_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = rep("PF3D7_11_V3", length(locus_40KO)),
                                                     "start" = locus_40KO - 4999,
                                                     "end"= locus_40KO + 5000,
                                                     "p.val"=PfAP2.P_KO_40hrs$p_value.x), keep.extra.columns = TRUE)

locus_16WT_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = rep("PF3D7_11_V3", length(locus_16WT)),
                                                     "start" = locus_16WT - 4999,
                                                     "end"= locus_16WT + 5000,
                                                     "p.val"=PfAP2.P_WT_16hrs$p_value.x), keep.extra.columns = TRUE)

locus_16KO_GR <- makeGRangesFromDataFrame(data.frame("seqnames" = rep("PF3D7_11_V3", length(locus_16KO)),
                                                     "start" = locus_16KO - 4999,
                                                     "end"= locus_16KO + 5000,
                                                     "p.val"=PfAP2.P_KO_16hrs$p_value.x), keep.extra.columns = TRUE)
#annotate IRs
locus_40WT_GR_Ann <- annotateIRs2(queryGRanges = locus_40WT_GR, outFileName = "locus_40WT_GR_Ann")
locus_40KO_GR_Ann <- annotateIRs2(queryGRanges = locus_40KO_GR, outFileName = "locus_40KO_GR_Ann")
locus_16WT_GR_Ann <- annotateIRs2(queryGRanges = locus_16WT_GR, outFileName = "locus_16WT_GR_Ann")
locus_16KO_GR_Ann <- annotateIRs2(queryGRanges = locus_16KO_GR, outFileName = "locus_16KO_GR_Ann")

#get genes
geneList_40WT <- get_genes2(GRangesDF = locus_40WT_GR_Ann)
geneList_40KO <- get_genes2(GRangesDF = locus_40KO_GR_Ann)
geneList_16WT <- get_genes2(GRangesDF = locus_16WT_GR_Ann)
geneList_16KO <- get_genes2(GRangesDF = locus_16KO_GR_Ann)
