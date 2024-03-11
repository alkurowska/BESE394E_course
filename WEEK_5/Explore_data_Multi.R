# Load necessary libraries for single-cell analysis, motif analysis, and plotting
library(Seurat)
library(TFBSTools)
library(motifmatchr)
library(Signac)
library(ggplot2)
library(Azimuth)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)

# Multiome https://www.10xgenomics.com/datasets/mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-saltyez-protocol-and-10x-complex-tissue-dp-ct-sorted-and-ct-unsorted-1-standard
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.tar.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz.tbi

# Create a function to plot quality control features of the Seurat object

plot_QC_features <- function(seurat_object, pdf_path){
  pdf(pdf_path, width =10) # Start a PDF output to save plots
  # shows the distribution of the transcripts per cells
  # Generate violin plots for a range of QC metrics including counts of RNA and ATAC, features, mitochondrial DNA percentage, etc.
  print(VlnPlot(
    object = seurat_object,
    features = c("nCount_RNA", "nCount_ATAC", 'nFeature_RNA',  'mitoPct', 
                 'FRiP', 'log10GenesPerUMI', "TSS.enrichment", 
                 "nucleosome_signal"),
    pt.size = 0, ncol =4))
  DefaultAssay(seurat_object) <- "ATAC" # Switch default assay to ATAC and plot fragment histogram and TSS enrichment plots
  print(FragmentHistogram(object = seurat_object, region = 'chr1-1-10000000', group.by = 'nucleosome_group'))  # Plot fragment histogram for ATAC data showing nucleosome banding pattern
  print(TSSPlot(seurat_object, group.by = 'high.tss') + NoLegend()) # TSS enrichment plot to assess ATAC signal around transcription start sites
  print(DensityScatter(seurat_object, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)) # Scatter plot showing correlation between ATAC count and TSS enrichment
  print(ggplot(seurat_object@meta.data, aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
          geom_density(alpha = 0.2) + 
          theme_classic() +
          scale_x_log10() + 
          geom_vline(xintercept = 300) + ggtitle('GENESpercell'))
  # Visualize the distribution of genes detected per cell via boxplot
  print(ggplot(seurat_object@meta.data, aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
          geom_boxplot() + 
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
          theme(plot.title = element_text(hjust=0.5, face="bold")) +
          ggtitle("NCells vs NGenes"))
  # correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  print(ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, colour=mitoPct, group=orig.ident)) + 
          geom_point() + 
          scale_colour_gradient(low = "gray90", high = "black", limits=c(0,100)) +
          stat_smooth(method=lm) +
          scale_x_log10() + 
          scale_y_log10() + 
          theme_classic() +
          geom_vline(xintercept = 500) +
          geom_hline(yintercept = 6000) +
          geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(seurat_object))))
  dev.off() # Close the PDF file
  
}

# Load 10x Genomics data using the Read10X function and specify the paths to the the ATAC-seq fragments
counts <- Read10X("./filtered_feature_bc_matrix")
fragpath <- './M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz'

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation)) # Set chromosome names

# create a Seurat object containing the RNA adata
mouse_brain <- CreateSeuratObject(
  counts = counts[['Gene Expression']],
  assay = "RNA"
)

# create ATAC assay and add it to the object
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes

# Preparing ATAC-seq data for analysis
atac_counts <- counts$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts) # Filter with standards chromosomes
atac_counts <- atac_counts[as.vector(grange.use), ]

# Add ATAC-seq data to the Seurat object with chromatin assay
mouse_brain[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

################################################################################
# Calculate metrics for QC
################################################################################

DefaultAssay(mouse_brain) <- "ATAC"
# Calculating nucleosome signal, TSS enrichment, and fragment counts for ATAC

mouse_brain <- NucleosomeSignal(mouse_brain)
mouse_brain <- TSSEnrichment(mouse_brain, fast=FALSE)
total_fragments <- CountFragments(fragments = fragpath)
rownames(total_fragments) <- total_fragments$CB
mouse_brain$fragments <- total_fragments[colnames(mouse_brain), "frequency_count"]
mouse_brain <- FRiP(object = mouse_brain, assay = 'ATAC', total.fragments = 'fragments')
mouse_brain$nucleosome_group <- ifelse(mouse_brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
mouse_brain$high.tss <- ifelse(mouse_brain$TSS.enrichment > 3, 'High', 'Low')

DefaultAssay(mouse_brain) <- "RNA"
# Setting up various QC metrics for RNA data including mitochondrial gene percentage, ribosomal protein gene percentage, and ratio of features to UMIs

mouse_brain$mitoPct <- PercentageFeatureSet(mouse_brain, pattern = "^mt-")
mouse_brain$RPSPct  <- PercentageFeatureSet(object = mouse_brain, pattern = "^Rp[sl]")
mouse_brain$log10GenesPerUMI <- log10(mouse_brain$nFeature_RNA) / log10(mouse_brain$nCount_RNA)

# Plot quality control features before applying filters using a previously defined function
plot_QC_features(mouse_brain, '/home/kurowsaa/5k_mouse_brain_GEX_QC_Pre.pdf')

# Save the pre-QC Seurat object 
saveRDS(mouse_brain, './mouse_brain_multiome_Pre_QC.rds')

################################################################################
# Filter cells based on QC metrics
################################################################################

mouse_brain <- subset(
  x = mouse_brain,
  subset = nCount_ATAC < 100000 & # Filter cells with ATAC read counts less than 100,000
    nCount_RNA < 25000 & # Filter cells with RNA read counts less than 25,000
    nCount_ATAC > 1000 & # Exclude cells with ATAC read counts below 1,000
    nCount_RNA > 1000 & # Exclude cells with RNA read counts below 1,000
    nucleosome_signal < 2 & # Exclude cells with high nucleosome signal (indicative of poor ATAC quality)
    TSS.enrichment > 1 & # Select cells with TSS enrichment above 1 (higher is better)
    FRiP > 0.2 # Select cells with Fragment Reads in Peaks (FRiP) score above 0.2
)

# Save the filtered Seurat object for post-QC analysis
saveRDS(mouse_brain, './mouse_brain_multiome_Post_QC.rds')
# mouse_brain <- readRDS('/ibex/user/serrang/Projects_data/MultiomeCourse/Data/mouse_brain_multiome_Post_QC.rds')

# Plot quality control features after applying filters using a previously defined function
plot_QC_features(mouse_brain, './5k_mouse_brain_GEX_QC_Post.pdf')



################################################################################
# Process the data (normalization, scaling, PCA, clustering, UMAP, ....)
################################################################################


setwd("/home/kurowsaa/")
# Load the Seurat object post-QC for further analysis
mouse_brain <- readRDS("mouse_brain_multiome_Post_QC.rds")

# Normalize and identify variable features within the RNA assay
DefaultAssay(mouse_brain) <- "RNA"
mouse_brain <- NormalizeData(mouse_brain)
mouse_brain <- FindVariableFeatures(mouse_brain, 
                                    selection.method = "vst", 
                                    nfeatures = 2000)
mouse_brain <- ScaleData(mouse_brain)
mouse_brain <- RunPCA(mouse_brain) # Principal Component Analysis for dimensionality reduction

# Perform clustering based on PCA results
mouse_brain <- FindNeighbors(mouse_brain, dims = 1:30)
mouse_brain <- FindClusters(mouse_brain, 
                            resolution = 0.4, 
                            algorithm = 3, 
                            cluster.name="RNA_clusters_03")

mouse_brain <- RunPCA(mouse_brain)
# Run UMAP for visualization of the RNA-based clustering
mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
                       reduction = "pca", 
                       reduction.name = "rna_umap")

# Use Azimuth to map single-cell data to reference datasets for annotation
mouse_brain <- Azimuth::RunAzimuth(mouse_brain, 
                                   reference = "mousecortexref")


# Plot UMAP visualizations based on RNA data and predicted cell subclasses
pdf('./UMAP_RNA.pdf')
DimPlot(mouse_brain, reduction = "rna_umap", group.by='RNA_clusters_03', 
        label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "rna_umap", group.by = "predicted.subclass", 
        label = TRUE, label.size = 3) + NoLegend()
dev.off()

# Switch to the ATAC assay and perform analysis specific to chromatin accessibility data
DefaultAssay(mouse_brain) <- "ATAC"
mouse_brain <- FindTopFeatures(mouse_brain, min.cutoff = 5) # Identify top accessible regions
mouse_brain <- RunTFIDF(mouse_brain) # Perform Term Frequency-Inverse Document Frequency (TF-IDF) normalization 
mouse_brain <- RunSVD(mouse_brain) # Singular Value Decomposition for dimensionality reduction on ATAC data

# Plot depth correlation as a quality metric
pdf('./Correlation.pdf')
DepthCor(mouse_brain)
dev.off()

# Run UMAP for visualization based on ATAC data
mouse_brain <- RunUMAP(mouse_brain, dims=1:30, 
                       reduction = "lsi", 
                       reduction.name = "atac_umap")

# Plot UMAP visualizations based on ATAC data and predicted cell subclasses
pdf('./UMAP_ATAC.pdf')
DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", group.by = "predicted.subclass",
        label = TRUE, label.size = 3) + NoLegend()
dev.off()

# Save the Seurat object with RNA and ATAC data
saveRDS(mouse_brain, './mouse_brain_multiome.rds')

################################################################################
# Integrate the RNA and ATAC data
################################################################################

# Integrate RNA and ATAC assays using multimodal neighbors to leverage both datasets
mouse_brain <- FindMultiModalNeighbors(
  object = mouse_brain,
  reduction.list = list("pca", "lsi"), # PCA for RNA, LSI for ATAC
  dims.list = list(1:50, 1:40),  # Dimensions to use from each reduction
  verbose = TRUE
)

# build a joint UMAP visualization
mouse_brain <- RunUMAP(
  object = mouse_brain, 
  reduction.name = "wnn.umap", dims = 1:30,
  assay = "RNA",
  verbose = TRUE
)

# Save the integrated Seurat object
saveRDS(mouse_brain, './mouse_brain_multiome_Integrated.rds')

# Plot UMAP visualizations for integrated data and predicted subclasses
pdf('./Umap_integrated.pdf')
DimPlot(mouse_brain, reduction = "rna_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "atac_umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
DimPlot(mouse_brain, reduction = "wnn.umap", label = TRUE, group.by = "predicted.subclass",
        repel = TRUE) + NoLegend()
dev.off()

# Plot violin plots for RNA and ATAC weights ("contribution" to the integrated space)
pdf('./Umap_integrated_weights.pdf')
Idents(mouse_brain) <- "predicted.subclass"
VlnPlot(mouse_brain, features = c("RNA.weight", "ATAC.weight"),  pt.size = 0, ncol = 1)
dev.off()

################################################################################
# Perform a Differential accessibility analysis
################################################################################

mouse_brain <- readRDS('./mouse_brain_multiome_Integrated.rds')

# Set the default assay to ATAC for differential accessibility analysis
DefaultAssay(mouse_brain) <- 'ATAC'

# Define cell identities for differential analysis
Idents(mouse_brain) <- "predicted.subclass"

# Perform differential accessibility analysis between Oligo and Astro subclasses using a logistic regression model, correcting for total ATAC read counts
da_peaks <- FindMarkers(
  object = mouse_brain,
  ident.1 = c("Oligo"), 
  ident.2 = c("Astro"),
  test.use = 'LR', # Logistic Regression
  latent.vars = 'nCount_ATAC'
)
# Save the differential analysis results
saveRDS(da_peaks, './DA_peaks.rds')

# Sort differential peaks by log fold change in descending order
da_peaks <- da_peaks[order(da_peaks$avg_log2FC, decreasing = TRUE), ]

# Generate plots for differentially accessible peaks and associated genes
pdf('./Peaks_DA.pdf', width=12, height=20)
cowplot::plot_grid(
  # Violin plots for top differentially accessible peaks in ATAC data
  VlnPlot(
    object = mouse_brain,
    assay = 'ATAC',
    features = rownames(da_peaks)[1:3],
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  # Violin plots for gene expression levels of genes closest to top differentially accessible peaks
  VlnPlot(
    object = mouse_brain,
    assay = 'RNA',
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name,
    pt.size = 0.1,
    group.by='predicted.subclass', 
    ncol=3
  ),
  # Feature plots showing the spatial distribution of top differentially accessible peaks across cells in UMAP space
  FeaturePlot(
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = rownames(da_peaks)[1:3],
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
  # Feature plots showing the spatial distribution of gene expression for genes closest to top DA peaks across cells in UMAP space
  FeaturePlot(
    object = mouse_brain,
    reduction = 'wnn.umap', 
    order=TRUE,
    features = ClosestFeature(mouse_brain, rownames(da_peaks)[1:3])$gene_name,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol=3
  ) & NoLegend(),
  nrow=4)
dev.off()


# Annotate peaks with his closest feature
# Identify peaks with significantly different accessibility (log2 fold change > 3 for Oligo, < 3 for Astro) and get their closest gene
open_Oligo <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_Astro <- rownames(da_peaks[da_peaks$avg_log2FC < 3, ])
closest_Oligo <- ClosestFeature(mouse_brain, open_Oligo)
closest_Astro <- ClosestFeature(mouse_brain, open_Astro)

# https://www.cellsignal.com/pathways/neuronal-and-glial-cell-markers
# Visualize the coverage of the peaks
frags <- Fragments(mouse_brain)  # get list of fragment objects
Fragments(mouse_brain) <- NULL  # remove fragment information from assay
newpath <- "./M_Brain_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]], new.path = newpath)  # update path. Do this for any/all fragment objects in the list
Fragments(mouse_brain) <- frags  # assign update list of fragment objects back to the assay

pdf('./Coverage_selected.pdf', height=12)
CoveragePlot(
  object = mouse_brain,
  region = c('Olig1', # Select genes representing oligodendrocytes and astrocytes
             'Gfap'),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)
dev.off()


################################################################################
# Motif analysis
################################################################################

# Motif analysis with the DA peaks
pfm <- getMatrixSet(
  x = JASPAR2020, # Specify the JASPAR database
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE) # Filter options for core collection of vertebrate motifs
)


# add motif information sing the Position Frequency Matrices (PFM) obtained from JASPAR
mouse_brain <- AddMotifs(
  object = mouse_brain, # Seurat object to add motifs to
  genome = BSgenome.Mmusculus.UCSC.mm10, # Specify the genome version for motif analysis
  pfm = pfm # The matrix set from JASPAR used for motif analysis
)

# Identify top differentially accessible peaks based on a p-value threshold
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

# Set the default assay in Seurat object to ATAC for subsequent ATAC-specific analyses
DefaultAssay(mouse_brain) <- "ATAC"

# Obtain Position Weight Matrices (PWM) set for mouse from JASPAR2020
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 'Mus musculus', all_versions = FALSE))

# Create a motif matrix using the PWM set and the genomic ranges of peaks in the Seurat object
motif.matrix <- CreateMotifMatrix(features = granges(mouse_brain), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)

# Create a motif object from the motif matrix and the PWM set
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)

# Update the Seurat object with the motif data
mouse_brain <- SetAssayData(mouse_brain, assay = 'ATAC', layer = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
# Run ChromVAR to analyze variability in motif accessibility, which can indicate transcription factor activity
mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10 # Specify the genome version
)

# Find motifs enriched in the top differentially accessible peaks
enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top.da.peak # Peaks identified as significantly differentially accessible
)


mouse_brain <- readRDS('./mouse_brain_multiome_Motif.rds')

# Plot the enriched motifs to visualize motif enrichment
pdf('./Motif_enrichment.pdf')
MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

saveRDS(mouse_brain, './mouse_brain_multiome_Motif.rds')

################################################################################
# Generate a RNA activity matrix based on the ATAC-seq data
################################################################################

# We can create a proxy of the gene expression from the ATAC-seq data using  GeneActivity function. 
# We can also create a proxy of the gene expression from the ATAC-seq data using the chromVAR package. 
# This package uses the motif accessibility to infer the gene expression. 
# We can use the motif models from JASPAR2020 to perform this analysis.

# Create a proxy of gene expression from ATAC-seq data using the GeneActivity function
gene.activities <- GeneActivity(mouse_brain)
mouse_brain[['RNA_ACTIVITY']] <- CreateAssayObject(counts = gene.activities) # Add the gene activity data as a new assay in the Seurat object

# Normalize the gene activity data using log normalization
mouse_brain <- NormalizeData(
  object = mouse_brain,
  assay = 'RNA_ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(mouse_brain$nCount_RNA) # Use median RNA count as scaling factor
)

# Set default assay to 'RNA' for gene expression analysis
DefaultAssay(mouse_brain) <- 'RNA'

# Plot feature plots for selected genes based on RNA data
RNA_plot <- FeaturePlot(
  object = mouse_brain,
  order=TRUE,
  features =c("Sema5a","Dennd4a","Nkain1"), # Genes of interest
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)& NoLegend()

# Set default assay to 'RNA_ACTIVITY' to analyze gene activity based on ATAC-seq
DefaultAssay(mouse_brain) <- 'RNA_ACTIVITY'

# Plot feature plots for the same selected genes based on RNA activity inferred from ATAC-seq
RNA_Activity_plot <- FeaturePlot(
  object = mouse_brain,
  order=TRUE,
  features =c("Sema5a","Dennd4a","Nkain1"), # Genes of interest
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)& NoLegend()

# Save comparison plots between RNA and RNA activity to PDF
pdf('./RNA_comparison.pdf', height=12)
cowplot::plot_grid(
  RNA_plot,
  RNA_Activity_plot,
  nrow=2)
dev.off()

saveRDS(mouse_brain, './mouse_brain_multiome_RNA_Activity.rds')


# ASSIGNMENT 3 
# Generate motifs for the top 5 differentially accessible peaks per cell type. 
# Use the JASPAR database to obtain the motifs.


# Analyze all of the 4 cell lines in a loop
cell_lines <- c('Pvalb', 'Lamp5','Endo', 'OPC')

for (i in 1:length(cell_lines)) {
  # Generate motifs for the top 5 differentially accessible peaks per cell type
  mouse_brain <- readRDS('./mouse_brain_multiome_Integrated.rds')
  
  # Set the default assay to ATAC for differential accessibility analysis
  DefaultAssay(mouse_brain) <- 'ATAC'
  
  # Define cell identities for differential analysis
  Idents(mouse_brain) <- "predicted.subclass"
  
  da_peaks <- FindMarkers(
    object = mouse_brain,
    ident.1 = cell_lines[i], 
    ident.2 = NULL, # Against all other cells
    #test.use = 'LR', # using default as LR took 6h ! :o 
    # However, based on existing knowledge and common practices in the field, the Wilcoxon test 
    # is often used for differential analysis in scATAC-seq data due to its non-parametric nature.
    # making it suitable for the sparse and binary-like distribution of accessibility data. 
    #latent.vars = 'nCount_ATAC'
  )
  
  # Save the differential analysis results for each cell type
  saveRDS(da_peaks, paste0('./DA_peaks_',cell_lines[i],'.rds'))
  
  # Sort differential peaks by log fold change in descending order
  da_peaks <- da_peaks[order(da_peaks$avg_log2FC, decreasing = TRUE), ]
  
  # Identify top differentially accessible peaks based on a p-value threshold
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
  
  da.peaks <- da_peaks[da_peaks$p_val < 0.005, ]
  da.peaks <- da.peaks[order(da.peaks$p_val, decreasing = FALSE), ]
  # Generate plots for differentially accessible peaks and associated genes
  pdf(paste0('./Peaks_DA_',cell_lines[i],'.pdf'), width=12, height=20)
  cowplot::plot_grid(
    # Violin plots for top differentially accessible peaks in ATAC data
    VlnPlot(
      object = mouse_brain,
      assay = 'ATAC',
      features = rownames(da.peaks)[1:3],
      pt.size = 0.1,
      group.by='predicted.subclass', 
      ncol=3
    ),
    # Violin plots for gene expression levels of genes closest to top differentially accessible peaks
    VlnPlot(
      object = mouse_brain,
      assay = 'RNA',
      features = ClosestFeature(mouse_brain, rownames(da.peaks)[1:3])$gene_name,
      pt.size = 0.1,
      group.by='predicted.subclass', 
      ncol=3
    ),
    # Feature plots showing the spatial distribution of top differentially accessible peaks across cells in UMAP space
    FeaturePlot(
      object = mouse_brain,
      reduction = 'wnn.umap', 
      order=TRUE,
      features = rownames(da.peaks)[1:3],
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol=3
    ) & NoLegend(),
    # Feature plots showing the spatial distribution of gene expression for genes closest to top DA peaks across cells in UMAP space
    FeaturePlot(
      object = mouse_brain,
      reduction = 'wnn.umap', 
      order=TRUE,
      features = ClosestFeature(mouse_brain, rownames(da.peaks)[1:3])$gene_name,
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol=3
    ) & NoLegend(),
    nrow=4)
  dev.off()
  
  
  # Top 5 differentially accessible peaks
  top.da.peak <- top.da.peak[1:5]
  
  # Motif analysis with the DA peaks
  pfm <- getMatrixSet(
    x = JASPAR2020, # Specify the JASPAR database
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE) # Filter options for core collection of vertebrate motifs
  )
  
  # add motif information sing the Position Frequency Matrices (PFM) obtained from JASPAR
  mouse_brain <- AddMotifs(
    object = mouse_brain, # Seurat object to add motifs to
    genome = BSgenome.Mmusculus.UCSC.mm10, # Specify the genome version for motif analysis
    pfm = pfm # The matrix set from JASPAR used for motif analysis
  )
  
  # Obtain Position Weight Matrices (PWM) set for mouse from JASPAR2020
  pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 'Mus musculus', all_versions = FALSE))
  
  # Create a motif matrix using the PWM set and the genomic ranges of peaks in the Seurat object
  motif.matrix <- CreateMotifMatrix(features = granges(mouse_brain), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
  
  # Create a motif object from the motif matrix and the PWM set
  motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
  
  # Update the Seurat object with the motif data
  mouse_brain <- SetAssayData(mouse_brain, assay = 'ATAC', layer = 'motifs', new.data = motif.object)
  
  # Note that this step can take 30-60 minutes 
  # Run ChromVAR to analyze variability in motif accessibility, which can indicate transcription factor activity
  mouse_brain <- RunChromVAR(
    object = mouse_brain,
    genome = BSgenome.Mmusculus.UCSC.mm10 # Specify the genome version
  )
  
  # Find motifs enriched in the top differentially accessible peaks
  enriched.motifs <- FindMotifs(
    object = mouse_brain,
    features = top.da.peak # Peaks identified as significantly differentially accessible
  )
  
  saveRDS(enriched.motifs, paste0('./Motif_enrichment_',cell_lines[i],'.rds'))
  
  # Plot the enriched motifs to visualize motif enrichment
  pdf(paste0('./Motif_enrichment_',cell_lines[i],'.pdf'))
  MotifPlot(
    object = mouse_brain,
    motifs = head(rownames(enriched.motifs))
  )
  dev.off()
}
