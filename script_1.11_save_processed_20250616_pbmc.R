# ~/R/scCHANTS/script_1.11_save_processed_20250616_pbmc.R

# Load packages

#! check correct file is loaded, saved, that quality control metrics and dims have been changed, and PBMC vs T cell
library("R.utils")    # need for decompressing files
library("dplyr")
library("Seurat")
library("ggplot2")
library("clustree")
library("patchwork")


# Create file path for saving plots
file_path <- "/cephfs/volumes/hpc_data_prj/id_hill_sims_wellcda/c1947608-5b3a-4d60-8179-b8e0779d7319/scratch_tmp/scCHANTS/20250616_benchmark/plots"


# Load data
# load the unzipped files
# gene.column = 2 gives gene symbol/name instead of ENSEMBL ID
scCHANTS_data <- Read10X(data.dir = "/cephfs/volumes/hpc_data_prj/id_hill_sims_wellcda/c1947608-5b3a-4d60-8179-b8e0779d7319/scratch_tmp/scCHANTS/20250616_benchmark/20250616_benchmark/", gene.column=2)

# create seurat object
# could add arguments for min.cells or min.features to be included
scCHANTS <- CreateSeuratObject(counts = scCHANTS_data, project = "scCHANTS")
Assays(scCHANTS)


## Rename HTO and ADT ids
# add HTO and ADT assay to the seurat object as it is by default ignored by CreateSeuratObject
# extract feature names
HTO_names <- rownames(scCHANTS_data$`Antibody Capture`) [1:14]
ADT_names <- rownames(scCHANTS_data$`Antibody Capture`) [15:18]

# extract count matrices
HTO_counts <- scCHANTS_data$`Antibody Capture` [HTO_names, , drop=FALSE]
ADT_counts <- scCHANTS_data$`Antibody Capture` [ADT_names, , drop=FALSE]

# check names and order
rownames(HTO_counts)
rownames(ADT_counts)

# create mapping from long to short names. using - not _ as cannot have underscores in seurat object names
name_to_id <- c(
  "TotalSeq™-C0251 anti-human Hashtag 1 Antibody"  = "HTO-1",
  "TotalSeq™-C0252 anti-human Hashtag 2 Antibody"  = "HTO-2",
  "TotalSeq™-C0253 anti-human Hashtag 3 Antibody"  = "HTO-3",
  "TotalSeq™-C0254 anti-human Hashtag 4 Antibody"  = "HTO-4",
  "TotalSeq™-C0255 anti-human Hashtag 5 Antibody"  = "HTO-5",
  "TotalSeq™-C0256 anti-human Hashtag 6 Antibody"  = "HTO-6",
  "TotalSeq™-C0257 anti-human Hashtag 7 Antibody"  = "HTO-7",
  "TotalSeq™-C0258 anti-human Hashtag 8 Antibody"  = "HTO-8",
  "TotalSeq™-C0259 anti-human Hashtag 9 Antibody"  = "HTO-9",
  "TotalSeq™-C0260 anti-human Hashtag 10 Antibody" = "HTO-10",
  "TotalSeq™-C0296 anti-human Hashtag 11 Antibody" = "HTO-11",
  "TotalSeq™-C0262 anti-human Hashtag 12 Antibody" = "HTO-12",
  "TotalSeq™-C0263 anti-human Hashtag 13 Antibody" = "HTO-13",
  "TotalSeq™-C0264 anti-human Hashtag 14 Antibody" = "HTO-14",
  
  # ADTs
  "TotalSeq-C anti human CD14 (clone M5E2, cat 301859)"             = "ADT-CD14",
  "TotalSeq-C anti human CD16 (clone 3G8 cat 302065)"               = "ADT-CD16",
  "TotalSeq™-C0045 anti-human CD4 Antibody (clone SK3, cat 344651)" = "ADT-CD4",
  "TotalSeq™-C0046 anti-human CD8 Antibody (clone SK1, cat 344753)" = "ADT-CD8"
)

# rename the rows (features) using the mapping
rownames(HTO_counts) <- unname(name_to_id[rownames(HTO_counts)])
rownames(ADT_counts) <- unname(name_to_id[rownames(ADT_counts)])

# check names have changed
rownames(HTO_counts)
rownames(ADT_counts)


## Add HTO and ADT assays to seurat object
# add HTO assay
scCHANTS[["HTO"]] <- CreateAssayObject(counts = HTO_counts)

# add ADT assay
scCHANTS[["ADT"]] <- CreateAssayObject(counts = ADT_counts)

# validate assays were added
Assays(scCHANTS)

# validate layers
Layers(scCHANTS)

# check default assay
DefaultAssay(scCHANTS)

# if needed, can change default assay
DefaultAssay(scCHANTS) <- "RNA"

# check number of cells (samples) and features looks right
scCHANTS@assays$RNA
scCHANTS@assays$HTO
scCHANTS@assays$ADT

# free space
rm(scCHANTS_data)


# Demultiplex HTOs
DefaultAssay(scCHANTS) <- "HTO"
    
# normalise HTO counts, standard is centered log-ratio (CLR) transformation
# can normalise across features (margin = 1) or across cells (margin = 2)
# comment on github form satija lab that recommend normalising across tags if variation between hash performances but since have very low signals, will try across cells
scCHANTS <- NormalizeData(scCHANTS, assay = "HTO", normalization.method = "CLR", margin = 2)
    
# using default threshold (quantile of inferred negative distribution for each hashtag), default kfunc is clara (for fast k-medoids clustering on large applications, it is faster and more memory-efficient as it repeatedly samples subsets and clusters those)
scCHANTS <- HTODemux(scCHANTS, assay = "HTO", positive.quantile = 0.99, kfunc = "clara")

# Reassign scCHANTS to only singlets
length(Cells(scCHANTS))
    
scCHANTS <- subset(scCHANTS, subset = HTO_classification.global == "Singlet")
    
# check cell numbers
length(Cells(scCHANTS))
length(Cells(scCHANTS[["RNA"]]))
length(Cells(scCHANTS[["HTO"]]))
length(Cells(scCHANTS[["ADT"]]))
## not subsetting RNA assay

    
### Remove demultiplexing plots to save memory
# remove objects
rm(ADT_counts, HTO_counts, ADT_names, HTO_names, hto_info, name_to_id)

    
# Add metadata
# examine existing metadata
metadata <- scCHANTS@meta.data
    
dim(metadata)
head(metadata)
summary(metadata$nCount_RNA)
    
# read in additional sample metadata
metadata_samples <- read.csv("/cephfs/volumes/hpc_data_prj/id_hill_sims_wellcda/c1947608-5b3a-4d60-8179-b8e0779d7319/scratch_tmp/scCHANTS/20250616_benchmark/sample_metadata.csv")
    
head(metadata_samples)
    
# steps to add to seurat object
    
# align each cell's HTO id to row in the metadata_samples, rownames here represent HTOs but can't have duplicate rownames hence 4.1 etc
hto_info <- metadata_samples[match(scCHANTS$HTO_classification, metadata_samples$HTO), ]
head(hto_info)
    
# add each row to new variable in metadata
scCHANTS$sample_id <- hto_info$sample
scCHANTS$timepoint <- hto_info$timepoint
scCHANTS$treatment <- hto_info$treatment
scCHANTS$PBMC_or_T <- hto_info$PBMC_or_T
    
# check new metadata
metadata <- scCHANTS@meta.data

dim(metadata)
colnames(metadata)
metadata[14:17, 14:17]

    
# Separate PBMC and sorted T cell samples
print("Full scCHANTS")
DefaultAssay(scCHANTS) <- "RNA"
length(Cells(scCHANTS@assays$RNA))
length(Cells(scCHANTS@assays$HTO))
length(Cells(scCHANTS@assays$ADT))
    
print("PBMC filtered scCHANTS")
scCHANTS_pbmc <- subset(scCHANTS, subset = PBMC_or_T == "PBMC")
length(Cells(scCHANTS_pbmc@assays$RNA))
length(Cells(scCHANTS_pbmc@assays$HTO))
length(Cells(scCHANTS_pbmc@assays$ADT))
    
print("T cells filtered scCHANTS")
scCHANTS_t <- subset(scCHANTS, subset = PBMC_or_T == "T")
length(Cells(scCHANTS_t@assays$RNA))
length(Cells(scCHANTS_t@assays$HTO))
length(Cells(scCHANTS_t@assays$ADT))
 
    
# PBMC pre-processing
    
## Quality control on PBMCs
# calculate mt contamination and add column to metadata
scCHANTS_pbmc[["percent_mt"]] <- PercentageFeatureSet(scCHANTS_pbmc, pattern = "^MT-")
    
# calculate mt contamination and add column to metadata
scCHANTS_pbmc[["percent_rb"]] <- PercentageFeatureSet(scCHANTS_pbmc, pattern = "^RP[SL]")
    
    
## Do subsetting
# can subset straightaway based on QC metrics e.g. :
#scCHANTS_pbmc <- subset(scCHANTS_pbmc, subset = nFeature_RNA > 130  & nFeature_RNA < 800 & percent_mt < 6) 
    
# or assign new metadata column
# start with all cells labeled as "initial"
scCHANTS_pbmc$QC <- "initial"
    
# label low nFeature_RNA (< 130)
## i.e. if nFeature_RNA < 130 (and previously passed QC as cautionary step), call it low_nFeature, otherwise leave as initial
    
scCHANTS_pbmc$QC <- ifelse(
  scCHANTS_pbmc$nFeature_RNA < 130 & scCHANTS_pbmc$QC == "initial",
  "low_nFeature",
  scCHANTS_pbmc$QC
  )
    
# label high nFeature_RNA (> 800)
scCHANTS_pbmc$QC <- ifelse(
  scCHANTS_pbmc$nFeature_RNA > 800 & scCHANTS_pbmc$QC == "initial",
  "high_nFeature",
  scCHANTS_pbmc$QC
  )
    
    # label high mitochondrial content (> 6%)
scCHANTS_pbmc$QC <- ifelse(
  scCHANTS_pbmc$percent_mt > 6 & scCHANTS_pbmc$QC == "initial",
  "high_mt",
  scCHANTS_pbmc$QC
  )

# label combined: low_nFeature + high_mt
scCHANTS_pbmc$QC <- ifelse(
  scCHANTS_pbmc$nFeature_RNA < 130 & scCHANTS_pbmc$percent_mt > 6 & scCHANTS_pbmc$QC != "pass",
  "low_nFeature_high_mt",
  scCHANTS_pbmc$QC
  )

    # label combined: high_nFeature + high_mt
scCHANTS_pbmc$QC <- ifelse(
  scCHANTS_pbmc$nFeature_RNA > 800 & scCHANTS_pbmc$percent_mt > 6 & scCHANTS_pbmc$QC != "pass",
  "high_nFeature_high_mt",
  scCHANTS_pbmc$QC
)

# finally, label cells that pass all filters
scCHANTS_pbmc$QC <- ifelse(
  scCHANTS_pbmc$nFeature_RNA >= 130 &
scCHANTS_pbmc$nFeature_RNA <= 800 &
scCHANTS_pbmc$percent_mt <= 6 &
scCHANTS_pbmc$QC == "initial",
  "pass",
  scCHANTS_pbmc$QC
)


# subset
scCHANTS_pbmc <- subset(scCHANTS_pbmc, subset = QC == "pass")
length(Cells(scCHANTS_pbmc@assays$RNA))
length(Cells(scCHANTS_pbmc@assays$HTO))
length(Cells(scCHANTS_pbmc@assays$ADT))

  

## Normalisation of PBMCs
DefaultAssay(scCHANTS_pbmc) <- "RNA"
DefaultLayer(scCHANTS_pbmc[["RNA"]])

# normalise
scCHANTS_pbmc <- NormalizeData(scCHANTS_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# find variable features
scCHANTS_pbmc <- FindVariableFeatures(scCHANTS_pbmc, selection.method = "vst", nfeatures = 2000, layer = "counts.Gene Expression")

# identify 10 most highly variable genes
top10_features <- head(VariableFeatures(scCHANTS_pbmc), 10)

## NAs introduced are all for HTO or ADT rows
hvf_info <- HVFInfo(scCHANTS_pbmc)
summary(is.na(hvf_info))
hvf_info[!complete.cases(hvf_info), ]

# scale data
all.genes <- rownames(scCHANTS_pbmc)
scCHANTS_pbmc <- ScaleData(scCHANTS_pbmc, features = all.genes)


## Dim red on PBMC
# run PCA
scCHANTS_pbmc <- RunPCA(scCHANTS_pbmc, features = VariableFeatures(object = scCHANTS_pbmc))

 
## Deciding resolution for clustering
scCHANTS_pbmc <- FindNeighbors(scCHANTS_pbmc, dims = 1:11) 

scCHANTS_pbmc <- FindClusters(scCHANTS_pbmc, resolution = 0.6)
print("Finding clusters")
print("Metadata columns are:")
colnames(scCHANTS_pbmc@meta.data)

# set resolution as default
Idents(scCHANTS_pbmc) <- "RNA_snn_res.0.6"


## Run umap and visualise
scCHANTS_pbmc <- RunUMAP(scCHANTS_pbmc, dims = 1:11)


# Save PBMC pre-processed RDS (takes a long time), save before joining layers to find markers
saveRDS(object=scCHANTS_pbmc, file ="/cephfs/volumes/hpc_data_prj/id_hill_sims_wellcda/c1947608-5b3a-4d60-8179-b8e0779d7319/scratch_tmp/scCHANTS/20250616_benchmark/20250616_scCHANTS_pbmc_processed.Rds")

print("RDS saved :)")


# Find markers for clusters
# must join layers to do find markers
scCHANTS_pbmc <- JoinLayers(scCHANTS_pbmc)
    
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc_markers <- FindAllMarkers(scCHANTS_pbmc, only.pos = TRUE)
    
write.csv(pbmc_markers,"/cephfs/volumes/hpc_data_prj/id_hill_sims_wellcda/c1947608-5b3a-4d60-8179-b8e0779d7319/scratch_tmp/scCHANTS/20250616_benchmark/20250616_scCHANTS_pbmc_markers.csv")

print("CSV saved :)")

pbmc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
    
    
    