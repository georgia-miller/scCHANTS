
# Script to do all processing of 20250805_run1_pbmcs
  
# combines: script_1.Rmd : HTO demultiplexing, add metadata, split into PBMC and T cells, do QC on PBMCs
# and script_2a.Rmd which also : SCTransform normalisaiton, dimensional reduction (11 PCs), clustering (at resolution 0.7), saves the Rds, finds cluster biomarkers and saves them as a csv

library("dplyr")
library("tidyr")
library("Seurat")
library("sctransform")


# before doing sctransform for first time, allows faster running
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("glmGamPoi")

# install.packages("sctransform")


demux_qc <- function(date, run, donor, pbmc_or_t = c("pbmc", "t"), remote = "no", clr_norm = 2, demux_positive = 0.99, nFeature_min, nFeature_max, mt_max) {
  
  # check all args are valid
  stopifnot(
    remote %in% c("yes", "no"),
    pbmc_or_t %in% c("pbmc", "t"),
    clr_norm %in% c(1, 2),
    demux_positive >= 0 & demux_positive <= 0.99,
    is.numeric(nFeature_min), length(nFeature_min) == 1,
    is.numeric(nFeature_max), length(nFeature_max) == 1,
    is.numeric(mt_max), length(mt_max) == 1,
    # check donor has 3 digits
    grepl("^[0-9]{3}$", donor)
  )
  
  ## Load data
  # create local or remote file path
  date_run <- paste0(date, "_run", run)
  file_path <- if(remote == "no") {
    "/Users/k2477939/Library/CloudStorage/OneDrive-King\'sCollegeLondon/General\ -\ Laboratory\ of\ Molecular\ Immunotherapy\ and\ Antibiotics/Data/Georgia\ Miller/scCHANTS/"
  } else if(remote == "yes") {
    "/scratch/prj/id_hill_sims_wellcda/scCHANTS/"
  }
  data_path <- paste0(file_path, date_run, "/", date_run, "_", donor)
  
  # load 10x file
  scCHANTS_data <- Read10X(data.dir = data_path, gene.column=2)
  
  # create seurat object
  scCHANTS <- CreateSeuratObject(counts = scCHANTS_data, project = "scCHANTS")
  
  
  ## Rename HTO and ADT IDs
  # extract feature names
  HTO_names <- rownames(scCHANTS_data$`Antibody Capture`) [1:13]
  ADT_names <- rownames(scCHANTS_data$`Antibody Capture`) [14:17]
  
  # extract count matrices
  HTO_counts <- scCHANTS_data$`Antibody Capture` [HTO_names, , drop=FALSE]
  ADT_counts <- scCHANTS_data$`Antibody Capture` [ADT_names, , drop=FALSE]
  
  # to check names and order
  # message("New HTO names are:")
  # print(rownames(HTO_counts))
  # message("New ADT names are:")
  # print(rownames(ADT_counts))
  
  # create mapping from long to short names
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
  message("New HTO names are:")
  print(rownames(HTO_counts))
  message("New ADT names are:")
  print(rownames(ADT_counts))
  
  
  ## Add HTO and ADT assays to seurat object
  # add HTO assay
  scCHANTS[["HTO"]] <- CreateAssayObject(counts = HTO_counts)
  
  # add ADT assay
  scCHANTS[["ADT"]] <- CreateAssayObject(counts = ADT_counts)
  
  # validate assays and layers were added
  message("Assays are:")
  print(Assays(scCHANTS))
  
  # check number of cells (samples) and features looks right
  message("Features and cell numbers for RNA assay:")
  print(scCHANTS@assays$RNA)
  message("Features and cell numbers for HTO assay:")
  print(scCHANTS@assays$HTO)
  message("Features and cell numbers for ADT assay:")
  print(scCHANTS@assays$ADT)
  
  
  ## Demultiplex HTOs
  # remove empty cells, add total RNA and HTO counts to metadata then filter using them
  scCHANTS$totRNA <- Matrix::colSums(scCHANTS@assays$RNA@layers$`counts.Gene Expression`)
  scCHANTS$totHTO <- Matrix::colSums(scCHANTS@assays$HTO@counts)
  
  scCHANTS <- subset(scCHANTS, subset = totRNA > 0 | totHTO > 0)
  
  DefaultAssay(scCHANTS) <- "HTO"
  
  # normalise HTO counts, standard is centered log-ratio (CLR) transformation
  # can normalise across features (margin = 1) or across cells (margin = 2)
  # comment on github forum satija lab says they recommend normalising across tags if variation between hash performances but since have very low signals, will try across cells
  scCHANTS <- NormalizeData(scCHANTS, assay = "HTO", normalization.method = "CLR", margin = clr_norm)
  
  scCHANTS <- HTODemux(scCHANTS, assay = "HTO", positive.quantile = demux_positive, kfunc = "clara")
  
  DefaultAssay(scCHANTS) <- "HTO"
  
  message("Original number of cells:")
  print(length(Cells(scCHANTS)))
  
  scCHANTS <- subset(scCHANTS, subset = HTO_classification.global == "Singlet")
  
  message("Number of cells after subsetting for singlets in RNA, HTO and ADT assays (should all be the same):")
  print(length(Cells(scCHANTS[["RNA"]])))
  print(length(Cells(scCHANTS[["HTO"]])))
  print(length(Cells(scCHANTS[["ADT"]])))
  
  
  ## Add metadata
  # read in additional sample metadata
  metadata_samples <- read.csv(paste0(file_path, "sample_metadata.csv"))
  
  # align each cell's HTO id to row in the metadata_samples
  # rownames here represent HTOs but can't have duplicate rownames hence e.g. 4.1
  hto_info <- metadata_samples[match(scCHANTS$HTO_classification, metadata_samples$HTO), ]
  head(hto_info)
  
  # add each row to a new variable in metadata
  scCHANTS$sample_id <- hto_info$sample
  scCHANTS$timepoint <- hto_info$timepoint
  scCHANTS$treatment <- hto_info$treatment
  scCHANTS$PBMC_or_T <- hto_info$PBMC_or_T
  
  # here also reorder HTO levels (get all unique levels, order numerically and factorise HTO_classification metadata column)
  hto_levels <- unique(scCHANTS$HTO_classification)
  hto_levels <- hto_levels[order(as.numeric(sub("HTO-", "", hto_levels)))]
  scCHANTS$HTO_classification <- factor(
    scCHANTS$HTO_classification,
    levels = hto_levels
  )
  
  # print number of cells per sample
  message("Number of cells per sample:")
  print(table(scCHANTS$HTO_classification))
  
  # order timepoints
  timepoint_levels <- c("0", "7", "28", "benchmark")
  scCHANTS$timepoint <- factor(
    scCHANTS$timepoint,
    levels = timepoint_levels
  )
  
  ## Split into PBMCs and T cells
  DefaultAssay(scCHANTS) <- "RNA"
  message("Cell number before splitting:")
  print(length(Cells(scCHANTS@assays$RNA)))
  
  message("PBMC cell number in each assay incl benchmark sample 13:")
  scCHANTS_pbmc <- subset(scCHANTS, subset = PBMC_or_T == "PBMC")
  print(length(Cells(scCHANTS_pbmc@assays$RNA)))
  print(length(Cells(scCHANTS_pbmc@assays$HTO)))
  print(length(Cells(scCHANTS_pbmc@assays$ADT)))
  
  message("T cell number in each assay incl benchmark sample 13:")
  scCHANTS_t <- subset(scCHANTS, subset = PBMC_or_T == "T" | sample_id == "13")
  print(length(Cells(scCHANTS_t@assays$RNA)))
  print(length(Cells(scCHANTS_t@assays$HTO)))
  print(length(Cells(scCHANTS_t@assays$ADT)))
  
  
  ## Rename scCHANTS to designated pbmc or t cell subset
  if(pbmc_or_t == "pbmc") {
    scCHANTS <- scCHANTS_pbmc
    pbmc_or_t_up <- "PBMCs"
  } else if (pbmc_or_t == "t") {
    scCHANTS <- scCHANTS_t
    pbmc_or_t_up <- "T cells"
  }
  
  message(paste0("Subsetted to ", pbmc_or_t_up, " only"))
  
  
  ## Quality control
  # calculate mt contamination and add column to metadata
  scCHANTS[["percent_mt"]] <- PercentageFeatureSet(scCHANTS, pattern = "^MT-")
  
  # calculate mt contamination and add column to metadata
  scCHANTS[["percent_rb"]] <- PercentageFeatureSet(scCHANTS, pattern = "^RP[SL]")
  
  scCHANTS <- subset(scCHANTS, subset = nFeature_RNA >= nFeature_min  & nFeature_RNA <= nFeature_max & percent_mt <= mt_max) 
  
  message(paste0("Number of cells in each assay post-QC using metrics:   ", nFeature_min, " <= nFeature_RNA <= ", nFeature_max, "    and percent_mt   <= ", mt_max))
  print(length(Cells(scCHANTS@assays$RNA)))
  print(length(Cells(scCHANTS@assays$HTO)))
  print(length(Cells(scCHANTS@assays$ADT)))
  
  message("Donor, data and run number added to metadata")
  
  scCHANTS$date <- date
  scCHANTS$run <- run
  scCHANTS$donor <- donor
  
  # remove counts.Antibody Capture layer from RNA assay
  scCHANTS[["RNA"]] <- CreateAssayObject(counts = scCHANTS@assays$RNA$`counts.Gene Expression`)
  message("counts.Antibody Capture removed from RNA assay. RNA assay is now:")
  print(scCHANTS@assays$RNA)
  
  message("...........................................................................................................................................................")
  message(paste0("HTO demux and Quality Control completed for ", pbmc_or_t_up, " for donor ", donor, " in run ", run, " sequenced on ", date))
  
  # return the processed object
  return(scCHANTS)
  
}

## Load data for each donor
scCHANTS_055_pbmc <- demux_qc(date = "20250805", run = "1", donor = "055", pbmc_or_t = "pbmc", remote = "yes", clr_norm = 2, demux_positive = 0.99, nFeature_min = 1100, nFeature_max = 4000, mt_max = 7)

scCHANTS_060_pbmc <- demux_qc(date = "20250805", run = "1", donor = "060", pbmc_or_t = "pbmc", remote = "yes", clr_norm = 2, demux_positive = 0.99, nFeature_min = 1100, nFeature_max = 4000, mt_max = 7)

# for donor 90, demux_positive is 0.95 and nFeature_min is 1000 unlike for other donors
scCHANTS_090_pbmc <- demux_qc(date = "20250805", run = "1", donor = "090", pbmc_or_t = "pbmc", remote = "yes", clr_norm = 2, demux_positive = 0.95, nFeature_min = 1000, nFeature_max = 4000, mt_max = 7)

scCHANTS_096_pbmc <- demux_qc(date = "20250805", run = "1", donor = "096", pbmc_or_t = "pbmc", remote = "yes", clr_norm = 2, demux_positive = 0.99, nFeature_min = 1100, nFeature_max = 4000, mt_max = 7)


## Merge the donors

# merge the four seurat objects, add suffix to cell names
scCHANTS_pbmc <- merge(scCHANTS_055_pbmc, c(scCHANTS_060_pbmc, scCHANTS_090_pbmc, scCHANTS_096_pbmc), add.cell.ids = c("D055", "D060", "D090", "D096"))

# check number of cells is correct
n055 <- Cells(scCHANTS_055_pbmc) %>%  length()
n060 <- Cells(scCHANTS_060_pbmc) %>%  length()
n090 <- Cells(scCHANTS_090_pbmc) %>%  length()
n096 <- Cells(scCHANTS_096_pbmc) %>%  length()
Cells(scCHANTS_pbmc) %>%  length()
n055 + n060 + n090 + n096

# check donor data has merged, not kept as separate layer
Layers(scCHANTS_pbmc@assays$RNA)

# change default assay to RNA
DefaultAssay(scCHANTS_pbmc) <- "RNA"


## When ensured is fully merged etc, remove to save memory
rm(scCHANTS_055_pbmc, scCHANTS_060_pbmc, scCHANTS_090_pbmc, scCHANTS_096_pbmc)



## Normalisation and dimensional reduction


# Do SCTransform on the cluster (takes a few mins) but crashes on normal computer
# vars.to.regress removes the effect of mitochondrial gene percentage as a tehcnical covariate from normalised expression values. Thus each genes expression is adjusted to remove the influence of mitochondrial percentage (as high mt content cells often have biased gene expression profiles) vst v2 is the default, updated version (as of August 2025)

scCHANTS_pbmc <- SCTransform(scCHANTS_pbmc, vars.to.regress = "percent_mt", verbose = FALSE, vst.flavor = "v2")

rna_genes <- rownames(scCHANTS_pbmc[["RNA"]])
length(rna_genes)
sct_genes <- rownames(scCHANTS_pbmc[["SCT"]])
length(sct_genes)

dropped_genes <- setdiff(rna_genes, sct_genes)

# how many genes were filtered out
length(dropped_genes)


# run PCA
scCHANTS_pbmc <- RunPCA(scCHANTS_pbmc, features = VariableFeatures(object = scCHANTS_pbmc))


# if not clear then do this calculation to find where drop in variance explained between two successive PCs no longer significant (is <0.1%)
pct <- scCHANTS_pbmc[["pca"]]@stdev / sum(scCHANTS_pbmc[["pca"]]@stdev) * 100
# calculate cumulative percents for each PC
cumu <- cumsum(pct)
# determine the difference between variation of PC and subsequent PC
# give the last point where difference in % of variation is more than 0.1%
sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
## 11


scCHANTS_pbmc <- FindNeighbors(scCHANTS_pbmc, dims = 1:11) 

scCHANTS_pbmc <- FindClusters(scCHANTS_pbmc, resolution = 0.7)

colnames(scCHANTS_pbmc@meta.data)
Idents(scCHANTS_pbmc) <- "SCT_snn_res.0.7"


## save PBMC pre-processed RDS (takes a long time)

saveRDS(object = scCHANTS_pbmc, file ="/scratch/prj/id_hill_sims_wellcda/scCHANTS/20250805_run1/20250805_scCHANTS_pbmc_processed.Rds")


## Find markers for clusters

# must join layers to do find markers
scCHANTS_pbmc <- JoinLayers(scCHANTS_pbmc)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc_markers <- FindAllMarkers(scCHANTS_pbmc, only.pos = TRUE)

write.csv(pbmc_markers,"/scratch/prj/id_hill_sims_wellcda/scCHANTS/20250805_run1/20250805_scCHANTS_pbmc_markers.csv")

pbmc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
