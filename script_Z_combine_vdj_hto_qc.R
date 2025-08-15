# extract hto - cell barcode links from seurat object and merge with vdj csv

library("dplyr")
library("Seurat")



# function to do processing, without qc, just hto demultiplexing, does not filter singlets
demux <- function(date, run, donor, remote = "no", clr_norm = 2, demux_positive = 0.99) {
  
  # check all args are valid
  stopifnot(
    remote %in% c("yes", "no"),
    clr_norm %in% c(1, 2),
    demux_positive >= 0 & demux_positive <= 0.99,
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
  #scCHANTS$totRNA <- Matrix::colSums(scCHANTS@assays$RNA@layers$`counts.Gene Expression`)
  #scCHANTS$totHTO <- Matrix::colSums(scCHANTS@assays$HTO@counts)
  
  #scCHANTS <- subset(scCHANTS, subset = totRNA > 0 | totHTO > 0)
  
  DefaultAssay(scCHANTS) <- "HTO"
  
  # normalise HTO counts, standard is centered log-ratio (CLR) transformation
  # can normalise across features (margin = 1) or across cells (margin = 2)
  # comment on github forum satija lab says they recommend normalising across tags if variation between hash performances but since have very low signals, will try across cells
  scCHANTS <- NormalizeData(scCHANTS, assay = "HTO", normalization.method = "CLR", margin = clr_norm)
  
  scCHANTS <- HTODemux(scCHANTS, assay = "HTO", positive.quantile = demux_positive, kfunc = "clara")
  
  DefaultAssay(scCHANTS) <- "HTO"
  
  message("Original number of cells:")
  print(length(Cells(scCHANTS)))
  
  #scCHANTS <- subset(scCHANTS, subset = HTO_classification.global == "Singlet")
  
  #message("Number of cells after subsetting for singlets in RNA, HTO and ADT assays (should all be the same):")
  #print(length(Cells(scCHANTS[["RNA"]])))
  #print(length(Cells(scCHANTS[["HTO"]])))
  #print(length(Cells(scCHANTS[["ADT"]])))
  
  
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
  
  message("...........................................................................................................................................................")
  message(paste0("HTO demux completed for donor ", donor, " in run ", run, " sequenced on ", date))
  
  # return the processed object
  return(scCHANTS)
  
}

scCHANTS_055 <- demux(date = "20250805", run = "1", donor = "055", remote = "no", clr_norm = 2, demux_positive = 0.99)

# extract hto classification
# can add extra metadata columns if you want. may have to adjust the is.nas for that column as well when joining to vdj
hto_calls <- data.frame(
  barcode = colnames(scCHANTS_055),
  HTO_classification.global = scCHANTS_055$HTO_classification.global,
  HTO_classification = scCHANTS_055$HTO_classification,
  sample_id = scCHANTS_055$sample_id
)
head(hto_calls)



# load vdj data

vdj <- read.csv("/Users/k2477939/Library/CloudStorage/OneDrive-King\'sCollegeLondon/General\ -\ Laboratory\ of\ Molecular\ Immunotherapy\ and\ Antibiotics/Data/Georgia\ Miller/scCHANTS/filtered_contig_annotations.csv")
colnames(vdj)
head(vdj)



# add hto to vdj data
vdj_annot <- vdj %>%
  left_join(hto_calls, by = "barcode") %>% 
  # doublet sample_ids will be set as NA, change to "doublet"
  mutate(sample_id = ifelse(is.na(sample_id), "doublet", sample_id))

colnames(vdj_annot)
head(vdj_annot)




# now extract qc status from pbmcs and t cells
# function to do processing, split into pbmc or t cells and qc as well. does not filter by singlets or qc
demux_qc <- function(date, run, donor, remote = "no", clr_norm = 2, demux_positive = 0.99, nFeature_min, nFeature_max, mt_max, do_QC = FALSE) {
  
  # check all args are valid
  stopifnot(
    remote %in% c("yes", "no"),
    clr_norm %in% c(1, 2),
    demux_positive >= 0 & demux_positive <= 0.99,
    is.numeric(nFeature_min), length(nFeature_min) == 1,
    is.numeric(nFeature_max), length(nFeature_max) == 1,
    is.numeric(mt_max), length(mt_max) == 1,
    # check donor has 3 digits
    grepl("^[0-9]{3}$", donor),
    do_QC %in% c(TRUE, FALSE)
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
  
  
  ## Quality control
  # calculate mt contamination and add column to metadata
  scCHANTS[["percent_mt"]] <- PercentageFeatureSet(scCHANTS, pattern = "^MT-")
  
  # calculate mt contamination and add column to metadata
  scCHANTS[["percent_rb"]] <- PercentageFeatureSet(scCHANTS, pattern = "^RP[SL]")
  
  # add qc metrics to metadata
  # or assign new metadata column
  # start with all cells labeled as "initial"
  scCHANTS$QC <- "initial"
  
  # label low nFeature_RNA (< nFeature_min)
  ## i.e. if nFeature_RNA < nFeature_min (and previously passed QC as cautionary step), call it low_nFeature, otherwise leave as initial
  
  scCHANTS$QC <- ifelse(
    scCHANTS$nFeature_RNA < nFeature_min & scCHANTS$QC == "initial",
    "low_nFeature",
    scCHANTS$QC
  )
  
  # label high nFeature_RNA (> nFeature_max)
  scCHANTS$QC <- ifelse(
    scCHANTS$nFeature_RNA > nFeature_max & scCHANTS$QC == "initial",
    "high_nFeature",
    scCHANTS$QC
  )
  
  # label high mitochondrial content (> mt_max%)
  scCHANTS$QC <- ifelse(
    scCHANTS$percent_mt > mt_max & scCHANTS$QC == "initial",
    "high_mt",
    scCHANTS$QC
  )
  
  # label combined: low_nFeature + high_mt
  scCHANTS$QC <- ifelse(
    scCHANTS$nFeature_RNA < nFeature_min & scCHANTS$percent_mt > mt_max & scCHANTS$QC != "pass",
    "low_nFeature_high_mt",
    scCHANTS$QC
  )
  
  # label combined: high_nFeature + high_mt
  scCHANTS$QC <- ifelse(
    scCHANTS$nFeature_RNA > nFeature_max & scCHANTS$percent_mt > mt_max & scCHANTS$QC != "pass",
    "high_nFeature_high_mt",
    scCHANTS$QC
  )
  
  # finally, label cells that pass all filters
  scCHANTS$QC <- ifelse(
    scCHANTS$nFeature_RNA >= nFeature_min &
      scCHANTS$nFeature_RNA <= nFeature_max &
      scCHANTS$percent_mt <= mt_max &
      scCHANTS$QC == "initial",
    "pass",
    scCHANTS$QC
  )
  
  
  if(do_QC == TRUE) {
    scCHANTS <- subset(scCHANTS, subset = QC == "pass") 
    message("Data subsetted to QC = pass")
    message(paste0("Number of cells in each assay post-QC using metrics:   ", nFeature_min, " <= nFeature_RNA <= ", nFeature_max, "    and percent_mt   <= ", mt_max))
    print(length(Cells(scCHANTS@assays$RNA)))
    print(length(Cells(scCHANTS@assays$HTO)))
    print(length(Cells(scCHANTS@assays$ADT)))
  }
  
  message("QC metric, donor, data and run number added to metadata")
  
  scCHANTS$date <- date
  scCHANTS$run <- run
  scCHANTS$donor <- donor
  
  # remove counts.Antibody Capture layer from RNA assay
  scCHANTS[["RNA"]] <- CreateAssayObject(counts = scCHANTS@assays$RNA$`counts.Gene Expression`)
  message("counts.Antibody Capture removed from RNA assay. RNA assay is now:")
  print(scCHANTS@assays$RNA)
  
  message("...........................................................................................................................................................")
  message(paste0("HTO demux completed and Quality Control added for donor ", donor, " in run ", run, " sequenced on ", date))
  
  # return the processed object
  return(scCHANTS)
  
}

scCHANTS_055_qc <- demux_qc(date = "20250805", run = "1", donor = "055", remote = "no", clr_norm = 2, demux_positive = 0.99, nFeature_min = 1100, nFeature_max = 4000, mt_max = 7, do_QC = FALSE)

# to look at the metadata
metadata <- scCHANTS_055_qc@meta.data
colnames(metadata)
head(metadata)

# extract QC classification
# can add extra metadata columns if you want. may have to adjust the is.nas for that column as well when joining to vdj
qc_calls <- data.frame(
  barcode = colnames(scCHANTS_055_qc),
  QC = scCHANTS_055_qc$QC
)
head(qc_calls)



# add qc to vdj data
vdj_annot2 <- vdj_annot %>%
  left_join(qc_calls, by = "barcode") %>% 
  # doublet sample_ids will be set as NA, change to "doublet"
  mutate(QC = ifelse(is.na(QC), "prev_filtered", QC))
colnames(vdj_annot2)
head(vdj_annot2)


# filter vdj
vdj_filtered <- vdj_annot2 %>% 
  filter(HTO_classification.global != "Doublet",
         QC == "pass")

# export
#write.csv(vdj_filtered, "/Users/k2477939/Library/CloudStorage/OneDrive-King\'sCollegeLondon/General\ -\ Laboratory\ of\ Molecular\ Immunotherapy\ and\ Antibiotics/Data/Georgia\ Miller/scCHANTS/singlet_qc_filtered_contig_annotations_055.csv")

