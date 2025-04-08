######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "V1 & ALM - SMART-seq (2018)"  (https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq).

# This code was run within the scrattch docker environment using docker://alleninst/scrattch:1.1.2 and may produde different results if run in any other environment.


######################################################################################
print("========== Prepare our working environment ==========")

## Prior to running the code, navigate to the desired working directory
taxonomyDir <- getwd()

## Load the libraries
suppressPackageStartupMessages({
  library(scrattch.taxonomy)
  library(reticulate)
})
reticulate::use_python("/usr/bin/python3")
cell_type_mapper <- import("cell_type_mapper") # For hierarchical mapping
set.seed(42)


######################################################################################
print("========== Download and read in the reference taxonomy from Allen Brain Map ==========")
  
## Download the data
# Prior to running the code, download this file and upzip in your working directory: https://celltypes.brain-map.org/api/v2/well_known_file_download/694416044

# ***NOTE***: we also need the dendrogram and t-SNE coordinates, neither of which I could find publicly. I took these from internal sources and have imported them into the AIT file below and have also posted them in the GitHub data folder here: 

# ***NOTE***: two of the cluster names have changed between publication of the data and what was used in this paper. For this AIT file we CHANGE the cluster names to align with what is in the paper.

## Read in the VISp data
exons_VISp    <- as.matrix(fread("data/mouse_VISp_2018-06-14_exon-matrix.csv"),rownames=1)
introns_VISp  <- as.matrix(fread("data/mouse_VISp_2018-06-14_intron-matrix.csv"),rownames=1)
geneInfo      <- read.csv("data/mouse_VISp_2018-06-14_genes-rows.csv",row.names=1)
sampInfo_VISp <- read.csv("data/mouse_VISp_2018-06-14_samples-columns.csv",row.names=1)
kp_VISp       <- !is.element(sampInfo_VISp$class,c("Low Quality","No Class"))  ## Identify cells that failed QC.

## Read in the VISp data
exons_ALM    <- as.matrix(fread("data/mouse_ALM_2018-06-14_exon-matrix.csv"),rownames=1)
introns_ALM  <- as.matrix(fread("data/mouse_ALM_2018-06-14_intron-matrix.csv"),rownames=1)
sampInfo_ALM <- read.csv("data/mouse_ALM_2018-06-14_samples-columns.csv",row.names=1)
kp_ALM       <- !is.element(sampInfo_ALM$class,c("Low Quality","No Class"))  ## Identify cells that failed QC.


## Format and subset counts, merging VISp and ALM
taxonomy.counts <- cbind((introns_VISp + exons_VISp),(introns_ALM + exons_ALM))
rownames(taxonomy.counts) <- rownames(geneInfo)
taxonomy.counts <- taxonomy.counts[, c(kp_VISp,kp_ALM)]  # Omit cells from outlier clusters as above

## Format and subset metadata
taxonomy.metadata <- rbind(sampInfo_VISp[kp_VISp,],sampInfo_ALM[kp_ALM,])

## Read in and format the t-SNE from internal sources. 
# NOTE: There are 23 cells with missing tsne coordinates, which we will set to 0 for this file

tsne_feather      <- feather::feather("Mouse_VISp_ALM_SMART_seq_04042025_tsne.feather")
tsne_feather      <- as.data.frame(tsne_feather)
kp_tsne           <- !is.na(tsne_feather[,1])
tsne              <- tsne_feather[kp_tsne,2:3]
rownames(tsne)    <- tsne_feather[kp_tsne,1]
tsne              <- tsne[match(rownames(tsne),(taxonomy.metadata$seq_name)),]
rownames(tsne)    <- rownames(taxonomy.metadata)
tsne[is.na(tsne)] <- 0

## Read the dendrogram from internal sources
dend <- readRDS("Mouse_VISp_ALM_SMART_seq_04042025_dend.RData")

## Rename cluster to cluster_id to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"] = "cluster_id"
taxonomy.metadata$cluster_id <- as.character(taxonomy.metadata$cluster_id)

## Rename two clusters in the taxonomy.metadata to match the dendrogram (and paper).
taxonomy.metadata$cluster_id[taxonomy.metadata$cluster_id=="L6 CT ALM Nxph2 Sla"]  = "L6 CT Nxph2 Sla"
taxonomy.metadata$cluster_id[taxonomy.metadata$cluster_id=="L6b VISp Col8a1 Rprm"] = "L6b Col8a1 Rprm"
taxonomy.metadata$cluster_id <- factor(taxonomy.metadata$cluster_id,levels=labels(dend))


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## -- See note about subclass above
hierarchy = list("class", "subclass", "cluster_id")

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(taxonomy.counts), ncbi.taxid = 10090)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="sex"]         = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor"]       = "donor_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="sample_type"] = "suspension_type"
taxonomy.metadata$load_id           = taxonomy.metadata$seq_batch  # Setting load_id as seq_batch to match Hodge et al
taxonomy.metadata$assay             = "Smart-seq v4"  
taxonomy.metadata$anatomical_region = taxonomy.metadata$brain_region
taxonomy.metadata$is_primary_data   = TRUE
taxonomy.metadata$disease           = "control" 

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE) 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms


### NOTE I NEED TO UPDATE THE ABOVE CALL TO LOOK IN THE MOUSE BRAIN ONOLOGY AND NOT THE HUMAN ONE!!!


## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Mouse_VISp_ALM_SMART_seq_04042025",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1200, ## Select top 1200 binary genes for consistency with Hodge et al 2019
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
							gene.meta.data = geneInfo, # Additional gene information
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = tsne, # Use pre-existing t-SNE coordinates
                            ##
                            dend = dend, ## Pre-computed dendrogram
                            taxonomyDir = taxonomyDir, ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Mouse_VISp_ALM_SMART_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Mouse_VISp_ALM_SMART_seq_04042025.h5ad")

