######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "Whole Cortex & Hippocampus - SMART-seq (2020) with 10x-SMART-seq taxonomy (2021)" (https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x).

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

# ***NOTE***: The published dendrogram is not a binary tree and therefore is incompatible with scrattch.taxonomy at this time. Therefore we are creating a new tree for the AIT file.  Note that this tree will also NOT match the tree created in Mouse_cortex_hippocampus_10X_seq_04042025

## Read in the data and metadata
sampInfo        <- read.csv("metadata.csv",row.names=1)
taxonomy.counts <- as.matrix(fread("matrix.csv"),rownames=1)
taxonomy.counts <- taxonomy.counts[rownames(sampInfo),]

## Read in the UMAP from website. 
umap <- read.csv("tsne.csv",row.names=1)
umap <- umap[rownames(sampInfo),]

## Format and subset metadata, and omit the _id, _order, and _color columns (we are NOT retaining orders, except in cluster)
cn <- colnames(sampInfo)
cn <- !(grepl("_color",cn)|grepl("_order",cn)|grepl("_id",cn))
taxonomy.metadata <- sampInfo[,cn]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))
taxonomy.metadata <- taxonomy.metadata[,!is.element(colnames(taxonomy.metadata),c("cell_type_alt_alias","injection_type","platform"))]


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## -- NOTE: We exclude neighborhood from the hierarchy; otherwise, mapmycells cannot be set up correctly (I don't know why...)
hierarchy = list("class", "subclass", "cluster_id")  


## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(taxonomy.counts), ncbi.taxid = 10090)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"]             = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor_sex"]           = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="external_donor_name"] = "donor_id"
taxonomy.metadata$load_id           = "Not reported" 
taxonomy.metadata$assay             = "Smart-seq v4"  
taxonomy.metadata$anatomical_region = taxonomy.metadata$region
taxonomy.metadata$is_primary_data   = TRUE
taxonomy.metadata$disease           = "control" 
taxonomy.metadata$suspension_type   = "cells" 
taxonomy.metadata$organism          = "mus musculus"

## Five brain regions do not line up with MBA (mostly due to them being combinations of regions).  
#  For these we rename with the lowest-resolution term I can figure out
old_name <- c("ALM","PL-ILA","TEa-PERI-ECT","PAR-POST-PRE","SUB-ProS")
new_name <- c("MO","Isocortex","Isocortex","RHP","RHP")
for (i in 1:5){
 taxonomy.metadata$anatomical_region[taxonomy.metadata$anatomical_region==old_name[i]] = new_name[i]
}

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE, compute.brain.atlas.terms = "MBA") 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata

## Several cell types in MBA do not align with UBERON terms in the current UBERON version I am using
# ---- Can be revisited in future updates
taxonomy.anno <- taxonomy.anno[,colnames(taxonomy.anno)!="anatomical_region_ontology_term_id"]


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Mouse_cortex_hippocampus_SMART_seq_04042025",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1200, ## Select top 1200 binary genes for consistency with Hodge et al 2019
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
							gene.meta.data = NULL, # No additional gene information
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = umap, # Use pre-existing UMAP coordinates
                            ##
                            dend = "highly_variable_genes_standard", ## Create new dendrogram for this file
							reorder.dendrogram = TRUE,  ## Reorder leaf nodes to try and match original dendrogram
                            taxonomyDir = taxonomyDir, ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Mouse_cortex_hippocampus_SMART_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Mouse_cortex_hippocampus_SMART_seq_04042025.h5ad")

