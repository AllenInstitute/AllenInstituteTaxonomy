######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "Human Multiple Cortical Areas SMART-seq (2019)" (https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq).  Note that this 

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

# ***NOTE***: we also need the dendrogram and I cannot figure out how to read the kind of json on our website into R. I therefore took these from internal sources and have imported them into the AIT file below and have also posted them in the GitHub data folder.

## Read in the data and metadata
sampInfo        <- read.csv("metadata.csv",row.names=1)
taxonomy.counts <- as.matrix(fread("matrix.csv"),rownames=1)
taxonomy.counts <- taxonomy.counts[rownames(sampInfo),]

## Read in the json dendrogram and convert to R dendrogram format
dend <- readRDS("Human_neocortex_SMART_seq_04042025_dend.RData")

## Remove low quality cells from counts
kp = sampInfo$class_label!=""
taxonomy.counts <- counts[kp,]

## Format and subset metadata, and omit the _id, _order, and _color columns (we are NOT retaining orders, except in cluster)
cn <- colnames(sampInfo)
cn <- !(grepl("_color",cn)|grepl("_order",cn)|grepl("_id",cn))
taxonomy.metadata <- sampInfo[kp,cn]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))
taxonomy.metadata <- taxonomy.metadata[,!is.element(colnames(taxonomy.metadata),c("cell_type_alt_alias","full_genotype","outlier_type", "outlier_call","cell_type_alias"))]

## Read in the t-SNE from website. 
tsne <- read.csv("tsne.csv",row.names=1)
tsne <- tsne[rownames(taxonomy.metadata),]

######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## cluster changes to cluster_id below
hierarchy = list("class", "subclass", "cluster_id")  

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = colnames(taxonomy.counts), ncbi.taxid = 9606)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"]             = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor_sex"]           = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="external_donor_name"] = "donor_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="specimen_type"]       = "suspension_type"
taxonomy.metadata$self_reported_ethnicity = "unknown"
taxonomy.metadata$load_id           = "Not reported" 
taxonomy.metadata$assay             = "Smart-seq v4"  
taxonomy.metadata$is_primary_data   = TRUE
taxonomy.metadata$disease           = "control" 
taxonomy.metadata$organism          = "homo sapiens"
taxonomy.metadata$anatomical_region = taxonomy.metadata$region

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE) 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata

## Updating M1 and S1 terms because lower limb region and upper limb region are not in UBERON
old_name <- c("M1lm","M1ul","S1lm","S1ul")
new_name <- c("UBERON:0001384","UBERON:0001384","UBERON:0008933","UBERON:0008933")
taxonomy.anno$anatomical_region_ontology_term_id <- as.character(taxonomy.anno$anatomical_region_ontology_term_id)
for (i in 1:4){
  taxonomy.anno$anatomical_region_ontology_term_id[taxonomy.anno$region==old_name[i]] = new_name[i]
}


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Human_neocortex_SMART_seq_04042025",
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
                            embeddings = tsne, # Use pre-existing t-SNE coordinates
                            ##
                            dend = dend, ## Use the original dendrogram from the website.
							reorder.dendrogram = TRUE,  ## Reorder leaf nodes to try and match original dendrogram
                            taxonomyDir = taxonomyDir, ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_neocortex_SMART_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Human_neocortex_SMART_seq_04042025.h5ad")

