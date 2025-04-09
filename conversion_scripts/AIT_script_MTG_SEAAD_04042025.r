######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "MTG - 10x SEA-AD (2022)  (https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad).

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
  
# Read data in but do not subset (since we are making data public).  We note that "cluster_label", "subclass_label", and "class_label" correspond to SEA-AD supertype, subclass, and class, respectively, and are used for defining the hierarchy.  

## Download the reference data to the working directory and read it in
seaad_url  <- "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/MTG/RNAseq/Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad"
dend_url   <- "https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0f/37/0f3755cb-3acb-4b93-8a62-5d6adc74c673/dend.rds"
#download.file(seaad_url,"Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")  # NOTE: we recommend downloading via the web browser, as this command may fail
#download.file(dend_url,"Reference_MTG_dend.rds")
seaad_data <- read_h5ad("Reference_MTG_RNAseq_final-nuclei.2022-06-07.h5ad")
seaad_dend <- readRDS("Reference_MTG_dend.rds")

## Read in the UMAP from website. 
umap <- read.csv("umap.csv",row.names=1)

## Get data and annotations
taxonomy.counts = seaad_data$X
cn <- c("sample_name","cluster_label","cluster_confidence","subclass_label","class_label",
        "external_donor_name_label","age_label","donor_sex_label")
taxonomy.metadata = seaad_data$obs[,cn]

## Ensure count matrix and annotations are in the same order (this shouldn't be needed)
taxonomy.metadata = taxonomy.metadata[match(rownames(taxonomy.counts), taxonomy.metadata$sample_name),]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## -- This MUST be from broadest to most specific types, and NOT vice versa
## -- Note that "cluster" is the SEAAD supertypes and will be named "cluster_id" below
hierarchy = list("class", "subclass", "cluster_id")

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = colnames(taxonomy.counts), ncbi.taxid = 9606)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"]             = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor_sex"]           = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="external_donor_name"] = "donor_id"
taxonomy.metadata$load_id           = "Not reported"
taxonomy.metadata$assay             = "10x 3' v3"  
taxonomy.metadata$organism          = "Homo sapiens"
taxonomy.metadata$anatomical_region = "Middle temporal gyrus"
taxonomy.metadata$suspension_type   = "nucleus"
taxonomy.metadata$is_primary_data   = TRUE
taxonomy.metadata$self_reported_ethnicity = "unknown"
taxonomy.metadata$disease           = "control" 

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE) 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Human_MTG_SEAAD_04042025",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1000, ## Select top 1000 binary genes
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = umap, ## Use the precomputed UMAP coordinates
                            ##
                            dend = seaad_dend, ## Pre-computed dendrogram
                            taxonomyDir = getwd(), ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_MTG_SEAAD_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Human_MTG_SEAAD_04042025.h5ad")

