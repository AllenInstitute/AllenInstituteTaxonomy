######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "Human MTG SMART-seq (2018)"  (https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq).

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
# Prior to running the code, download several files from https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x to your working directory ("human_dendrogram.rds","matrix.csv","metadata.csv", "tsne.csv")

## ***NOTE***: we are omitting some of the CCN information related to higher levels of the taxonomy (subclass and cluster), which are available as part of the other taxonomy files at https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x.  These may be rescued in later versions of the file.

## Read in the data
taxonomy.counts   <- as.matrix(fread("matrix.csv"),rownames=1)
taxonomy.metadata <- read.csv("metadata.csv",row.names=1)

## Format metadata
# -- Note, we are removing the colors form _color and embedding the factor information from _order, and then dropping those columns
cn_in <- colnames(taxonomy.metadata)[grepl("_label",colnames(taxonomy.metadata))][1:11] # Omit full_genotype
for(i in 1:8) {  # No need to factorize all the CCN stuff
  ids   <- taxonomy.metadata[,gsub("_label","_order",cn_in[i])]
  level <- taxonomy.metadata[,cn_in[i]][match(sort(unique(ids)),ids)]
  taxonomy.metadata[,cn_in[i]] <- factor(taxonomy.metadata[,cn_in[i]],levels=level)
}
taxonomy.metadata <- taxonomy.metadata[,cn_in]
colnames(taxonomy.metadata) <- gsub("_label","",cn_in)

## Read in and format the t-SNE from website. 
tsne <- read.csv("tsne.csv",row.names=1)

## Read the dendrogram from website
dend <- readRDS("human_dendrogram.rds")

## Omit "Oligo L3-6 OPALIN LRP4-AS1_Excluded" from the dendrogram
dend <- prune(dend, "Oligo L3-6 OPALIN LRP4-AS1_Excluded")
sum(labels(dend)!=levels(taxonomy.metadata$cluster))
# [1] 0  # Confirm matching taxonomy and dendrogram levels and names


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
hierarchy = list("class", "subclass", "cluster_id")  # cluster --> cluster_id below

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = colnames(taxonomy.counts), ncbi.taxid = 9606)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"]             = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor_sex"]           = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="external_donor_name"] = "donor_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="region"]              = "anatomical_region"
taxonomy.metadata$load_id           = "Not reported"
taxonomy.metadata$assay             = "10x 3' v3"  
taxonomy.metadata$suspension_type   = "nucleus"
taxonomy.metadata$is_primary_data   = TRUE
taxonomy.metadata$self_reported_ethnicity = "unknown"
taxonomy.metadata$disease           = "control" 
taxonomy.metadata$organism          = "Homo sapiens"

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE) 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata

## Update UBERON and standard terms for "primary motor cortex" since the translation is incorrect
taxonomy.anno$anatomical_region_ontology_term_id = "UBERON:0001384"
taxonomy.anno$brain_region      =  taxonomy.anno$anatomical_region
taxonomy.anno$anatomical_region = "primary motor cortex"


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Human_M1_10X_seq_04042025",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1000, ## Select top 1000 binary genes for consistency with Hodge et al 2019
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
							gene.meta.data = NULL, # Additional gene information
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = tsne, # Use pre-existing t-SNE coordinates
                            ##
                            dend = dend, ## Pre-computed dendrogram
                            taxonomyDir = getwd(), ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_M1_10X_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Human_M1_10X_seq_04042025.h5ad")

