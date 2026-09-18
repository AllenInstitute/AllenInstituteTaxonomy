######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "Human Comparative LGN (2018)"  (https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn).

# This code was run within the scrattch docker environment using docker://jeremyinseattle/scrattch:1.1.2 and may produde different results if run in any other environment.


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

## Read in the data
exons    <- as.matrix(fread("human_LGN_2021_exon-matrix.csv"),rownames=1)
introns  <- as.matrix(fread("human_LGN_2021_intron-matrix.csv"),rownames=1)
geneInfo <- read.csv("human_LGN_2021_genes-rows.csv",row.names=1)
sampInfo <- read.csv("human_LGN_2021_metadata.csv",row.names=2)

## Format and subset counts
taxonomy.counts <- introns + exons
rownames(taxonomy.counts) <- rownames(geneInfo)

## Omit the _id and _color columns (we are NOT retaining orders)
cn <- colnames(sampInfo)
cn <- (!(grepl("_color",cn)|grepl("_id",cn)|is.element(cn,c("X","fastq_file_name","dissection_roi"))))|(cn=="sample_id")
taxonomy.metadata <- sampInfo[,cn]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## -- Two level hierachy; note that cluster --> cluster_id below
hierarchy = list("class", "cluster_id")

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(taxonomy.counts), ncbi.taxid = 9606)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"] = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="gender"]  = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor"]   = "donor_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="batch"]   = "load_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="species"] = "organism"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cell_prep_type"] = "suspension_type"
taxonomy.metadata$assay           = "Smart-seq v4"  
taxonomy.metadata$is_primary_data = TRUE
taxonomy.metadata$disease         = "control" 
taxonomy.metadata$organism        = "homo sapiens" 
taxonomy.metadata$self_reported_ethnicity = "unknown"
taxonomy.metadata$anatomical_region  = "lateral geniculate nucleus"
taxonomy.metadata$brain_region       = "lateral geniculate nucleus"

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE) 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata

## Update UBERON numbers and anatomical_region for different cellular layers
roi    = c("KC","MC","PC")
name   = c("koniocellular layer of dorsal nucleus of lateral geniculate body",
           "magnocellular layer of dorsal nucleus of lateral geniculate body",
		   "parvocellular layer of dorsal nucleus of lateral geniculate body")
uberon = c("UBERON:0013615","UBERON:0013606","UBERON:0002479")
taxonomy.anno$anatomical_region_ontology_term_id <- as.character(taxonomy.anno$anatomical_region_ontology_term_id)
taxonomy.anno$anatomical_region <- as.character(taxonomy.anno$anatomical_region)
for (i in 1:3){
  kp = taxonomy.anno$roi==roi[i]
  taxonomy.anno$anatomical_region_ontology_term_id[kp] = uberon[i]
  taxonomy.anno$anatomical_region[kp] = name[i]
}


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Human_LGN_SMART_seq_04042025",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1000, ## Select top 1000 binary genes for consistency with Hodge et al 2019
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
							gene.meta.data = geneInfo, # Additional gene information
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = "highly_variable_genes_standard", # Calculate UMAP based on binary genes
                            ##
                            dend = "highly_variable_genes_standard", # Calculate dendrogram based on binary genes
                            taxonomyDir = taxonomyDir, ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## NOTE: This taxonomy FAILS to create statistics for MapMyCells, which means that only correlation and Seurat mapping will work.  Given that this is not a taxonomy folks will likely use anymore, I'm not going to worry about this.							
							
## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))

######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_LGN_SMART_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Human_LGN_SMART_seq_04042025.h5ad")

