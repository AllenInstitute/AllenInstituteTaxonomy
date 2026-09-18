######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "Mouse Comparative LGN (2018)"  (https://portal.brain-map.org/atlases-and-data/rnaseq/comparative-lgn).

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

# *** NOTE ***: we are losing a few cells because the data and metadata do not perfectly match!

## Read in the data
exons    <- as.matrix(fread("mouse_LGN_2021_exon-matrix.csv"),rownames=1)
introns  <- as.matrix(fread("mouse_LGN_2021_intron-matrix.csv"),rownames=1)
geneInfo <- read.csv("mouse_LGN_2021_genes-rows.csv",row.names=1)
sampInfo <- read.csv("mouse_LGN_2021_metadata.csv",row.names=2)

## Identify and remove cells from the "Low Quality" cluster.
kp <- sampInfo$cluster_label!="Low Quality"

## Format and subset counts
taxonomy.counts <- introns + exons
taxonomy.counts <- taxonomy.counts[,kp]  # Omit cells from outlier clusters as above
rownames(taxonomy.counts) <- rownames(geneInfo)

## Format and subset metadata, and omit the _id and _color columns (we are NOT retaining orders)
cn <- colnames(sampInfo)
cn <- (!(grepl("_color",cn)|grepl("_id",cn)|is.element(cn,c("X","species_label","gender.y","roi"))))|(cn=="sample_id")
taxonomy.metadata <- sampInfo[kp,cn]
colnames(taxonomy.metadata) <- gsub("_label","",colnames(taxonomy.metadata))

## Align metadata and counts (there are different cells in both sets of files!!!)
kpCells           <- intersect(rownames(taxonomy.metadata),colnames(taxonomy.counts))
taxonomy.metadata <- taxonomy.metadata[kpCells,]
taxonomy.counts   <- taxonomy.counts[,kpCells]


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## -- Two level hierachy; note that cluster --> cluster_id below
hierarchy = list("class", "cluster_id")

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
# -- Note, we are using rhesus macaque here because the genome is better annotated.
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(taxonomy.counts), ncbi.taxid = 10090)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"] = "cluster_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="gender.x"]= "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor"]   = "donor_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="batch"]   = "load_id"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cell_prep_type"] = "suspension_type"
taxonomy.metadata$assay           = "Smart-seq v4"  
taxonomy.metadata$is_primary_data = TRUE
taxonomy.metadata$disease         = "control" 
taxonomy.metadata$organism        = "Mus musculus" 
taxonomy.metadata$anatomical_region  = "Dorsal part of the lateral geniculate complex"  # This will get overwritten below
taxonomy.metadata$brain_region       = "Dorsal part of the lateral geniculate complex"
taxonomy.metadata$self_reported_sex[taxonomy.metadata$self_reported_sex==""] = "unknown"

## Check taxonomy metadata aligns with AIT standard and perform minor error corrections
## Also add ontology terms corresponding to the above schema elements (and can also correct misspellings, etc.)
full.taxonomy.anno <- computeOntologyTerms(taxonomy.metadata, standardize.metadata=TRUE, print.messages=TRUE, compute.brain.atlas.terms = "MBA") 
# NOTE: We encourage reviewing the messages from this function CAREFULLY, as some assumptions are made when calculating ontology terms

## Save final metadata data frame
taxonomy.anno <- full.taxonomy.anno$metadata

## Update UBERON numbers and anatomical_region for different cellular layers
roi    = c("dLGN - Core","dLGN - Shell", "LGv", "LP")
name   = c("dorsal lateral geniculate nucleus",
           "dorsal lateral geniculate nucleus",
		   "ventral lateral geniculate nucleus",
		   "lateral posterior nucleus of thalamus")
uberon = c("UBERON:0002479","UBERON:0002479","UBERON:0002480","UBERON:0002983")
taxonomy.anno$anatomical_region_ontology_term_id <- as.character(taxonomy.anno$anatomical_region_ontology_term_id)
taxonomy.anno$anatomical_region <- as.character(taxonomy.anno$anatomical_region)
for (i in 1:4){
  kp = taxonomy.anno$roi==roi[i]
  taxonomy.anno$anatomical_region_ontology_term_id[kp] = uberon[i]
  taxonomy.anno$anatomical_region[kp] = name[i]
}


######################################################################################
print("========== Create the (parent) AIT Taxonomy ==========")

## Build Allen Insitute Taxonomy, for large taxonomies you can pass in tpm and cluster_stats if pre-computed.
AIT.anndata = buildTaxonomy(title="Mouse_LGN_SMART_seq_04042025",
                            meta.data = taxonomy.anno,
                            hierarchy = hierarchy,
                            ## --- Optional parameters ---
                            counts = taxonomy.counts,
                            normalized.expr = NULL,
                            highly_variable_genes = 1000, ## Select top 1000 binary genes for consistency with Hodge et al 2019
                            marker_genes = NULL,
                            ensembl_id = ensembl_id,
							gene.meta.data = NULL, # Not including additional gene information here
                            cluster_stats = NULL, ## Pre-computed cluster stats
                            embeddings = "highly_variable_genes_standard", # Calculate UMAP based on binary genes
                            ##
                            dend = "highly_variable_genes_standard", # Calculate dendrogram based on binary genes
                            taxonomyDir = taxonomyDir, ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Mouse_LGN_SMART_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Mouse_LGN_SMART_seq_04042025.h5ad")

