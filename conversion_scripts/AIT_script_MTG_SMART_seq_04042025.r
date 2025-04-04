######################################################################################
## Overview

# This script describes how we created an AIT taxonomy (scrattch v1.1) of the taxonomy for "Human MTG SMART-seq (2018)"  (https://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-smart-seq).

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

# ***NOTE***: we also need the dendrogram and t-SNE coordinates, neither of which I could find publicly. I took these from internal sources and have imported them into the AIT file below and have also posted them in the GitHub data folder here: 

# ***NOTE***: some of the cluster names have changed between publication of the data and what was used in this paper. For this AIT file we CHANGE the cluster names to align with what is in the paper.

# ***NOTE***: We have slotted "pCL name" into subclass since this project predates the formal concept of >2 taxonomy levels.

## Read in the data
exons    <- as.matrix(fread("human_MTG_2018-06-14_exon-matrix.csv"),rownames=1)
introns  <- as.matrix(fread("human_MTG_2018-06-14_intron-matrix.csv"),rownames=1)
geneInfo <- read.csv("human_MTG_2018-06-14_genes-rows.csv",row.names=1)
sampInfo <- read.csv("human_MTG_2018-06-14_samples-columns.csv",row.names=1)

## Identify cells with no class. These failed QC.
kp <- sampInfo$cluster!="no class"

## Format and subset counts
taxonomy.counts <- introns + exons
rownames(taxonomy.counts) <- rownames(geneInfo)
taxonomy.counts <- taxonomy.counts[,kp]  # Omit cells from outlier clusters as above

## Format and subset metadata
taxonomy.metadata <- sampInfo[kp,]

## Read in and format the t-SNE from internal sources. 
# NOTE: There are 337 cells with missing tsne coordinates, which we will set to 0 for this file

tsne <- read.csv("Human_MTG_SMART_seq_04042025_tsne.csv",row.names=1)
tsne <- tsne[match(rownames(tsne),(taxonomy.metadata$seq_name)),]
rownames(tsne)    <- rownames(taxonomy.metadata)
tsne[is.na(tsne)] <- 0

## Read the dendrogram from internal sources
load("Human_MTG_SMART_seq_04042025_dend.RData")

## Extract subclass and original cluster names from Supplementary Table 3 in Hodge et al 2019, available here:  https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1506-7/MediaObjects/41586_2019_1506_MOESM5_ESM.xlsx
# I have downloaded this and converted it into a csv called "Human_MTG_SMART_seq_04042025_cluster_info.csv" which is also now shared in the above GitHub folder

cluster_info <- read.csv("Human_MTG_SMART_seq_04042025_cluster_info.csv", row.names=1)
all_clusters <- cluster_info$transcriptome_data_cluster[1:75]
subclass     <- cluster_info$`pCL_name..or.CL_name.`[match(cluster_info[1:75,'is_a..CL.or.pCL_id.'],rownames(cluster_info))]

## Change cluster names in annotations and dendrograms to align with paper
old_names <- setdiff(labels(dend),all_clusters)
new_names <- setdiff(all_clusters,labels(dend))

colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="cluster"] = "cluster_id"
taxonomy.metadata$cluster_id <- as.character(taxonomy.metadata$cluster_id)
for (i in 1:length(old_names)){
  labels(dend)[labels(dend)==old_names[i]] = new_names[i]
  taxonomy.metadata$cluster_id[taxonomy.metadata$cluster_id==old_names[i]] = new_names[i]
}
taxonomy.metadata$cluster_id <- factor(taxonomy.metadata$cluster_id,levels=all_clusters)

## Add the subclass
subclass <- setNames(subclass,all_clusters)
taxonomy.metadata$subclass <- factor(as.character(subclass[taxonomy.metadata$cluster_id]), levels=unique(subclass))


######################################################################################
print("========== Align metadata to AIT Taxonomy ==========")

## Set up the levels of hierarchy for all mapping functions later
## -- See note about subclass above
hierarchy = list("class", "subclass", "cluster_id")

## Identify Ensembl IDs 
# Common NCBI taxIDs: Human = 9606; Mouse = 10090; Macaque (rhesus) = 9544; Marmoset = 9483
ensembl_id <- geneSymbolToEnsembl(gene.symbols = rownames(taxonomy.counts), ncbi.taxid = 9606)

## Update the metadata to align with AIT schema
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="sex"]     = "self_reported_sex"
colnames(taxonomy.metadata)[colnames(taxonomy.metadata)=="donor"]   = "donor_id"
taxonomy.metadata$load_id           = taxonomy.metadata$seq_batch  # Setting load_id as seq_batch since it's the lowest resolution batch that I can find
taxonomy.metadata$assay             = "Smart-seq v4"  
taxonomy.metadata$anatomical_region = "middle temporal gyrus"
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
AIT.anndata = buildTaxonomy(title="Human_MTG_SMART_seq_04042025",
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
                            taxonomyDir = getwd(), ## This is where our taxonomy will be created
							addMapMyCells = TRUE, 
                            ##
                            subsample=1000000)  ## Using a huge number since we don't want to subsample

## Check whether the taxonomy file is valid (This happens within buildTaxonomy and is not strictly necessary)
print(paste("Is taxonomy valid?", AIT.anndata$uns$valid))


######################################################################################
## Upload to AWS

## The file created above is now uploaded to AWS and can be downloaded using the URL https://allen-cell-type-taxonomies.s3.us-west-2.amazonaws.com/Human_MTG_SMART_seq_04042025.h5ad


######################################################################################
## To access the taxonomy, 

# First download into your working directory from the link above

# Second, follow the steps from "========== Prepare our working environment ==========" above

# Third, load the taxonomy by uncommenting the line of code below
# AIT.anndata <- loadTaxonomy("Human_MTG_SMART_seq_04042025.h5ad")

