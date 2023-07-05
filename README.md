# AllenInstituteTaxonomy

To distribute Allen Institute Taxonomies (AIT) we define an [`anndata`](https://anndata.readthedocs.io/en/latest/index.html) .h5ad file which encapsulates the essential components of a taxonomy required for downstream analysis such as [cell type mapping](https://github.com/AllenInstitute/scrattch-mapping/tree/main).

### Schematic

![Schematic](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/AIT_anndata_schema.png)

### Anndata components

#### X

The `X` component contains logCPM normalized expression data (cell x gene).

#### layers

The `layers` component contains the count matrix (cell x gene). (**Optional**)

* `counts`: The count matrix from which `X` was derived.

#### obsm

The `obsm` component contains all dimensionality reductions of the taxonomy (cell x dim). 

* `umap`: 2-dimensional representation of cells.

#### obs

The `obs` component contains cell level metadata including cell type annotation, brain region and species. Standardized nomenclature for cell metadata to be decided.

* `cluster`
* `brain_region`
* `species`
* `age`
* ...

#### var

The `var` component contains gene level metadata including which genes are highly variable or marker genes. Additional information may enclude Ensemble IDs.

* `gene`: A vector of gene symbols
* `ensemble_id`: A vector of corresponding ensemble ids for each gene symbol. (**Optional**)
* `highly_variable_genes`: A logical vector (T/F) indicating which genes are highly variable.
* `marker_genes`: A logical vector (T/F) indicating which genes are markers used to build dendrogram. (**Optional**)


#### uns

The `uns` component contains taxonomy associated files useful for reproducing analysis or mapping against the taxonomy.

* `dend`: A file path to the dendrogram .RData object. Will be replaced by language agnostic encoding of dendrogram.
* `filter`: A logical vector indicating which cells to filter out for a specific taxonomy version.
* `QC_markers`: Marker gene expression in on-target and off-target cell populations, useful for patchseq analysis.
* `mode`: Taxonomy mode that determins which `filter` to use. See scrattch.mapping documentation.
* `clustersUse`: A vector of cluster names to use for taxonomy.
* `clusterInfo`: A data.frame of cluster information.
* `taxonomyName`: Taxonomy name and AIT identifier.
* `taxonomyDir`: Directory where taxonomy related files can be located. (**Optional**)
