# Allen Institute Taxonomy schema

We have developed a compartmenatlized schema for storing all required aspects of a taxonomy. The fields in the AIT schema are associated to a broad category term (described below) which form a piece of the whole AIT file format. 

*(Note: A pervious version of this standard is available **[as a Google Doc](https://docs.google.com/document/d/1nj6LHUPoo3JnNwZ7PTdniT9pBPsoJr1B/edit?usp=sharing&ouid=113573359044104089630&rtpof=true&sd=true)**).*

### Schema category terms

Described here are the broad categories that all fields are associated with.

* **[Data](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#data)**: Includes anything critical for understanding the cell by gene matrix and to link it with other components.  Includes data (raw and processed), gene information, and cell identifiers.  *Note that for the purposes of this schema, we are excluding raw data (fastq, bam files, etc.) from consideration and are starting from the count matrix.*
  
* **[Assigned metadata](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#assigned-metadata)**: Includes cell-level metadata that is assigned at some point in the process between when a cell goes from the donor to a value in the data, and (in theory) can be ENTIRELY captured by values in Allen Institute, BICAN, or related standardized pipelines.  It includes fields that describe: donor metadata, experimental protocols, dissection information, RNA QC metrics, and sequencing metadata.
  
* **[Calculated metadata](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#calculated-metadata)**: Includes any cell-level or cluster-level metadata that can be calculated explicitly from the **Data** and **Assigned Metadata** without the need for human intervention. It includes fields that describe: # reads detected/cell, # UMI/cell, fraction of cells per cluster derived from each anatomic dissections, expressed neurotransmitter genes (quantitatively defined), standard quality control metrics (e.g., doublet score) per cluster.
  
* **[Annotations](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#annotations)**: Includes fields related to the annotation of clusters or groups of clusters (collectively called "cell sets").  It includes fields that describe cluster levels, cluster relationships, canonical marker genes, links to existing ontologies (e.g., CL, UBERON), expert annotations, and dendrograms.
  
* **[Analysis](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#analysis)**: Includes fields which are required for specific analysis for example: latent spaces (e.g., UMAP), cluster level gene summaries (e.g., cluster means, proportions), and variable genes.
  
* **[Tooling](https://github.com/AllenInstitute/scrattch.taxonomy/blob/KL_div/schema/aligned_schema.md#tooling)**: Includes fields required for specific tools (e.g., cellxgene, TDT, CAS, CAP, and cell type annotation) that are not strictly part of the taxonomy but are required to inter-operate between various tools.

<!-- Here is a graphical representation of these terms in the context of data, metadata, and taxonomies:
![image](https://github.com/AllenInstitute/scrattch.taxonomy/assets/25486679/eaf6b3d3-0b5f-49fc-9a49-2b7168605964) -->

## Schema

Within each broad categorical term, fields are ordered by their location in the anndata object: X, raw, obsm, obs, var, uns.

<!-- ![Schematic](https://github.com/AllenInstitute/scrattch.taxonomy/blob/main/schema/AIT_anndata_schema.png) -->

The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "NOT RECOMMENDED" "MAY", and "OPTIONAL" in this document are to be interpreted as described in [BCP 14](https://tools.ietf.org/html/bcp14), [RFC2119](https://www.rfc-editor.org/rfc/rfc2119.txt), and [RFC8174](https://www.rfc-editor.org/rfc/rfc8174.txt) when, and only when, they appear in all capitals, as shown here.

## `X` 

| Key | X |
| :-- | :-- |
| Annotator | Curator |
| Value | `X` component contains normalized expression data (cell x gene) in [scipy.sparse.csr_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html) matrix format. |
| Type| `Index`, `str` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

## `raw`

The `raw` component contains the unfiltered anndata object containing a count matrix in `raw.X`.

| Key | raw.X |
| :-- | :-- |
| Annotator | Curator |
| Value | `raw.X` component contains a count matrix in [scipy.sparse.csr_matrix](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html) matrix format.|
| Type| `Index`, `str` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |


## `obs`

`obs` is a [pandas.Dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)

The `obs` component contains **cell-level metadata** summarized at the cell level. 

#### index of pandas.Dataframe
| Key | index of pandas.Dataframe |
| :-- | :-- |
| Annotator | Curator |
| Value | Unique identifier corresponding to each individual cell. |
| Type| `Index`, `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### cluster_id
| Key | cluster_id |
| :-- | :-- |
| Annotator | Curator |
| Value | Identifier for cell set computed from a clustering algorithm. |
| Type| `str` |
| Required | MUST |
| Tags | Annotations |

<br>

#### [cellannotation_setname]
| Key | [cellannotation_setname] |
| :-- | :-- |
| Annotator | Curator |
| Value | Column name in `obs` is the string [cellannotation_setname] and the values are the strings describing an annotation level of the taxonomy.
| Type| `Categorical` |
| Required | RECOMMENDED |
| Tags | Annotations |

Examples: `Neuronal`, `Inhibitory`, `LHX6 (MGE)`, `PVALB`, `Inh L5-6 PVALB LGR5` from [Hodge et al. 2019](https://www.nature.com/articles/s41586-019-1506-7)

<br>

#### cell_type_ontology_term_id
| Key | cell_type_ontology_term_id |
| :-- | :-- |
| Annotator | Curator |
| Value | This MUST be a CL term. If no appropriate high-level term can be found or the cell type is unknown, then it is STRONGLY RECOMMENDED to use "CL:0000003" for native cell.
| Type| `Categorical` |
| Required | RECOMMENDED |
| Tags | Annotations |

<br>

#### load_id
| Key | load_id |
| :-- | :-- |
| Annotator | Curator |
| Value | Identifier for the sequencing library for which molecular measurements from a specific set of cells is derived. |
| Type| `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### donor_id
| Key | donor_id |
| :-- | :-- |
| Annotator | Curator |
| Value | Identifier for the unique individual, ideal from the specimen portal (or other upstream source). This is called `donor_label` in the **BKP**. Should converge on a standard term. More than one identifier may be needed, but ideally for the analysis only a single one is retained and stored here. |
| Type| `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### assay
| Key | assay |
| :-- | :-- |
| Annotator | Curator |
| Value | Human-readable sequencing modality which should have a corresponding EFO ontology term. e.g., 'Smart-seq2'corresponds to 'EFO:0008931', '10x 3' v3'corresponds to 'EFO:0009922'.
| Type| `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### organism
| Key | organism |
| :-- | :-- |
| Annotator | Curator |
| Value | Species from which cells were collected. This MUST be the human-readable name assigned to the value of organism_ontology_term_id
| Type| `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### organism_ontology_term_id
| Key | organism_ontology_term_id |
| :-- | :-- |
| Annotator | Computed |
| Value | NCBITaxon identifier which MUST be a child of NCBITaxon:33208 for Metazoa. |
| Type| `str` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

<br>

#### donor_age
| Key | donor_age |
| :-- | :-- |
| Annotator | Curator |
| Value | Currently a free text field for defining the age of the donor. In **CELLxGENE** this is recorded in `development_stage_ontology_term_id` and is HsapDv if human, MmusDv if mouse.  I'm not sure what this means, but more generally, we should align with BICAN on how to deal with this value. |
| Type| `Categorical` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### anatomical_region
| Key | anatomical_region |
| :-- | :-- |
| Annotator | Curator |
| Value | Human readable name assigned to the value of `anatomical_region_ontology_term_id` |
| Type| `Categorical` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### anatomical_region_ontology_term_id
| Key | anatomical_region_ontology_term_id |
| :-- | :-- |
| Annotator | Computed |
| Value | UBERON terms for the `anatomical_region` field that we have (e.g., 'brain': 'UBERON_0000955'). |
| Type| `Categorical` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

<br>

#### self_reported_sex
| Key | self_reported_sex |
| :-- | :-- |
| Annotator | Curator |
| Value | Placeholder for donor reported sex. Called `sex_ontology_term_id` (e.g., PATO:0000384/383 for male/female) in **CELLxGENE** and called "donor_sex" in BKP. We should align on a single term. |
| Type| `Categorical` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### self_reported_ethnicity_ontology_term_id
| Key | self_reported_ethnicity_ontology_term_id |
| :-- | :-- |
| Annotator | Computed |
| Value | This MUST be a child of PATO:0001894 for phenotypic sex or "unknown" if unavailable. |
| Type| `Categorical` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

<br>

#### disease_ontology_term_id
| Key | disease_ontology_term_id |
| :-- | :-- |
| Annotator | Curator |
| Value | This MUST be a MONDO term or "PATO:0000461" for normal or healthy. |
| Type| `Categorical` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

<br>

#### suspension_type
| Key | suspension_type |
| :-- | :-- |
| Annotator | Curator |
| Value | Either "cell", "nucleus", or "na". |
| Type| `Categorical` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### is_primary_data
| Key | is_primary_data |
| :-- | :-- |
| Annotator | Curator |
| Value | This MUST be True if this is the canonical instance of this cellular observation and False if not. This is commonly False for meta-analyses reusing data or for secondary views of data. |
| Type| `bool` |
| Required | MUST |
| Tags | Assigned metadata |

## `var`

`var` is a [pandas.Dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)

The `var` component contains gene level information.

#### index of pandas.Dataframe
| Key | index of pandas.Dataframe |
| :-- | :-- |
| Annotator | Curator |
| Value | If the feature is a gene then this MUST be an human readable SYMBOL term. The index of the pandas.DataFrame MUST contain unique identifiers for features. If present, the index of raw.var MUST be identical to the index of var. |
| Type| `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### ensembl_id
| Key | ensembl_id |
| :-- | :-- |
| Annotator | Curator |
| Value | If the feature is a gene then this MUST be a gene ID from ensembl. Each index of the pandas.DataFrame MUST map to a unique emsembl_id identifiers for features. If present, the raw.var.ensembl_id MUST be identical to the var.ensembl_id. |
| Type| `str` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

<br>

#### highly_variable_genes
| Key | highly_variable_genes |
| :-- | :-- |
| Annotator | Curator |
| Value | A logical vector indicating which genes are highly variable. |
| Type| `bool` |
| Required | RECOMMENDED |
| Tags | Analysis |

<br>

#### marker_genes
| Key | marker_genes_[set_name] |
| :-- | :-- |
| Annotator | Curator |
| Value | A logical vector indicating which genes are markers. Multiple marker gene sets can be specified. |
| Type| `bool` |
| Required | RECOMMENDED |
| Tags | Analysis |

<br>

## `uns`

The `uns` component contains more general information and fields with formatting incompatible with the above components.

#### title
| Key | title |
| :-- | :-- |
| Annotator | Curator |
| Value | This text describes and differentiates the dataset from other datasets in the same collection. It is STRONGLY RECOMMENDED that each dataset title in a collection is unique and does not depend on other metadata such as a different assay to disambiguate it from other datasets in the collection. |
| Type| `str` |
| Required | MUST |
| Tags | Tooling |

<br>

#### dataset_purl
| Key | dataset_purl |
| :-- | :-- |
| Annotator | Curator |
| Value | Link to molelcular data (cell x gene) if not present in X or raw.X. This can be an AWS S3 bucket or other permanent URL for the taxonomy expression data. |
| Type| `str` |
| Required | RECOMMENDED |
| Tags | Data |

<br>

#### batch_condition
| Key | batch_condition |
| :-- | :-- |
| Annotator | Curator |
| Value | Together, these keys define the batches that a normalization or integration algorithm should be aware of. Values MUST refer to cell metadata keys in `obs`. |
| Type| `list[str]` |
| Required | RECOMMENDED |
| Tags | Tooling |

<br>

#### reference_genome
| Key | reference_genome |
| :-- | :-- |
| Annotator | Curator |
| Value | Reference genome used to align molecular measurements. |
| Type| `str` |
| Required | RECOMMENDED |
| Tags | Assigned metadata |

<br>

#### gene_annotation_version
| Key | gene_annotation_version |
| :-- | :-- |
| Annotator | Curator |
| Value | Genome annotation version used during alignment. e.g. .gtf or .gff file.  |
| Type| `str` |
| Required | MUST |
| Tags | Assigned metadata |

<br>

#### dend
| Key | dend |
| :-- | :-- |
| Annotator | Curator |
| Value | A json formatted dendrogram encoding the taxonomy heirarchy. Either computed or derived from cluster groupings.  |
| Type| `json` |
| Required | RECOMMENDED |
| Tags | Annotations |

<br>

#### heirarchy
| Key | heirarchy |
| :-- | :-- |
| Annotator | Curator |
| Value | An ordering of `cluster_id` and higher level groupings from `[cellannotation_setname]`. E.g. ["Class", "Subclass", "cluser_id"]  |
| Type| `list[str]` |
| Required | MUST |
| Tags | Annotations |

<br>

#### mode
| Key | mode |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | Indicator of which child of the parent taxonomy to utilize. Mode determines cells to remove based on `filter` as well as switching to relevant analysis components of the `uns` related to child taxonomy specific analysis tooling. |
| Type| `str` |
| Required | MUST |
| Tags | Tooling |

<br>

#### filter
| Key | filter |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | Indicator of which cells to use for a given child taxonomy saved as a list of booleans for each cell. TRUE indicates a cell should be removed and FALSE indicates the cell should not be removed. Each entree in this list is named for the relevant "mode" and has TRUE/FALSE calls indicating whether a cell is filtered out. e.g., the "standard" taxonony is all FALSE. |
| Type| `list[[mode]][bool]` |
| Required | MUST |
| Tags | Tooling |

<br>

#### qualty_control_markers
| Key | qualty_control_markers |
| :-- | :-- |
| Annotator | Curator |
| Value | Marker gene expression in on-target and off-target cell populations, useful for patchseq analysis.  Also includes information about KL divergence calculations and associated QC calls. Defined by buildPatchseqTaxonomy. |
| Type| `data.frame` |
| Required | RECOMMENDED |
| Tags | Analysis |

<br>

#### cluster_info
| Key | cluster_info |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | A data.frame of cluster information. |
| Type| `data.frame` |
| Required | MUST |
| Tags | Annotations |

<br>

#### cluster_id_median_expr
| Key | cluster_id_median_expr |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | Marker gene expression in on-target and off-target cell populations, useful for patchseq analysis.  Also includes information about KL divergence calculations and associated QC calls. Defined by buildPatchseqTaxonomy. |
| Type| `numpy.ndarray` |
| Required | MUST |
| Tags | Annotations |

<br>

#### default_embedding
| Key | default_embedding |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | The value MUST match a key to an embedding in `obsm` for the embedding to display by default. |
| Type| `str` |
| Required | RECOMMENDED |
| Tags | Tooling |

<br>

#### schema_version
| Key | schema_version |
| :-- | :-- |
| Annotator | Computed |
| Value | Allen Institute Taxonomy schema version. e.g. "1.0.0" |
| Type| `str` |
| Required | MUST |
| Tags | Tooling |

<br>

#### cellannotation_schema
| Key | cell_annotation_schema |
| :-- | :-- |
| Annotator | Computed |
| Value | A json storing the entire cell annotation schema (CAS) information. |
| Type| `json` |
| Required | RECOMMENDED |
| Tags | Tooling |

<br>

#### cellannotation_schema_version
| Key | cellannotation_schema_version |
| :-- | :-- |
| Annotator | Computed |
| Value | A version id for the CAS schema. |
| Type| `str` |
| Required | RECOMMENDED |
| Tags | Tooling |

<br>

* `cell_annotation_schema`: extended `calculated metadata` about annotations and labelsets stores in `uns` as in [CAS - BICAN extension](https://github.com/cellannotation/cell-annotation-schema/blob/main/build/BICAN_schema.md) format under `labelsets`.  

## `obsm` (Embeddings)

The `obsm` component contains all dimensionality reductions of the taxonomy (cell x dim). To display a dataset Curators MUST annotate one or more embeddings of at least two-dimensions (e.g. tSNE, UMAP, PCA, spatial coordinates) as numpy.ndarrays in obsm.

#### X_[embedding]
| Key | X_[embedding] |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | An n-dimensional embedding (cell x dim) of the high dimensional expression data. |
| Type| `numpy.ndarray` |
| Required | MUST |
| Tags | Analysis |

<br>
