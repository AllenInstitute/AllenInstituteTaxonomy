## What are taxonomy ‘modes’
AIT files include the concept of a taxonomy mode. A taxonomy mode represents a specific **subset** of cells included in the AIT file, along with all of the associated variables required for performing mapping and patch-seq analysis. Upon creation, an AIT file will default to the ‘standard’ mode (sometimes called the ‘parent taxonomy’), which includes all or a subsampled selection of cells from every `cluster_id` in the taxonomy.  Additional taxonomy modes (or ‘child taxonomies’) represent specific subsets of cells in the ‘standard’ mode and can represent a combination of filtered cell types and/or additional subsampling.  

## Why build taxonomy modes
Taxonomy modes provide flexibility to perform many different analyses using the same base taxonomy, rather than needing to save multiple copies of the same taxomy.  Common use cases for taxonomy modes include: 

* **Regional taxonomies**: for example, the mouse basal ganglia taxonomy represents all cells from a subset of clusters from the whole mouse brain taxonomy from Yao et al 2023.
* **Patch-seq analysis**: In both human and mouse cortex, we subset taxonomies to include only neuronal cell types since patch-seq data is often noisy and cells will erroneously map to glial types if all clusters are included
* **Quality control**: When all cells are included by default, modes can be created that omit all clusters and cells failing QC
* **Algorithm testing**: To quickly testing algorithms, it is often useful to define a mode with just a handful of clusters for testing.

## What mode-specific variables are stored in an AIT file
Mode-specific variables can either be user-curated (U) or computed (C) and fall in a few general categories:
* **Data organization (U)**: Different embeddings (e.g., UMAP, t-SNE, scVI dimensions) and dendrograms can be stored for different subsets of the taxonomy.
* **Genes of importance (U/C)**: Various parameters for storing high variance genes, marker genes, and genes used for mapping algorithms can be stored separately for each mode.
* **Mapping statistics (C)**:  Mapping algorithms require pre-calculated statistics and variables that can be stored separately for each mode.
* **Patch-seq analysis metrics (C)**: Specific patch-seq QC metrics can be saved for different taxonomy modes (e.g., for mapping to glutamatergic vs. GABAergic types).

## Schema
All variables are tagged as “Curator” and/or “Computed” for compatibility with the general [AIT schema documentation](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/schema/README.md). Rather than tagging terms are REQUIRED or RECOMMENDED, information is provided about the specific use case for each variable.

## `var`

`var` is a [pandas.Dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html)

The `var` component contains gene level information.  These variables either match or supplement variables from the [`var` component of the general schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema#var).

#### highly_variable_genes_[mode]
| Key | highly_variable_genes |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | A logical vector indicating which genes are highly variable. Multiple highly variable gene sets can be specified. By default a set of genes based on binary scores can be automatically calculated and stored in the `highly_variable_genes_[mode]` for use with mapping algorithms. |
| Type| `bool` |
| Used for | `corr.map` and `seurat.map` in `scrattch.mapping`; defining embeddings in `scrattch.taxonomy` |

<br>

#### marker_genes_[mode]
| Key | marker_genes_[set_name] |
| :-- | :-- |
| Annotator | Curator |
| Value | A logical vector indicating which genes are markers. Multiple marker gene sets can be specified.  |
| Type| `bool` |
| Used for | (Not currently used in `scrattch` packages) |

<br>

### `obsm` (Embeddings)

The `obsm` component contains all dimensionality reductions of the taxonomy (cell x dim). To display a dataset Curators MUST annotate one or more embeddings of at least two-dimensions (e.g. tSNE, UMAP, PCA, spatial coordinates) as numpy.ndarrays in obsm.  These variables either match or supplement variables from the [`obsm` component of the general schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema#obsm-embeddings).

#### X_[embedding]
| Key | X_default_[mode] |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | An n-dimensional embedding (cell x dim) of the high dimensional expression data. Can be provided or computed based off `highly_variable_genes_[mode]` |
| Type| `numpy.ndarray` |
| Used for | (Not currently used in `scrattch` packages, but planned future usage for visualization and constellation diagram creation.)  |

### `uns`

The `uns` component contains more general information and fields with formatting incompatible with the above components. These variables either match or supplement variables from the [`uns` component of the general schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema#uns).

#### mode
| Key | mode |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | Indicator of which child of the parent taxonomy to utilize. Mode determines cells to remove based on `filter` as well as switching to relevant analysis components of the `uns` related to child taxonomy specific analysis tooling. |
| Type| `str` |
| Used for | All `scrattch` functions involving modes. |

<br>

#### filter
| Key | filter |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | Indicator of which cells to use for a given taxonomy mode saved as a list of booleans for each cell. TRUE indicates a cell should be removed and FALSE indicates the cell should not be removed. Each entree in this list is named for the relevant "mode" and has TRUE/FALSE calls indicating whether a cell is filtered out. e.g., the "standard" taxonony is all FALSE, modulo subsampling. |
| Type| `list[[mode]][bool]` |
| Used for | Preparing taxonomy modes for associated downstream analyses.  |

<br>

#### dend
| Key | dend |
| :-- | :-- |
| Annotator | Curator |
| Value | A json formatted dendrogram encoding the taxonomy hierarchy (see notes). Either computed or derived from cluster groupings.  |
| Type| `list[[mode]][json]` |
| Used for | `tree.map` mapping in `scrattch.mapping` and defining quality control metrics and cluster membership in `scrattch.patchseq` |
| Notes | `dend` must include four components (derived in R from the `hclust` function in the `stats` library) and can optionally include any other components. <br> **merge**: A [[list]] that describes the sequential merging [steps] of clusters at each step of the hierarchical clustering process into individual tree "nodes". Each [row] of the merge list represents a merging step. Negative values steps indicate the indices of individual observations (e.g., clusters), while positive values indicate nodes. <br> **node_heights**: A vector containing the heights (or distances) at which the clusters were merged. Typically, these heights represent the dissimilarity between the merged clusters. <br> **labels**: Cluster names corresponding to every leaf node in order (after reordering by the "order" value). More generally this is a vector of labels for the observations that were clustered.   <br> **order**: This is a vector that specifies the order in which the observations should be arranged to produce a dendrogram without crossing branches. If labels are ordered from left to right on the tree, then order would be [1,2,3,...,N] <br> **markers(?)** Lists of marker genes for each node of the tree (*NOTE: THIS NEEDS TO BE DESCRIBED BETTER*). |

<br>

#### clusterStatsColumns
| Key | clusterStatsColumns |
| :-- | :-- |
| Annotator | Computed |
| Value | Character vector of `cluster_id` values included in [mode].  |
| Type| `list[[mode]][str]` |
| Used for | `corr.map`, `mapmycells.hierarchical.map`, and `mapmycells.flat.map` mapping in `scrattch.mapping`. |

<br>

#### default_embedding
| Key | default_embedding_[mode] |
| :-- | :-- |
| Annotator | Curator/Computed |
| Value | Each value MUST match a key to an embedding in `obsm` for the embedding to display by default for [mode]. *Note: this is not yet implemented as of 18 Apr 2025.* |
| Type| `list[[mode]][str]` |
| Used for | (Not currently used in `scrattch` packages, but planned future usage for visualization and constellation diagram creation.) |

<br>

#### mapmycells
| Key | mapmycells |
| :-- | :-- |
| Annotator | Computed |
| Value | Sets of [mode]-specific precomputed statistics (`mapmycells[[mode]][["precomp_stats"]]`) and associated marker genes (`mapmycells[[mode]][["query_markers"]]`) for use with [`call_type_mapper`](https://github.com/AllenInstitute/cell_type_mapper), the python package underlying hierarchical mapping algorithms in [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells). (*NOTE: THIS NEEDS TO BE DESCRIBED BETTER*).  |
| Type| `list[[mode]][list]` |
| Used for | `mapmycells.hierarchical.map` and `mapmycells.flat.map` mapping in `scrattch.mapping`. |

<br>

#### memb
| Key | memb |
| :-- | :-- |
| Annotator | Computed |
| Value | Two metrics for determining the confusion between cluster `memb`ership in a taxonomy [mode]. (See notes below) |
| Type| `list[[mode]][list]` |
| Used for | Defining quality control metrics and cluster membership in `scrattch.patchseq`, and a (not currently implemented) method for building constellation diagrams. |
| Notes | **`memb[[mode]][['memb.ref']]`**: matrix indicating how much confusion there is the mapping between each cell all of the nodes in the tree (including all cell types) when comparing clustering and mapping results with various subsamplings of the data. <br> **`memb[[mode]][['**map.df.ref**']]`**: Result of tree mapping for each cell in the reference against the clustering tree, including various statistics and marker gene evidence.  This is the same output that comes from `tree.map` when include AIT.anndata$X as input |

<br>

#### QC_markers_[mode]
| Key | memb |
| :-- | :-- |
| Annotator | Computed |
| Value | A list of several variables required for applying [`patchseqQC`](https://github.com/PavlidisLab/patchSeqQC/tree/master) to a set of query data for QCing patchseq data. Can be defined separately for each [mode]. (See notes below) |
| Type| `list[[mode]][list]` |
| Used for | Defining quality control metrics in `scrattch.patchseq`. |
| Notes | **`QC_markers_[[mode]][['markers']]`**: Output from the `defineClassMarkers` function. <br> **`QC_markers_[[mode]][['allMarkers']]`**: A character vector of all genes included in 'markers' above. <br> **`QC_markers_[[mode]][['countsQC']]`**: Count matrix of reference data set including the subset of genes and cells needed for `patchseqQC`. Required because X and raw.X are optional. <br> **`QC_markers_[[mode]][['cpmQC']]`**: Count per million of `countsQC` above. Required because all genes not saved in `countsQC`. <br> **`QC_markers_[[mode]][['classBr']]`**: Categorical vector (e.g., factor) of class-level calls, which mixes two [cellannotation_setname] levels of the hierarchy for on vs. off-target types. <br> **`QC_markers_[[mode]][['subclassF']]`**: Categorical vector (e.g., factor) of subclass-level calls. <br>  |
<br>
