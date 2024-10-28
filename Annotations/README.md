# Annotation Schema

This schema supports the encoding of cell set evidence from multiple modality types include spatial transcriptomics, patch-seq, viral genetic tools as well as literature.

# Storage of evidence for cell set annotations

Building evdience from multiple modalities for individual cell sets requires some harmonization of metadata and features that describes a particular cellular phenotype. Harmonizing evidence across cell sets and studies will enable the construction of cell type cards and a comprehensive knowledge base.

This document describes the AIT Annotation schema, a type of contract, that all datasets can adhere to for enabling various products and tooling around taxonomies.

# Schema:

## scRNA-seq

Controlled schema for scRNA-seq description of cell sets (clustering, and higher) is already covered by [Cell Annotation Schema](https://github.com/cellannotation/cell-annotation-schema) which is best communicated by the [Cell Annotation Platform](https://docs.google.com/document/d/1CqL_t2CMcZF257rvjObDR_2WejsQ019rgUoX0YPJggM/edit?tab=t.0#bookmark=id.90zuscqjyay).

## Literature

Here we describe the schema for a single piece of literature evidence. It's expected that multiple literature references will be recorded for each cell set.

| Key | cellannotation_setname--literature |
| :---  | :--- |
| Value | A pre-defined cell set from the taxonomy |
| Type | Categorical with `str` categories. `Controlled`: Pre-defined cell set annotations including clustering and higher level groupings from a taxonomy. |

<br>

| Key | reference |
| :---  | :--- |
| Value | DOI link to associated paper |
| Type | String. `Free text` |

<br>

| Key | synonym |
| :---  | :--- |
| Value | Nomenclature used to describe the associated cell set in *reference* |
| Type | String. `Free text` |

<br>

| Key | evidence_type |
| :---  | :--- |
| Value | Description of the evidence or phenotype being presented in *reference*. for the associated *cellannotation_setname*. |
| Type | Categorical with `str` categories. `Controlled`: Predefined set of evidence categories. |

## Spatial

Spatial transcriptomics records the x,y,z coordinate of a cell in sectioned tissue. 

| Key | cellannotation_setname--spatial |
| :---  | :--- |
| Value | A pre-defined cell set from the taxonomy |
| Type | Categorical with `str` categories. `Controlled`: Pre-defined cell set annotations including clustering and higher level groupings from a taxonomy. |

<br>

| Key | dissected_roi |
| :---  | :--- |
| Annotator | Automatic |
| Value | Computed frequence for region of interest as determined from targted region during dissection. |
| Type | Categorical with `str` categories. `Controlled`: Pre-defined set of region names from anatomical atlas. |

<br>

| Key | merscope_roi |
| :---  | :--- |
| Annotator | Curator |
| Value | Expert description for the enriched region occupied by the cell set based on visual inspection of spatial transcriptomics. |
| Type | Categorical with `str` categories. `Controlled`: Pre-defined set of region names from anatomical atlas. |


<br>

| Key | common_coordinate_framework_roi |
| :---  | :--- |
| Annotator | Automatic |
| Value | Computed frequence for region of interest as determined from spatial transcriptomics data which has been aligned to a common coordinate framework. |
| Type | Categorical with `str` categories. `Controlled`: Pre-defined set of region names from anatomical atlas. |

<br>

| Key | expert_description--spatial |
| :---  | :--- |
| Value | Expert description of the cell set in terms of the cellular phenotypes captured by spatial profiling. |
| Type | String. `Free text` |

## Patch-Seq

Patch-Seq records a tri-modality view of the cell includuing transcriptomics, morphology and elecrophysiology. 

| Key | cellannotation_setname--patchseq |
| :---  | :--- |
| Value | A pre-defined cell set from the taxonomy |
| Type | Categorical with `str` categories. `Controlled`: Pre-defined cell set annotations including clustering and higher level groupings from a taxonomy. |

<br>

| Key | expert_description--patchseq |
| :---  | :--- |
| Value | Expert description of the cell set in terms of the cellular phenotypes captured by patch-seq. |
| Type | String. `Free text` |

## Example: Cell type card stitching. 

STR FS PTHLH ST18 GABA corresponds to PVALB+, PTHLH+ (ref), { synonym } (ref).   

Supporting evidence:
Location:  …  This is supported by its location in the striatum (references) and merfish data, its morphology {some description} (reference) 
and patch-seq data and that it is fast spiking (references) - summary of e-phys data supporting fast spiking.
