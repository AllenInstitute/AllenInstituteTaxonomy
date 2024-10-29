# Allen Institute Taxonomy (AIT)

To distribute Allen Institute Taxonomies (AIT) we define an [`anndata`](https://anndata.readthedocs.io/en/latest/index.html) .h5ad file which encapsulates the essential components of a taxonomy required for downstream analysis such as [cell type mapping](https://github.com/AllenInstitute/scrattch-mapping/tree/main) with a [formalized schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/schema/README.md).

## Overview

*(Note: A pervious version of this standard is available **[as a Google Doc](https://docs.google.com/document/d/1nj6LHUPoo3JnNwZ7PTdniT9pBPsoJr1B/edit?usp=sharing&ouid=113573359044104089630&rtpof=true&sd=true)**).*

Several competing schema have been created for packaging of taxonomies, data sets, and associated metadata and annotations.  This document aims to align three such schema and propose a way of integrating them into the Allen Institute Taxonomies (AIT) .h5ad file format presented as part of this GitHub repository. The three standards are:

1. **AIT** (described herein)
2. **[Cell Annotation Schema](https://github.com/cellannotation/cell-annotation-schema/) (CAS)**: this schema is becoming more widely-used in the cell typing field as a whole because it is largely compatible with [the CZ CELLxGENE schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md). It is also compabible with [Cell Annotation Platform](https://celltype.info/) (CAP) and with [Taxonomy Development Tools](https://brain-bican.github.io/taxonomy-development-tools/) (TDT). CAS has both a general schema and a BICAN-associated schema, both of which are considered herein.  CAS can be embedded in the header (`uns`) of an AIT/Scraatch.taxonomy file, where it functions as a store of extended information about an annotation, including ontology term mappings, evidence for annotation (from annotation transfer and marker expression).
3. **Brain Knowledge Platform (BKP)**: this schema isn't publicly laid out anywhere that I can find, but this is the data model used for [Jupyter Notebooks](https://alleninstitute.github.io/abc_atlas_access/intro.html) associated with the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas).  More generally, any data sets to be included in ABC Atlas, [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells), or other related BKP resources will eventually need to conform to this format.

## Cell type taxonomy organization

One major challenge in creating a cell type taxonomy schema is in definition of terms such as "taxonomy",  "dataset", "annotation", "metadata", and "data".  It is becoming increasingly important to separate out the data from the other components, and compartmentalize all components to avoid the need to download, open, or upload huge and unweildy files. 

![Taxonomy_overview](https://github.com/AllenInstitute/scrattch.taxonomy/assets/25486679/9d36e6bc-db14-4d73-8011-23026756ec08)

That said, it is still important for many use cases to have an option of including all of the information listed above in a single h5ad file for use with CELLxGENE, [scrattch.mapping](https://github.com/AllenInstitute/scrattch.mapping), analysis tools, and for ease of sharing in a single file format. 
