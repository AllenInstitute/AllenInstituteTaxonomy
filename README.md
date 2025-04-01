# Allen Institute Taxonomy (AIT)

To distribute Allen Institute Taxonomies (AIT) we define an **[`anndata` .h5ad file](https://anndata.readthedocs.io/en/latest/index.html)** which encapsulates the essential components of a taxonomy required for downstream analysis with a **[formalized schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy/blob/main/schema/README.md)**. *For information on how to build and work with AIT files, see the companion **[scrattch R libraries](https://github.com/AllenInstitute/scrattch)***.

## Overview

One major challenge in creating a cell type taxonomy schema is in definition of terms such as "taxonomy",  "dataset", "annotation", "metadata", and "data".  It is becoming increasingly important to separate out the data from the other components, and compartmentalize all components to avoid the need to download, open, or upload huge and unweildy files.  *AIT addresses this challenge by extending and modifying the popular **[CELLxGENE schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.2.0/schema.md)** to better align with BICAN and Allen Institute needs.*

![Taxonomy_overview](https://github.com/AllenInstitute/scrattch.taxonomy/assets/25486679/9d36e6bc-db14-4d73-8011-23026756ec08)

(Brief description of AIT and it's difference from CELLxGENE to be entered here.)

## Related efforts

AIT is being developed alongside three complementary efforts for packaging of taxonomies, data sets, and associated metadata and annotations.

1. **[Cell Annotation Platform (CAP)](https://celltype.info/)**: CAP 'is a centralized, community-driven platform for the creation, exploration, and storage of cell annotations for single-cell RNA-sequencing (scRNA-seq) datasets.' The Allen Institute and BICAN are partnering with CAP for annotation of brain (including basal ganglia) and spinal cord taxonomies.
2. **[Cell Annotation Schema](https://github.com/cellannotation/cell-annotation-schema/) (CAS)**: Compatible with [Cell Annotation Platform] (CAP) and with [Taxonomy Development Tools (TdT)](https://brain-bican.github.io/taxonomy-development-tools/), functions as a store of extended information about cell sets, including ontology term mappings and evidence for annotation (from annotation transfer and marker expression). CAS complements other cell-centric and occasionally cluster-centric schema more commonly used. CAS has both a general schema and a BICAN-associated schema, and can be embedded in the header (`uns`) of an AIT file.
3. **[Brain Knowledge Platform (BKP)]**: While not publicly laid out anywhere that I can find, the BKP schema is the data model used for [Jupyter Notebooks](https://alleninstitute.github.io/abc_atlas_access/intro.html) associated with the [Allen Brain Cell (ABC) Atlas](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas) and will eventually power all novel content hosted on [Allen Brain Map](https://portal.brain-map.org/atlases-and-data).  Currently, any data sets to be included in ABC Atlas or [MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells) must be ingested into BKP.
