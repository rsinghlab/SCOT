---
layout: template
permalink: /data/
---

## Data
Here, you can download the data used in the examples in the paper describing SCOT:

#### Simulated Datasets:

|                                 | Domain 1 | Domain 2 | Notes |
| --------------------------------|----------|----------|----------
| **Simulation 1**: Bifurcating Tree  |[300 x 1000](data/s1_mapped1.txt)|[300 x 2000](data/s1_mapped2.txt)| Originally from [here](https://noble.gs.washington.edu/proj/mmd-ma/)|
| **Simulation 2**: Swiss Roll        |[300 x 1000](data/s2_mapped1.txt)|[300 x 2000](data/s2_mapped2.txt)| Originally from [here](https://noble.gs.washington.edu/proj/mmd-ma/)|
| **Simulation 3**: Circular Frustum  |[300 x 1000](data/s3_mapped1.txt)|[300 x 2000](data/s3_mapped2.txt)| Originally from [here](https://noble.gs.washington.edu/proj/mmd-ma/)|
| **Simulation 4**: Synthetic RNA-seq |[5000 x 50](data/s4_splatterX.txt)|[5000 x 500](data/s4_splattery.txt)|Generated using [Splatter](https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html)|

#### Real-world Sequencing Datasets:
**<ins>Note</ins>** that the files in "Domain 1" and "Domain 2" columns of the real sequencing datasets contain data pre-processed according to their original publications (linked in Notes), so they are dimensionality reduced. To get access to the original raw datasets, follow the "Raw data" links.

|                                       | Domain 1 | Domain 2 | Notes |
| --------------------------------------|----------|----------|-------|
| **SNAREseq** Cell Line Mixture            |[1047 x 19 (chromatin accessibility)](data/snare_chromatin.txt)|[1047 x 10 (gene expression)](data/snare_rna.txt)|Original publication [here](https://www.nature.com/articles/s41587-019-0290-0). Raw data [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074) |
| **scGEM** Dataset                         |[177 x 34 (gene expression)](data/scGEM_expression.txt)|[177 x 27 (DNA methylation)](data/scGEM_methylation.txt)|[Original publication [here](https://pubmed.ncbi.nlm.nih.gov/27525975/). Raw data [here](https://www.nature.com/articles/nmeth.3961#Sec11) |

##### Don't hesitate to [contact us](rsinghlab.github.io/SCOT/contact) if you have any questions about these datasets.