
![](assets/SCOT_logo.png)plo

**Note:** This GitHub page is a work in progress. It will contain information on the SCOT methodology, as well as demonstrated examples in various application cases.  

### Overview
SCOT is an unsupervised algorithm for cell-to-cell alignment of single-cell multi-omic data. <br>
For many single-cell technologies, it is difficult to apply them simultaneously on the same single-cell. As a result, to get multiple views of a genome at the single-cell resolution, we need to integrate data from different -omic sequencing datasets. When these single-cell sequencing methods are applied to the same cell population or different populations that are expected to share some underlying biological manifold (e.g. common cell types), we expect there to be some cell-to-cell correspondence/alignment to recovered. SCOT performs alignment by using Gromov-Wasserstein optimal transport.  
<br> 
SCOT works in three steps:  
**1.** It first checks pairwise correlations between samples in each sequencing dataset and constructs k nearest-neighbor graphs for each to computing intra-domain distances while accounting for local geometry.  
**2.** It then optimizes the Gromov-Wasserstein optimal transport formulation to look for a probabilistic sample-wise correspondence matrix between the input datasets.  
**3.** Using this coupling matrix, it projects one dataset onto the other one using barycentric projection.  

### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/rsinghlab/SCOT/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.
