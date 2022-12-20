# TCellPack v2.0
TCellPack is an R package used to visualize characteristics of T cell repertoires. It is possible to visualize T cell clonotype diversity and abundance in relation to its antigen specificity and phenotype in a single plot.

### About
T cell packs are circle packing plots representing T cell repertoires. They are able to visualize output data from [GLIPH](https://github.com/immunoengineer/gliph), [GLIPH2](http://50.255.35.37:8080), [Adaptive ImmunoSeq](https://www.adaptivebiotech.com/products-services/immunoseq), and [10X Chromium single cell immune profiling](https://www.10xgenomics.com/solutions/vdj). The GLIPH algorithm clusters T cell receptor CDR3 sequences by their sequence similarity which implies a high probability of sharing the same antigen specificity (i.e. shared CDR3 motifis may be a surrogate for shared antigen specificity). One can visualize GLIPH only, 10X only, or a combination of GLIPH and 10X or GLIPH and ImmunoSeq.

To plot data from GLIPH, simply provide the path to the convergence-groups.txt output table if using GLIPH version 1 or cluster.txt if using GLIPH version 2. Inner circles are represented by the T cell clonotypes and outer circles are presented by specificity groups. Additionally, if the frequency of each T cell clonotype is known (e.g. from Adaptive ImmunoSeq), that may be supplied as a 2 column data frame with column headers "clonotype" and "frequency". This will draw the clonotype circles in proportion to the frequency of the clones. If additional metadata is known about the T cell clone, then in place of or in addition to "frequency" a "data" column may be provided with discrete or continuous data and the T cell cell clonotypes will be colored accordingly.  If additional information is available at the level of an individual T cell (e.g. 10X Chromium single cell immune profiling), that may be supplied as a 3 column data frame with the column headers "clonotype", "cell", and "data". Data may be a discrete or continuous variable. In this case, individual T cell are represented by the smallest circles colored according to their data value.

The fill and line color can be adjusted and the legend may be hidden or displayed. Additionally, labels can be displayed for each cell, data value, clonotype or specificity group. Since TCellPack is based on the [ggraph](https://github.com/thomasp85/ggraph) package, plots can be modified further by adding additional layers.

### Installation instructions
```
install.packages("devtools")
devtools::install_github("davidcoffey/TCellPack")
```

### Basic usage
The package contains only one function, `PlotTCellPack()`.  The function may be run with our without GLIPH.  It will accept either GLIPH version 1 (convergence-groups.txt) or version 2 (cluster.txt) files and will detect the version automatically.  Since GLIPH version 2 allows a single clonotype to be assigned to multiple specificity groups, which is incompatible with this visualization technique, TCellPack will visualize the clonotype belonging to the largest specificity group.  In addition to or in place of GLIPH data, the user may supply clonotype.data (2-3 column data frame named "clonotype", "frequency" and/or "data") and/or cell.data (3 column data frame named "clonotype", "cell", and "data").  The data column may be discrete or numeric.

TCellPack pack comes with 4 example data frames (`gliph.example`, `clonotype.data.example`, `cell.data.continuous.example`, `cell.data.discrete.example`) that are loaded automatically. Below is the code used to create the five shown images.

```
library(TCellPack)

# Get help
?PlotTCellPack()

# T cell clone size proportional to frequency with GLIPH specificity groups
PlotTCellPack(gliph = gliph.example, clonotype.data = clonotype.data.example, legend = TRUE)

# T cell colored by discrete variable with GLIPH specificity groups
PlotTCellPack(gliph = gliph.example, cell.data = cell.data.discrete.example, legend = TRUE)

# T cell colored by continuous variable with GLIPH specificity groups
PlotTCellPack(gliph = gliph.example, cell.data = cell.data.continuous.example, legend = TRUE)

# T cell colored by discrete variable without GLIPH specificity groups
PlotTCellPack(cell.data = cell.data.discrete.example, legend = TRUE)

# T cell colored by continuous variable without GLIPH specificity groups
PlotTCellPack(cell.data = cell.data.continuous.example, legend = TRUE)
```

![](man/figures/example-plot.png)
