# Overview
`micInt` is an R package designed for analyzing microbial co-occurences. It takes in OTU tables was either pure dataFrames or as experiment-level `phyloseq` objects.

# Install
This package have the following dependencies which must be installed prior to installing this package:

From CRAN:
  * `infotheo`
  * `matrixStats`
  * `magrittr`
  * `deSolve`
  * `igraph`
  * `dplyr`
  * `rlang`
  * `glue`
  * `viridis`
  * `RhpcBLASctl`
  * `ggplot2`

[//]: # (Hello)

From Bioconductor:
* `phyloseq`

If all dependencies are satisfied, the package can we installed from GitHub as follows:
```
install.packages("devtools")
library(devtools)
install_github("AlmaasLab/micInt")
```

For a quit start automatically installing all dependencies, consider:
```
install.packages("BiocManager","devtools")
devtools::install_github("AlmaasLab/micInt",repos=BiocManager::BiocManager::repositories())
```

# Useage
Unfortunately, I have not yet taken the effort to write a proper manual. However, if you are interested in using the package, don`t hasitate asking me for help, either as a github issue or directly (jakob.p.pettersen@ntnu.no).
