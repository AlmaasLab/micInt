# Overview
`micInt` is an R package designed for analyzing microbial co-occurences. It takes in OTU tables was either pure dataFrames or as experiment-level `phyloseq` objects. Roughly, the package is composed of three major parts:
   * CoNet-based analysis of pairwise co-occurences. Primary functions: `runAnalysis` (high level) and `ccrepe` (low level)
   * Lotka-Volterra modelling of time dynamics. Primary functions: `integralSystem`, `cv.LV` and `ridge_fit`
   * Utility functions, includes `refine_data`, `subset_by_environment`, `similarity_measures` and `scale_by_column`.

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
  * `ggfortify`

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

Please read the article "Robust bacterial co-occurence community structures are independent of r- and K-selection history" (https://www.nature.com/articles/s41598-021-03018-z) where this package is used. The source code for the article is available on https://github.com/yaccos/Microbial-co-occurence.

Unfortunately, I have not yet taken the effort to write a proper vignette. However, if you are interested in using the package, don`t hasitate asking me for help, either as a github issue or directly (jakob.p.pettersen@ntnu.no).
