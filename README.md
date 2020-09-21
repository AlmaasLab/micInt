# Overview
`micInt` is an R package designed for analyzing microbial co-occurences. It takes in OTU tables was either pure dataFrames or as experiment-level `phyloseq` objects.

# Install
First, ensure the following packages are installed:

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

From Bioconductor:
	* `phyloseq`

Then, use:
```
# instal.packages("devtools")
library(devtools)
install_github("AlmaasLab/micInt")
```

# Useage
Unfortunately, I have not yet taken the effort to write a proper manual. However, if you are interested in using the package, don`t hasitate asking me for help, either as a github issue or directly (jakob.p.pettersen@ntnu.no).
