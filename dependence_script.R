# The purpose of this script is to investigate dependencies and removing
# unnecessary ones
library(DependenciesGraphs)
library(stringr)
library(micInt)
par(ask=TRUE)
for (object in getNamespaceExports('micInt')){
  if(stringr::str_detect(object,'\\._')){
    next()
  }
  print(object)
  deps <- DependenciesGraphs::allDepFunction('stringr',object)
  print(deps$Nomfun$label)
}
