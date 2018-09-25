install.packages(c("devtools", "RCurl", "XML", "reshape", "ggplot2", "RColorBrewer", "Matrix", "UpSetR", "snow", "dplyr", "gdata", "statmod", "methods", "rgl",
"Seurat", "amap", "rmarkdown", "plotly", "knitr", "spp", "gplots", "ggfortify", "DT"), dependencies=TRUE)
source("https://bioconductor.org/biocLite.R")
biocLite(c("devtools", "remotes"))
biocLite("rtracklayer")
library(rtracklayer)
biocLite(c("edgeR", "limma", "scran", "GSEABase", "GSVA", "ReactomePA", "clusterProfiler", "genefilter", "maftools", "org.Hs.eg.db"))
biocLite(c("COMBINE-lab/wasabi", "geneplotter", "pachterlab/sleuth", "DOSE", "DESeq2", "scater", "biomaRt", "ChIPseeker", "GenomicFeatures", "AnnotationDbi"))
