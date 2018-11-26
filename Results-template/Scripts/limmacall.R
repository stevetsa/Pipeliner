## Wrapper script to grab args
## EXAMPLE USAGE: Rscript limmacall.R 'sampletable.txt' 'RawCountFile_RSEM_genes_filtered.txt' 'KO WT' 'mm10' 'limmatest' 'Enter CCBR Project Description and Notes here.' 'RSEM_genes'
## Parse Args
args <- commandArgs(trailingOnly = TRUE)
# Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")
rmarkdown::render("LimmaReport.Rmd", params = list(
    sampleinfo = args[1],
    data = args[2],
    contrasts = args[3],
    species = args[4],
    projectId = args[5],
    projectDesc = args[6],
    dtype = args[7]
  )
)
