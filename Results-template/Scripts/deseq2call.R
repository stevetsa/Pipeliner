## Wrapper script to grab args
## Example USAGE: Rscript deseq2call.R 'sampletable.txt' 'RawCountFile_Subread_junctions_filtered.txt' 'TCWT TGWT TCWT HCWT TGWT HGWT TCWT TCKO TGWT TGKO TCKO TGKO TCKO HCKO TGKO HGKO' 'mm10' 'projectid' 'Project Description.' 'Subread_junctions' 
args <- commandArgs(trailingOnly = TRUE)
#Sys.setenv(RSTUDIO_PANDOC="/path/to/pandoc/installation") #You may need to add this to point to pandocs locatation
rmarkdown::render("Deseq2Report.Rmd", 
  params = list(
    sampleinfo = args[1],
    data = args[2],
    contrasts = args[3],
    species = args[4],
    projectId = args[5],
    projectDesc = args[6],
    dtype = args[7]))
