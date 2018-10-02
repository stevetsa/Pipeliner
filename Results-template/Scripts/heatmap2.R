###############################################################################################################################
# Create a hierarchical clustering heatmap from an input dataset. 
# The first column of the data should be a row identifier (e.g. gene symbols), 
# all other columns will be used for hierarchical clustering. 
# Needs packages Rcolorbrewer, dendsort, colorspace.
# USAGE:  Rscript heatmap2.R --data Genesets_heatmap.csv --row_id Geneset --color 'Rd Yl Bu' --display_row_names TRUE \
#                            --flip_color TRUE --z_scale FALSE --display_column_names TRUE --minimum -1 --cluster_rows TRUE \
#                            --image_height 2000 --flip_dendrogram TRUE --row_height 4 --column_width 15 --row_font 4
#/scratch/rnaseqVishalTesting/cmod_hg38pe/DEG_RSEM_genes$ Rscript heatmap2.R --data Genesets_heatmap.csv --row_id Geneset --color 'Rd Yl Bu' --display_row_names TRUE --flip_color TRUE --z_scale FALSE --display_column_names TRUE --minimum -1 --cluster_rows TRUE --image_height 2000 --flip_dendrogram TRUE --row_height 4 --column_width 15 --row_font 4
###############################################################################################################################
rm(list=ls())

library(colorspace)
library(dendsort)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(argparse)
#library(SparkR)

# create parser object
parser <- ArgumentParser()

# Adding command-line arguments
parser$add_argument("-d", "--data", type="character", required=TRUE,
                    help="Matrix of numerical data with row ids and column ids. Required.")

parser$add_argument("-i", "--row_id", type="character", required=TRUE,
                    help="The column to use as a row name. Required.")

parser$add_argument("-c", "--color", type="character", default="Default",
                    help="Heatmap colors: Required.")

parser$add_argument("-m", "--clustering_method", type="character", default="complete",
                    help="Clustering method: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid. Default value is complete.")

parser$add_argument("-k", "--cutree_row", type="integer", default=1,
                    help="Cut the row dendrogram to display this number of clusters. Default value is 1.")

parser$add_argument("-v", "--display_row_names", type="logical", default=FALSE,
                    help="Display row names. Default value is FALSE")

parser$add_argument("-l", "--cluster_columns", type="logical", default=FALSE,
                    help="Cluster Columns of the Heatmap Default value is FALSE")

parser$add_argument("-z", "--z_scale", type="logical", default=TRUE,
                    help="Perform Z-scaling before clustering. Default value is TRUE")

parser$add_argument("-f", "--flip_color", type="logical", default=FALSE,
                    help="Reverse Color Direction: Default value is FALSE")

parser$add_argument("-q", "--display_column_names", type="logical", default=FALSE,
                    help="Display column names: Default value is FALSE")

parser$add_argument("-n", "--minimum", type="integer", default=0,
                    help="Minimum: Default value is 0")

parser$add_argument("-s", "--maximum", type="integer", default=1,
                    help="Maximum: Default value is 1")

parser$add_argument("-b", "--breaks", type="logical", default=TRUE,
                    help="Get min max, or set min max for heatmap colors. Default value is TRUE.")

parser$add_argument("-e", "--distance_method", type="character", default='correlation',
                    help="Distance Method: euclidean, maximum, manhattan, canberra, binary, minkowski, correlation. Default value is correlation.")

parser$add_argument("-g", "--cluster_rows", type="logical", default=FALSE,
                    help="Cluster Rows of the Heatmap: Default value is FALSE.")

parser$add_argument("-j", "--tree_height_row", type="integer", default=50,
                    help="Row Dendrogram Height: Default value is 50.")

parser$add_argument("-a", "--tree_height_column", type="integer", default=50,
                    help="Column Dendrogram Height: Default value is 50.")

parser$add_argument("-o", "--image_height", type="integer", default=80,
                    help="Image Height: Default value is 80.")

parser$add_argument("-t", "--flip_dendrogram", type="logical",default=FALSE,
                    help="Image Height: Default value is FALSE.")

parser$add_argument("-x", "--row_height", type="integer", default=10,
                    help="Row Height: Default value is 10.")

parser$add_argument("-y", "--column_width", type="integer", default=15,
                    help="Column Width: Default value is 15.")

parser$add_argument("-w", "--row_font", type="integer", default=10,
                    help="Row Font Size: Default value is 10.")

args <- parser$parse_args()


pal=function (n, h = c(237, 43), c = 100, l = c(70, 90), power = 1, 
              fixup = TRUE, gamma = NULL, alpha = 1, ...) 
{
  if (!is.null(gamma)) 
    warning("'gamma' is deprecated and has no effect")
  if (n < 1L) 
    return(character(0L))
  h <- rep(h, length.out = 2L)
  c <- c[1L]
  l <- rep(l, length.out = 2L)
  power <- rep(power, length.out = 2L)
  rval <- seq(1, -1, length = n)
  rval <- hex(polarLUV(L = l[2L] - diff(l) * abs(rval)^power[2L], 
                       C = c * abs(rval)^power[1L], H = ifelse(rval > 0, h[1L], 
                                                               h[2L])), fixup = fixup, ...)
  if (!missing(alpha)) {
    alpha <- pmax(pmin(alpha, 1), 0)
    alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
                    width = 2L, upper.case = TRUE)
    rval <- paste(rval, alpha, sep = "")
  }
  return(rval)
}

#Color selections for heatmap:
np0=pal(100)
np1=diverge_hcl(100, c = 100, l = c(30, 80), power = 1)  #Blue to Red
np2=heat_hcl(100, c = c(80, 30), l = c(30, 90), power = c(1/5, 2))  #Red to Vanilla
np3=rev(heat_hcl(100, h = c(0, -100), c = c(40, 80), l = c(75, 40), power = 1)) #Violet to Pink
np4=colorRampPalette(brewer.pal(10, "RdYlBu"))(100)

np=list(np0,np1,np2,np3,np4)
names(np) = c("Default","Blue to Red","Red to Vanilla","Violet to Pink","Rd Yl Bu")

doheatmap <- function(dat,clus,clus2,ht,rn,cn,col){
  require(pheatmap)
  require(dendsort)
  if (args$z_scale){
    tmean.scale=t(scale(t(dat)))}
  else
  {tmean.scale=dat}
  #colnames(tmean.scale)=colnames(dat)
  col.pal <- np[[col]]
  if (args$flip_color){
    col.pal=rev(col.pal)
  }
  # define metrics for clustering
  drows1 <- args$distance_method
  dcols1 <- args$distance_method
  minx=min(tmean.scale)
  maxx=max(tmean.scale)
  
  if (args$breaks){
    breaks=seq(minx,maxx,length=100)
    legbreaks=seq(minx,maxx,length=5)}
  else {
    absmax = ceiling(max(abs(c(minx,maxx))))
    #tmean.scale %>% as.data.frame %>% 
    #  mutate_all(funs(ifelse(.>args$maximum, args$maximum, 
    #          (ifelse(.<args$minimum,args$minimum))))) -> tmean.scale
    breaks=c(-1*absmax,seq(args$minimum,args$maximum,length=98),absmax)
    #breaks=c(seq(args$minimum,args$maximum,length=100))
    #legbreaks=c(-1*absmax,seq(args$minimum,args$maximum, length=3),absmax)}
    legbreaks=c(-1*absmax,0, absmax)}
  #print ('defined metrics')
  #   annotation_col=data.frame(Group=factor(annot))
  #   print (dim(annotation_col))
  #   row.names(annotation_col) <- row.names(tmean.scale)
  #print ('set rownames of annotations')
  #Run cluster method using 
  hc=hclust(dist(t(tmean.scale)),method="average")
  hcrow=hclust(dist(tmean.scale),method="average")
  #print ('clustered the things')
  if (args$flip_dendrogram){
    sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))}
  else {
    sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))}
  if (clus2){
    rowclus <- sort_hclust(hcrow)}
  else {rowclus = FALSE}
  #print('sorted the clusters')
  hm.parameters <- list(tmean.scale, 
                        color = col.pal,
                        legend_breaks = legbreaks,
                        cellwidth = args$column_width, 
                        cellheight = args$row_height, 
                        scale = "none",
                        treeheight_col = args$tree_height_column,
                        treeheight_row = args$tree_height_row,
                        kmeans_k = NA,
                        breaks=breaks,
                        height=args$image_height,
                        fontsize_row = args$row_font,
                        show_rownames = rn, 
                        show_colnames = cn,
                        clustering_method = args$clustering_method,
                        cluster_rows = rowclus, 
                        cluster_cols = clus,
                        cutree_rows =args$cutree_row,
                        clustering_distance_rows = drows1, 
                        clustering_distance_cols = dcols1
                        #annotation_row = annotation_col, #Commenting out annotation column for now
                        #annotation_colors = annot_col
  )
  
  mat=t(tmean.scale)
  print('calculated mat')
  # print(dim(mat))
  callback = function(hc, mat){
    print('inside the callback')
    dend=rev(dendsort(as.dendrogram(hc)))
    print ('reversed the dendsorted hc')
    as.hclust(dend)
  }
  do.call("pheatmap", c(hm.parameters, list(clustering_callback = callback)))
}

df.orig=read.table(args$data, header=TRUE, sep=",")
rowid <- as.character(args$row_id)
df.orig %>% group_by(!!as.name(args$row_id)) %>% summarise_all(funs(mean)) -> df
print(df)
print(df[, rowid])
#print(colnames(df))
df.mat = data.frame(df, row.names=rowid)
print(dim(df))
print(dim(df.mat))
#rownames(df.mat) <- as.character(df[, args$row_id]) 
print(row.names(df.mat))
print(head(df.mat))

#example:
#doheatmap("408genes",finalres_heat,f.patient$label,ann_colors,TRUE,1,FALSE,1)
#print ('this is actually the first thing that should appear')
png("heatmap.png")
p=doheatmap(dat=df.mat, clus=args$cluster_columns,clus2=args$cluster_rows,ht=10,rn=args$display_row_names,cn=args$display_column_names,col=args$color)
if(args$cluster_rows == TRUE){
  print(row.names(df)[p$tree_row$order])
}


dev.off()


