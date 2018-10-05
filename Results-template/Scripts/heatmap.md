```R
install.packages("reticulate")
library(reticulate)
```
    Updating HTML index of packages in '.Library'
    Making 'packages.html' ... done
```R
library(reticulate)
#py_install("pandas")
#create environment
#conda_create("env")
# install SciPy
conda_install("env", "python=3.6")
```
```R
install.packages(c("dendsort", "argparse"), dependencies=TRUE)
```
    Updating HTML index of packages in '.Library'
    Making 'packages.html' ... done
```R
rm(list=ls())
library(colorspace)
library(dendsort)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)
library(argparse)
```
```R
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
```
```R
df.orig=read.table("/sbgenomics/project-files/Genesets_heatmap.csv", header=TRUE, sep=",")
#rowid <- as.character(args$row_id)
rowid <- "Geneset"
df.orig %>% group_by(!!as.name(rowid)) %>% summarise_all(funs(mean)) -> df
#print(df)
#print(df[, rowid])
#print(colnames(df))
df.mat = data.frame(df, row.names=rowid)
#print(dim(df))
#print(dim(df.mat))
##rownames(df.mat) <- as.character(df[, args$row_id]) 
#print(row.names(df.mat))
#print(head(df.mat))
dat<-df.mat
```
```R
z_scale<-"FALSE"
flip_color<-"TRUE"
distance_method<-"correlation"
breaks<-"TRUE"
flip_dendrogram<-"TRUE"
column_width<-10
row_height<-4
minimum<--1
maximum<-1
tree_height_column<-50
tree_height_row<-50
image_height<-2000
row_font<-4
clustering_method<-"complete"
png("heatmap.png")
clus<-"FALSE"
clus2<-"TRUE"
#ht=10
#rn<-"TRUE"
#cn<-"TRUE"
col<-"Rd Yl Bu"
if (z_scale){
    tmean.scale=t(scale(t(dat)))
} else {
    tmean.scale=dat
  }
  #colnames(tmean.scale)=colnames(dat)
col.pal <- np[[col]]
  
if (flip_color){
    col.pal=rev(col.pal)
}
drows1 <- distance_method
dcols1 <- distance_method
minx=min(tmean.scale)
maxx=max(tmean.scale)
if (breaks){
    breaks=seq(minx,maxx,length=100)
    legbreaks=seq(minx,maxx,length=5)
} else {
    absmax = ceiling(max(abs(c(minx,maxx))))
    #tmean.scale %>% as.data.frame %>% 
    #  mutate_all(funs(ifelse(.>args$maximum, args$maximum, 
    #          (ifelse(.<args$minimum,args$minimum))))) -> tmean.scale
    breaks=c(-1*absmax,seq(minimum,maximum,length=98),absmax)
    #breaks=c(seq(args$minimum,args$maximum,length=100))
    #legbreaks=c(-1*absmax,seq(args$minimum,args$maximum, length=3),absmax)}
    legbreaks=c(-1*absmax,0, absmax)
}
hc=hclust(dist(t(tmean.scale)),method="average")
hcrow=hclust(dist(tmean.scale),method="average")
  #print ('clustered the things')
if (flip_dendrogram){
    sort_hclust <- function(...) as.hclust(rev(dendsort(as.dendrogram(...))))
    } else {
      sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
    }
  
if (clus2){
    rowclus <- sort_hclust(hcrow)
} else {rowclus = FALSE
       }
```
```R
#png("heatmap.png", height = , width = 8, units  = "in", res=1200, pointsize = 4)
png("heatmap.png", height = 6, width = 6, units  = "in", res=600)
#ggsave( "ggtest.png", width = 3.25, height = 3.25, dpi = 1200)
p=pheatmap(df.mat,color = col.pal, height=10000, legend_breaks = legbreaks, fontsize_row =8, treeheight_col = 10, treeheight_row =10,  show_rownames =TRUE, 
                        cellwidth = column_width, 
                        cellheight = row_height, 
                        scale = "none",
                        kmeans_k = NA,
                        breaks=breaks,
                        clustering_method = clustering_method,
            #            cluster_rows = rowclus, 
            #            cluster_cols = clus,
                        cutree_rows = 1,
                        clustering_distance_rows = drows1, 
                        clustering_distance_cols = dcols1,
                        cex = 0.65
          )
dev.off()
#if(clus2 == "TRUE"){
# print(row.names(df)[p$tree_row$order])
#}
```
<strong>png:</strong> 2
