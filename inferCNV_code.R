
#load Seurat object
load("all_pure_cell_seurat.RData")

library(Seurat)
table(dat.seurat.integrated.final.pure@meta.data$tumour.type)
#Control Higher-grade glioma  Lower-grade glioma 
#56654               13815               32761 

GBM.cell<-row.names(dat.seurat.integrated.final.pure@meta.data)[which(dat.seurat.integrated.final.pure$tumour.type=='Higher-grade glioma')]
LGG.cell<-row.names(dat.seurat.integrated.final.pure@meta.data)[which(dat.seurat.integrated.final.pure$tumour.type=='Lower-grade glioma')]
Con.cell<-row.names(dat.seurat.integrated.final.pure@meta.data)[which(dat.seurat.integrated.final.pure$tumour.type=='Control')]

#Downsample Control and LGG 
set.seed(1234)
LGG.cells=sample(LGG.cell,10000)
set.seed(1234)
Con.cells=sample(Con.cell,10000)

GBM.mat<-as.data.frame(GetAssayData(subset(dat.seurat.integrated.final.pure, cells=GBM.cell)))
LGG.mat<-as.data.frame(GetAssayData(subset(dat.seurat.integrated.final.pure, cells=LGG.cells)))
Con.mat<-as.data.frame(GetAssayData(subset(dat.seurat.integrated.final.pure, cells=Con.cells)))

dat=cbind(GBM.mat,LGG.mat,Con.mat)
groupinfo=data.frame(v1=colnames(dat),
                     v2=c(rep('GBM',ncol(GBM.mat)),
                          rep('LGG',ncol(LGG.mat)),
                          rep('Con',ncol(Con.mat))))

#Get gene annotation
library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

#Prepare data for infer CNV
dat=dat[rownames(dat) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

#Run inferCNV pipeline
library(infercnv)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFiles,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names=c("Con"))


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="~/Desktop/inferCNV", 
                             cluster_by_groups=FALSE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             output_format="pdf")

library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
#  Import inferCNV dendrogram
infercnv.dend <- read.dendrogram(file = "infercnv.observations_dendrogram.txt")
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 6, order_clusters_as_data = FALSE)
table(infercnv.labels)
#infercnv.labels
#1    2    3    4    5    6 
#4393  934 9617 3423  761 4687 

# Color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

infercnv.dend %>% set("labels",rep("", nobs(infercnv.dend)) )  %>% plot(main="inferCNV dendrogram") %>%
  colored_bars(colors = as.data.frame(the_bars), dend = infercnv.dend, sort_by_labels_order = FALSE, add = T, y_scale=100 , y_shift = 0)

#Get cell ids with genome alteration
infercnv.clusters=data.frame(infercnv.labels)
tumor.cell.id<-subset(infercnv.clusters,infercnv.clusters$infercnv.labels==5)

tumor.cell.seurat<-subset(dat.seurat.integrated.final.pure,cells = rownames(tumor.cell.id))
DimPlot(tumor.cell.seurat,label = F)

library(reshape2)
tumor.cell.fraction<-data.frame(table(tumor.cell.seurat@meta.data$cluster.rev,tumor.cell.seurat@meta.data$orig.ident))
tumor.cell.fraction<-dcast(tumor.cell.fraction,Var1~Var2)

write.xlsx(tumor.cell.fraction,file = "Number_of_Tumor_cells.xlsx",rowNames=F,colNames=T)

#Fig1F 
library(ggplot2)
pdf("Maglinant_cell_UMAP.pdf",width = 6,height = 6)
DimPlot(dat.seurat.integrated.final.pure,cells.highlight =as.vector(rownames(tumor.cell.id)),sizes.highlight = 0.6)+
  scale_color_manual(labels = c("No-CNV", "CNV"), values = c("grey", "red"))
dev.off()

#session info
R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X  11.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] reshape2_1.4.4    Seurat_3.1.1      ggthemes_4.2.4    dendextend_1.16.0 gridExtra_2.3     phylogram_2.1.0   infercnv_1.2.1  

