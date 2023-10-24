library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
setwd("/Volumes/Elements/7human_sample_202111/13_sample_seurat")


#############################################################################################################
###################################### all cells filtering and clustering ###################################
#############################################################################################################
##load data
control1<-Read10X(data.dir = "N027/filtered_feature_bc_matrix")
tumour1<-Read10X(data.dir = "T027/filtered_feature_bc_matrix")

control2<-Read10X(data.dir = "NO36/filtered_feature_bc_matrix")

control3<-Read10X(data.dir = "N037/filtered_feature_bc_matrix")
tumour3<-Read10X(data.dir ="T037/filtered_feature_bc_matrix")

control4<-Read10X(data.dir ="N038/filtered_feature_bc_matrix")
tumour4<-Read10X(data.dir = "T038/filtered_feature_bc_matrix")

control5<-Read10X(data.dir ="N043/filtered_feature_bc_matrix")
tumour5<-Read10X(data.dir = "T043/filtered_feature_bc_matrix")

control6<-Read10X(data.dir ="N045/filtered_feature_bc_matrix")
tumour6<-Read10X(data.dir = "T045/filtered_feature_bc_matrix")

control7<-Read10X(data.dir ="N046/filtered_feature_bc_matrix")
tumour7<-Read10X(data.dir = "T046/filtered_feature_bc_matrix")

##make names specific
colnames(control1)<-paste0("control1_",colnames(control1))
colnames(control2)<-paste0("control2_",colnames(control2))
colnames(control3)<-paste0("control3_",colnames(control3))
colnames(control4)<-paste0("control4_",colnames(control4))
colnames(control5)<-paste0("control5_",colnames(control5))
colnames(control6)<-paste0("control6_",colnames(control6))
colnames(control7)<-paste0("control7_",colnames(control7))
colnames(tumour1)<-paste0("tumour1_",colnames(tumour1))
colnames(tumour3)<-paste0("tumour3_",colnames(tumour3))
colnames(tumour4)<-paste0("tumour4_",colnames(tumour4))
colnames(tumour5)<-paste0("tumour5_",colnames(tumour5))
colnames(tumour6)<-paste0("tumour6_",colnames(tumour6))
colnames(tumour7)<-paste0("tumour7_",colnames(tumour7))

##creat seurat object & Standard pre-processing
Rep1.dat=cbind(tumour1,control1)
Rep1 <- CreateSeuratObject(counts = Rep1.dat,min.cells = 3,min.features = 200)
Rep1[["percent.mt"]] <- PercentageFeatureSet(Rep1, pattern = "^MT-")
Rep1 <- subset(Rep1, subset =  percent.mt < 10)
Rep1 <- NormalizeData(Rep1, scale.factor = 10000)
Rep1$genotype <- "Rep1"
dim(GetAssayData(object = Rep1, slot = "counts"))
#[1] 21725 25527
save(Rep1,file = "Rep1.RData")


Rep2 <- CreateSeuratObject(counts = control2,min.cells = 3,min.features = 200)
Rep2[["percent.mt"]] <- PercentageFeatureSet(Rep2, pattern = "^MT-")
Rep2 <- subset(Rep2, subset =  percent.mt < 10)
Rep2 <- NormalizeData(Rep2, scale.factor = 10000)
Rep2$genotype <- "Rep2"
dim(GetAssayData(object = Rep2, slot = "counts"))
#[1] 21811 30572
save(Rep2,file = "Rep2.RData")

Rep3.dat=cbind(tumour3,control3)
Rep3 <- CreateSeuratObject(counts = Rep3.dat, min.cells = 3,min.features = 200)
Rep3[["percent.mt"]] <- PercentageFeatureSet(Rep3, pattern = "^MT-")
Rep3 <- subset(Rep3, subset =  percent.mt < 10)
Rep3 <- NormalizeData(Rep3, scale.factor = 10000)
Rep3$genotype <- "Rep3"
dim(GetAssayData(object = Rep3, slot = "counts"))
#[1] 24062 40022
save(Rep3,file = "Rep3.RData")

Rep4.dat=cbind(tumour4,control4)
Rep4 <- CreateSeuratObject(counts = Rep4.dat,min.cells = 3,min.features = 200)
Rep4[["percent.mt"]] <- PercentageFeatureSet(Rep4, pattern = "^MT-")
Rep4 <- subset(Rep4, subset =  percent.mt < 10)
Rep4 <- NormalizeData(Rep4, scale.factor = 10000)
Rep4$genotype <- "Rep4"
dim(GetAssayData(object = Rep4, slot = "counts"))
#[1] 24428 22256
save(Rep4,file = "Rep4.RData")

Rep5.dat=cbind(tumour5,control5)
Rep5 <- CreateSeuratObject(counts = Rep5.dat,min.cells = 3,min.features = 200)
Rep5[["percent.mt"]] <- PercentageFeatureSet(Rep5, pattern = "^MT-")
Rep5 <- subset(Rep5, subset =  percent.mt < 10)
Rep5 <- NormalizeData(Rep5, scale.factor = 10000)
Rep5$genotype <- "Rep5"
dim(GetAssayData(object = Rep5, slot = "counts"))
#[1] 22306 16641
save(Rep5,file = "Rep5.RData")

Rep6.dat=cbind(tumour6,control6)
Rep6 <- CreateSeuratObject(counts = Rep6.dat,min.cells = 3,min.features = 200)
Rep6[["percent.mt"]] <- PercentageFeatureSet(Rep6, pattern = "^MT-")
Rep6 <- subset(Rep6, subset =  percent.mt < 10)
Rep6 <- NormalizeData(Rep6, scale.factor = 10000)
Rep6$genotype <- "Rep6"
dim(GetAssayData(object = Rep6, slot = "counts"))
#[1] 23209 31229
save(Rep6,file = "Rep6.RData")

Rep7.dat=cbind(tumour7,control7)
Rep7 <- CreateSeuratObject(counts = Rep7.dat,min.cells = 3,min.features = 200)
Rep7[["percent.mt"]] <- PercentageFeatureSet(Rep7, pattern = "^MT-")
Rep7 <- subset(Rep7, subset =  percent.mt < 10)
Rep7 <- NormalizeData(Rep7, scale.factor = 10000)
Rep7$genotype <- "Rep7"
dim(GetAssayData(object = Rep7, slot = "counts"))
#[1] 23606 22146
save(Rep7,file = "Rep7.RData")

##combine 7samples
dat.combined=merge(x = Rep1, y = c(Rep2, Rep3, Rep4, Rep5, Rep6, Rep7))
rm(list=setdiff(ls(), "dat.combined"))

dat.combined <- NormalizeData(dat.combined, scale.factor = 10000)
dat.combined=FindVariableFeatures(dat.combined, selection.method = "vst", nfeatures = 2000)

dat.combined <- ScaleData(dat.combined)

dat.combined <- RunPCA(dat.combined, features = VariableFeatures(object = dat.combined))
#Clustering
dat.combined <- FindNeighbors(dat.combined, reduction = "pca", dims = 1:30)
dat.combined <- FindClusters(dat.combined, resolution = 0.5)
# UMAP after find clusters
dat.combined <- RunUMAP(dat.combined, reduction = "pca", dims = 1:30)

pdf(file = "UMAP_no_align_r0.5.pdf",width = 12,height = 8)
P1<-DimPlot(dat.combined,reduction = "umap",label=T)
P2=DimPlot(dat.combined,reduction = "umap",group.by = "genotype")
plot_grid(P1,P2)
dev.off()
save(dat.combined,file = "dat.combined.RData")


library(future)
library(future.apply)
plan("multiprocess", workers = 3)
options(future.globals.maxSize = 20000 * 1024^2)

#integrated 13 samples by rPCA
dat.seurat.list <- SplitObject(dat.combined, split.by = "genotype")
rm(list=setdiff(ls(), "dat.seurat.list"))
dat.seurat.list <- future_lapply(X = dat.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = dat.seurat.list)
dat.seurat.list <- future_lapply(X = dat.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
save(dat.seurat.list,file = "dat.seurat.list.RData")
anchors <- FindIntegrationAnchors(object.list = dat.seurat.list, reduction = "rpca", 
                                  dims = 1:30)
save(anchors, file="anchors.RData")

rm(list=setdiff(ls(), "anchors"))
dat.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
save(dat.seurat.integrated, file="dat.seurat.raw.integrated.RData")
rm(list=setdiff(ls(), "dat.seurat.integrated"))


dat.seurat.integrated<-ScaleData(dat.seurat.integrated,verbose=F)
dat.seurat.integrated<-RunPCA(dat.seurat.integrated,verbose=F)
dat.seurat.integrated<-RunUMAP(dat.seurat.integrated,dims=1:30)

dat.seurat.integrated<-FindNeighbors(dat.seurat.integrated,reduction="pca",dims=1:30)
dat.seurat.integrated<-FindClusters(dat.seurat.integrated,resolution=1)

dat.seurat.integrated@meta.data$stim<-substring(dat.seurat.integrated@meta.data$orig.ident,1,3)

Idents(dat.seurat.integrated.f)<-"seurat_clusters"
pdf(file = "UMAP_aligned_pc30_r1_rPCA.pdf",width = 12,height = 8)
p1<-DimPlot(dat.seurat.integrated,label = T)
p2<-DimPlot(dat.seurat.integrated,group.by = "stim")+scale_color_discrete(labels=c("Control","Tumor"))
plot_grid(p1,p2)
DimPlot(dat.seurat.integrated,split.by = "genotype")
dev.off()

save(dat.seurat.integrated, file="dat.seurat.integrated.pc30.r1.clustered.RData")

#Doublet filtering
library(DoubletFinder)
#work on individual dataset
dat.seurat.list <- SplitObject(dat.seurat.integrated, split.by = "orig.ident")

###
library(future)
library(future.apply)
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 20000 * 1024^2)

dat.seurat.list <- future_lapply(X = dat.seurat.list, FUN = function(x) {
  DefaultAssay(x)="RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:10)
  
  ## pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep_v3(x, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_value <- as.numeric(as.character(bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)]))
  
  ## Homotypic Doublet Proportion Estimate 
  annotations <- x@meta.data$integrated_snn_res.1
  homotypic.prop <- modelHomotypic(annotations)  
  
  ## Assuming 0.8%/2 doublet formation rate
  DoubletRate<-ncol(x)*4*1e-6
  
  nExp_poi <- round(DoubletRate*length(colnames(x)))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder 
  pN_value <- 0.25
  pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
  x <- doubletFinder_v3(x, PCs = 1:10, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
})

save(dat.seurat.list, file="seurat.list.doublet.analysis.RData")

####check doublet filtering result
table(dat.seurat.list[[1]]@meta.data[,12])
#Doublet Singlet 
#2497   22488 
VlnPlot(dat.seurat.list[[1]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.07_2497",pt.size = 0)

table(dat.seurat.list[[2]]@meta.data[,12])
#Doublet Singlet 
#1     541 

table(dat.seurat.list[[3]]@meta.data[,12])
#Doublet Singlet 
#3739   26833 
VlnPlot(dat.seurat.list[[3]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.21_3739",pt.size = 0)

table(dat.seurat.list[[4]]@meta.data[,12])
#Doublet Singlet 
#2116   20886
VlnPlot(dat.seurat.list[[4]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.29_2116",pt.size = 0)

table(dat.seurat.list[[5]]@meta.data[,12])
#Doublet Singlet 
#1159   15861 
VlnPlot(dat.seurat.list[[5]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.08_1159",pt.size = 0)

table(dat.seurat.list[[6]]@meta.data[,12])
#Doublet Singlet 
#163    6214 
VlnPlot(dat.seurat.list[[6]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.2_163",pt.size = 0)

table(dat.seurat.list[[7]]@meta.data[,12])
#Doublet Singlet 
#1009   14870 
VlnPlot(dat.seurat.list[[7]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.1_1009",pt.size=0)

table(dat.seurat.list[[8]]@meta.data[,12])
#Doublet Singlet 
#260    7796 
VlnPlot(dat.seurat.list[[8]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.27_260",pt.size=0)

table(dat.seurat.list[[9]]@meta.data[,12])
#Doublet Singlet 
#295    8290 
VlnPlot(dat.seurat.list[[9]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.005_295",pt.size=0)

table(dat.seurat.list[[10]]@meta.data[,12])
#Doublet Singlet 
#274    8002 
VlnPlot(dat.seurat.list[[10]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.02_274",pt.size=0)

table(dat.seurat.list[[11]]@meta.data[,12])
#Doublet Singlet 
#2107   20846 
VlnPlot(dat.seurat.list[[11]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.08_2107",pt.size=0)

table(dat.seurat.list[[12]]@meta.data[,12])
#Doublet Singlet 
#739   12852 
VlnPlot(dat.seurat.list[[12]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.05_739",pt.size=0)

table(dat.seurat.list[[13]]@meta.data[,12])
#Doublet Singlet 
#293    8262 
VlnPlot(dat.seurat.list[[13]],features = c("nCount_RNA"),group.by = "DF.classifications_0.25_0.09_293",pt.size=0)

###
dat.seurat.doublefinder.result=lapply(dat.seurat.list, function(x){
  data.frame(cell.id=rownames(x@meta.data), doublet.finder=x@meta.data[,12])
})

dat.seurat.doublefinder.result=do.call(rbind, dat.seurat.doublefinder.result)
dim(dat.seurat.doublefinder.result)
#[1] 188393      2
rownames(dat.seurat.doublefinder.result)=dat.seurat.doublefinder.result$cell.id

save(dat.seurat.doublefinder.result, file = "all13.dataset.doubletfinder.result.RData")

#get singlet cell ids
load("all13.dataset.doubletfinder.result.RData")
cell.id.single=with(dat.seurat.doublefinder.result, as.vector(cell.id)[doublet.finder=="Singlet"])

dat.seurat.integrated.f=subset(dat.seurat.integrated, cells = cell.id.single)
dim(dat.seurat.integrated.f)
#[1]  26197 173741

table(Idents(dat.seurat.integrated.f))
#0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15 
#16661 15728 11986 12571 11655 11076  8748  7435  9621  8662  6189  7001  5640  5633  5533  3101 
#16    17    18    19    20    21    22    23    24    25    26    27    28    29    30 
#4565  3216  2602  2568  2320  2026  1875  1913  1698   994   949   706   611   404    54 

DimPlot(dat.seurat.integrated.f,reduction = "umap", label = TRUE)

dat.seurat.integrated.f@meta.data[, "experiment"]=substring(dat.seurat.integrated.f@meta.data$orig.ident,1,3)
table(dat.seurat.integrated.f@meta.data$experiment)
#con   tum 
#95503 78238 

pdf(file = "UMAP_aligned_pc30_r1_rPCA.pdf",width = 12,height = 8)
p1<-DimPlot(dat.seurat.integrated.f,label = T)
p2<-DimPlot(dat.seurat.integrated.f,group.by = "stim")+scale_color_discrete(labels=c("Control","Tumor"))
plot_grid(p1,p2)
DimPlot(dat.seurat.integrated.f,split.by = "genotype")
dev.off()

save(dat.seurat.integrated.f, file="dat.seurat.integrated.r1.clustered.filtered.RData")

#QC
VlnPlot(dat.seurat.integrated.f,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)

#Find cluster markers 
all.markers<-FindAllMarkers(dat.seurat.integrated.f,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
nk=length(unique(all.markers$cluster))

all.marker.genes=sapply(0:(nk-1), function(i){
  t.k=all.markers[all.markers$cluster==i,]
  t.k.genes=with(t.k, gene[p_val_adj <0.05 ])
  t.k.genes=head(t.k.genes, 100)
  paste(t.k.genes, collapse = " ")
})
all.marker.genes=data.frame(cluster=0:(nk-1), genes=all.marker.genes)
write.xlsx(all.marker.genes,file = "top100_cluster_marker_genes.xlsx")

#cluster annotation
seurat.cluster.rev=as.numeric(as.vector(dat.seurat.integrated.f@meta.data$seurat_clusters))
#table(seurat.cluster.rev)
#seurat.cluster.rev
#0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15 
#16661 15728 11986 12571 11655 11076  8748  7435  9621  8662  6189  7001  5640  5633  5533  3101 
#16    17    18    19    20    21    22    23    24    25    26    27    28    29    30 
#4565  3216  2602  2568  2320  2026  1875  1913  1698   994   949   706   611   404    54 

seurat.cluster.rev[seurat.cluster.rev %in% c(23)]="Macrophage"
seurat.cluster.rev[seurat.cluster.rev %in% c(0,1,2,4,5,9,16,19,22)]="Endothelial"
seurat.cluster.rev[seurat.cluster.rev %in% c(3,11,13,14,17,20,26)]="Mural cell"
seurat.cluster.rev[seurat.cluster.rev %in% c(18)]="Perivascular Fibroblast"
seurat.cluster.rev[seurat.cluster.rev %in% c(28)]="High Matrix cells"
seurat.cluster.rev[seurat.cluster.rev %in% c(10)]="Meningeal Fibroblast"
seurat.cluster.rev[seurat.cluster.rev %in% c(8,12,21)]="Oligodendrocyte"
seurat.cluster.rev[seurat.cluster.rev %in% c(24)]="Astrocyte"
seurat.cluster.rev[seurat.cluster.rev %in% c(25)]="Proliferating cell"
seurat.cluster.rev[seurat.cluster.rev %in% c(27)]="T cell"
seurat.cluster.rev[seurat.cluster.rev %in% c(6,7,15,29,30)]="Contaminated and low quanlity cell"

table(seurat.cluster.rev)
#seurat.cluster.rev
#Astrocyte Contaminated and low quanlity cell                        Endothelial 
#1698                              19742                              84776 
#High Matrix cells                         Macrophage               Meningeal Fibroblast 
#611                               1913                               6189 
#Mural cell                    Oligodendrocyte            Perivascular Fibroblast 
#37223                              17287                               2602 
#Proliferating cell                             T cell 
#994                                706 

dat.seurat.integrated.f@meta.data[, "cluster.rev"]=seurat.cluster.rev
pdf(file="UMAP_aligned_r1_revised.pdf",width=8,height=6)
DimPlot(dat.seurat.integrated.f, reduction = "umap", group.by = "cluster.rev")
dev.off()
save(dat.seurat.integrated.f,file = "dat.seurat.integrated.pc30.r1.clustered.filtered.revised.RData")

#filter low quality and contamination cluster
dat.seurat.integrated.final<-subset(dat.seurat.integrated.f,idents = c(0:5,8:14,16:28))
DimPlot(dat.seurat.integrated.final, reduction = "umap", group.by = "cluster.rev")
#tumour type annotation
tumour.type.rev=as.vector(dat.seurat.integrated.final@meta.data$orig.ident)
table(tumour.type.rev)
#tumour.type.rev
#control1 control2 control3 control4 control5 control6 control7  tumour1  tumour3  tumour4  tumour5 
#496    21399    14906    13992     8106    16218     8156    21276    18438     6007     7358 
#tumour6  tumour7 
#7163    10484 


tumour.type.rev[tumour.type.rev %in% c("control1","control2","control3","control4","control5","control6","control7")]="Control"
tumour.type.rev[tumour.type.rev %in% c("tumour1","tumour3","tumour5")]="Lower-grade glioma"
tumour.type.rev[tumour.type.rev %in% c("tumour4","tumour6","tumour7")]="Higher-grade glioma"


table(tumour.type.rev)
#tumour.type.rev
#Control Higher-grade glioma  Lower-grade glioma 
#83273               23654               47072 

dat.seurat.integrated.final@meta.data[, "tumour.type"]=tumour.type.rev

#stimulation annotation
stim.rev=as.vector(dat.seurat.integrated.final@meta.data$stim)
table(stim.rev)
#stim.rev
#Control  Tumour 
#83273   70726 

stim.rev[stim.rev%in%c("con")]="Control"
stim.rev[stim.rev%in%c("tum")]="Tumour"
table(stim.rev)
#stim.rev
#Control  Tumour 
#83273   70726 

dat.seurat.integrated.final@meta.data[, "stim"]=stim.rev

#UMAP after filtering
pdf(file = "UMAP_plot_by_sample_filtered_final.pdf",width = 8,height = 6)
DimPlot(dat.seurat.integrated.final,group.by = "cluster.rev")
stim_type_cols <- c("#619CFF","#F8766D")
DimPlot(dat.seurat.integrated.final,group.by = "stim",cols = stim_type_cols)
DimPlot(dat.seurat.integrated.final,group.by = "stim",split.by = "stim",label = F,cols=stim_type_cols)+
  theme(legend.position = "none")
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
DimPlot(dat.seurat.integrated.final,group.by = "tumour.type",cols = tumour.type.color)
DimPlot(dat.seurat.integrated.final,group.by = "tumour.type",split.by = "tumour.type",cols=tumour.type.color)+
  theme(legend.position = "none")
dev.off()


save(dat.seurat.integrated.final,file = "dat.seurat.final.RData")



#############################################################################################################
######################################## Endothelial subcluster #############################################
#############################################################################################################
Idents(dat.seurat.integrated.final)<-"cluster.rev"
EC.seurat<-subset(dat.seurat.integrated.final,idents = c("Endothelial"))
save(EC.seurat,file = "EC_all_seurat.RData")

Idents(EC.seurat)<-"stim"
Control.EC.seurat<-subset(EC.seurat,idents = "Control")
dim(Control.EC.seurat)
#[1] 26197 45205
Tumour.EC.seurat<-subset(EC.seurat,idents = "Tumour")
dim(Tumour.EC.seurat)
#[1] 26197 39571
save(Control.EC.seurat,file = "EC_Control_seurat_raw.RData")
save(Tumour.EC.seurat,file = "EC_Tumour_seurat_raw.RData")

DefaultAssay(EC.seurat)<-"RNA"

EC.seurat = NormalizeData(EC.seurat, scale.factor =10000)
EC.seurat=FindVariableFeatures(EC.seurat, selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(EC.seurat)
EC.seurat <- ScaleData(EC.seurat, features = all.genes)

EC.seurat <- RunPCA(EC.seurat, features = VariableFeatures(object = EC.seurat))
EC.seurat <- FindNeighbors(EC.seurat, reduction = "pca", dims = 1:30)
EC.seurat <- FindClusters(EC.seurat, resolution = 0.5)
EC.seurat <- RunUMAP(EC.seurat, reduction = "pca", dims = 1:30)

DimPlot(EC.seurat, reduction = "umap", label = T)

##Do integration by rPCA
library(future)
library(future.apply)
plan("multiprocess", workers = 2)
options(future.globals.maxSize = 20000 * 1024^2)

EC.seurat.list <- SplitObject(EC.seurat, split.by = "genotype")
rm(list=setdiff(ls(), "EC.seurat.list"))

EC.seurat.list <- future_lapply(X = EC.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = EC.seurat.list)
EC.seurat.list <- future_lapply(X = EC.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
save(EC.seurat.list,file = "EC.seurat.list.RData")
anchors <- FindIntegrationAnchors(object.list = EC.seurat.list, reduction = "rpca", 
                                  dims = 1:50)
save(anchors, file="EC.anchors.RData")

rm(list=setdiff(ls(), "anchors"))
EC.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
dim(EC.seurat.integrated)
#[1] 26197 84776
save(EC.seurat.integrated, file="EC.seurat.raw.integrated.RData")
rm(list=setdiff(ls(), "EC.seurat.integrated"))

##Run UMAP and Clustering
EC.seurat.integrated<-ScaleData(EC.seurat.integrated,verbose=F)
EC.seurat.integrated<-RunPCA(EC.seurat.integrated,verbose=F)
EC.seurat.integrated<-RunUMAP(EC.seurat.integrated,dims=1:15)

EC.seurat.integrated<-FindNeighbors(EC.seurat.integrated,reduction="pca",dims=1:15)

##using clustree package to find best clustering resolution
seq <- seq(0.1, 1, by = 0.1)
library(clustree)
EC.seurat.integrated<-FindClusters(EC.seurat.integrated,resolution=seq)
clustree(EC.seurat.integrated, prefix = 'integrated_snn_res.') + coord_flip()
DimPlot(EC.seurat.integrated,label = T,pt.size = 0.1,group.by = "integrated_snn_res.0.6")
#set 0.8 as the final resolution
Idents(EC.seurat.integrated)<-"integrated_snn_res.0.8"

#visulization
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
pdf("EC_subcluster_UMAP.pdf",width = 8,height = 4)
p1<-DimPlot(EC.seurat.integrated,label = T,pt.size = 0.1)
p2<-DimPlot(EC.seurat.integrated,group.by = "stim",cols = c("#619CFF","#F8766D"),pt.size = 0.1)
print(p1+p2)
DimPlot(EC.seurat.integrated,split.by = "stim",pt.size = 0.1)
DimPlot(EC.seurat.integrated,split.by = "genotype",pt.size = 0.1)
DimPlot(EC.seurat.integrated,split.by = "tumour.type",pt.size = 0.1)
DimPlot(EC.seurat.integrated,split.by = "tumour.type",group.by = "tumour.type",cols = tumour.type.color,pt.size = 0.1)
dev.off()

##EC subcluster QC
pdf("EC_subcluster_QC.pdf",width = 6,height = 4)
VlnPlot(EC.seurat.integrated,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
dev.off()

##Find EC subcluster markers
all.markers<-FindAllMarkers(EC.seurat.integrated,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

top100.markers<-all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_logFC)

write.xlsx(top100.markers,file = "top100.EC.markers.data.xlsx",rowNames=T,colNames=T)
write.xlsx(all.markers,file = "all.EC.markers.xlsx",rowNames=T,colNames=T)

save(EC.seurat.integrated,file = "EC_rpca_subcluster_pc15_r0.8_clustered.RData")

##filter low quality EC sub-clusters
pure.EC.Seurat<-subset(EC.seurat.integrated,idents = c(0,1,3,5,6,9,10,11,12,13,14))
DefaultAssay(pure.EC.Seurat)<-"RNA"

dim(pure.EC.Seurat)
#[1] 26197 57324
DimPlot(pure.EC.Seurat,label = T,pt.size = 0.1,group.by = "cluster.rev")

##EC sub-cluster annotation
pure.EC.cluster.rev=as.numeric(as.vector(pure.EC.Seurat@meta.data$seurat_clusters))
#pure.EC.cluster.rev
#0    1    3    5    6    9   10   11   12   13   14 
#9895 9613 8291 7028 6852 4689 4474 2393 1736 1653  700 

pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(0,1,3,9)]="Capillary"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(5)]="Artery"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(10)]="Large Artery"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(6)]="Venule"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(11)]="Vein"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(12)]="GBM-EC type1"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(13)]="GBM-EC type2"
pure.EC.cluster.rev[pure.EC.cluster.rev %in% c(14)]="Inflamatory EC"

table(pure.EC.cluster.rev)
#pure.EC.cluster.rev
#Artery      Capillary Inflamatory EC   Large Artery           TEC1           TEC2 
#7028          32488            700           4474           1736           1653 
#Vein         Venule 
#2393           6852 

pure.EC.Seurat@meta.data[,"cluster.rev"]<-pure.EC.cluster.rev

##check cell fraction
library(reshape2)
patient.table<-data.frame(table(pure.EC.Seurat@meta.data$genotype,pure.EC.Seurat@meta.data$cluster.rev))
colnames(patient.table)<-c("Patient","CellType","Number")
patient.table <- dcast(data = patient.table,Patient~CellType)   
#Artery Capillary Inflamatory EC Large Artery TEC1 TEC2 Vein Venule
#Sample027   1858      5580             72         1240  133  210  503   1291
#Sample036    576      2978            158          242   32   39  108    509
#Sample037   1959      7905             54         1189    8   48  597   2689
#Sample038    869      4782            118          670  354  123  168    666
#Sample043    654      3414            104          399   14   64  150    486
#Sample045    602      4454            104          363  593 1028  453    709
#Sample046    510      3375             90          371  602  141  414    502
write.xlsx(patient.table,file = "EC_subcluster_patient_propotion.xlsx",rowNames=F,colNames=T)

sample.table<-data.frame(table(pure.EC.Seurat@meta.data$orig.ident,pure.EC.Seurat@meta.data$cluster.rev))
colnames(sample.table)<-c("Sample","CellType","Number")
sample.table <- dcast(data = sample.table,Sample~CellType)  
#Artery Capillary Inflamatory EC Large Artery TEC1 TEC2 Vein Venule
#control1     29       206              0           25    0    0   13     15
#control2    576      2978            158          242   32   39  108    509
#control3   1137      3613             17          617    2   13  355    679
#control4    729      4332            116          532    5    1  109    593
#control5    352      1811             37          240    2    1   38    197
#control6    441      3181             96          239   36  221  194    676
#control7    379      3106             83          188    0    8   31    473
#tumour1    1829      5374             72         1215  133  210  490   1276
#tumour3     822      4292             37          572    6   35  242   2010
#tumour4     140       450              2          138  349  122   59     73
#tumour5     302      1603             67          159   12   63  112    289
#tumour6     161      1273              8          124  557  807  259     33
#tumour7     131       269              7          183  602  133  383     29
write.xlsx(sample.table,file = "EC_subcluster_sample_propotion.xlsx",rowNames=F,colNames=T)

#data for Fig3A
tumour.type.table<-data.frame(table(pure.EC.Seurat@meta.data$tumour.type,pure.EC.Seurat@meta.data$cluster.rev))
colnames(tumour.type.table)<-c("TumourType","CellType","Number")
tumour.type.table <- dcast(data = tumour.type.table,TumourType~CellType)  
#Artery Capillary Inflamatory EC Large Artery  TEC1  TEC2  Vein Venule
#Control               3643     19227            507         2083    77   283   848   3142
#Higher-grade glioma    432      1992             17          445  1508  1062   701    135
#Lower-grade glioma    2953     11269            176         1946   151   308   844   3575
write.xlsx(tumour.type.table,file = "EC_subcluster_tumour_propotion.xlsx",rowNames=F,colNames=T)


table(pure.EC.Seurat@meta.data$stim,pure.EC.Seurat@meta.data$cluster.rev)
#Artery Capillary Inflamatory EC Large Artery  TEC1  TEC2  Vein Venule
#Control   3643     19227            507         2083    77   283   848   3142
#Tumour    3385     13261            193         2391  1659  1370  1545   3710

#Fig2A
Idents(pure.EC.Seurat)<-"cluster.rev"

pdf("EC_subcluster_revised_UMAP.pdf",width = 6,height = 4)

cluster.color<-c("#86AC34","#E77D70","#54BA6F","#C59733","#927CB5","#55BBC2","#539BD3","#C56FA2")
DimPlot(pure.EC.Seurat,label = F,pt.size = 0.1,cols = cluster.color)

dev.off()

pdf(file = "EC_subcluster_UMAP_by_tumour_type.pdf",width = 6,height = 6)

tumour.type<-data.frame(type=pure.EC.Seurat@meta.data$tumour.type,
                        row.names = rownames(pure.EC.Seurat@meta.data))
control.id<-rownames(subset(tumour.type,tumour.type$type=="Control"))
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
DimPlot(pure.EC.Seurat,cells = control.id,group.by = "tumour.type",cols = "#619CFF",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Control")

Lower.grade.glioma.id<-rownames(subset(tumour.type,tumour.type$type=="Lower-grade glioma"))
DimPlot(pure.EC.Seurat,cells = Lower.grade.glioma.id,group.by = "tumour.type",cols = "#00BA38",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Lower-grade Gioma")

Higher.grade.glioma.id<-rownames(subset(tumour.type,tumour.type$type=="Higher-grade glioma"))
DimPlot(pure.EC.Seurat,cells = Higher.grade.glioma.id,group.by = "tumour.type",cols = "#F8766D",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Higher-grade Glioma")

dev.off()

#Find markers of all pure EC subcluster
pure.EC.markers<-FindAllMarkers(pure.EC.Seurat,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

top100.EC.markers<-pure.EC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_logFC)
write.xlsx(top100.EC.markers,file = "top100_EC_subcluster_markers.xlsx",rowNames=T,colNames=T)

#Fig2B
g.select<-c("MYO1B","APCDD1","CA2","VWA1","COL4A1","ESM1","PLVAP","CX3CL1","CXCL1","VCAM1",
            "ACKR1","CXCL11","ADGRG6","TSHZ2","IL1R1","IGFBP3","MFSD2A","CAVIN2","CA4","SLCO1A2",
            "ALPL","GLUL","SEMA3G","KRT19","CAV1","ELN","GJA5")
Idents(pure.EC.Seurat)<-factor(Idents(pure.EC.Seurat),
                               levels = c("TEC2","TEC1","Inflamatory EC","Vein","Venule","Capillary","Artery","Large Artery"))
DotPlot(pure.EC.Seurat,features = g.select,cols = "RdYlBu")+RotatedAxis()

#Fig3B
g.select<-c("COL4A1","COL4A2","HSPG2","ESM1","ANGPT2","APLN","SPARC","CD93","INSR","PLVAP","MFSD2A","APCDD1","CLDN5")
Idents(pure.EC.Seurat)<-factor(Idents(pure.EC.Seurat),
                        levels = c("Large Artery","Artery","Capillary","Venule","Vein","Inflamatory EC","TEC1","TEC2"))
DotPlot(pure.EC.Seurat,features = g.select,cols = "RdYlBu")+RotatedAxis()+coord_flip()

save(pure.EC.Seurat,file = "EC_rpca_subcluster_pc15_r0.8_clustered_filtered_revised.RData")


pure.EC.id<-colnames(pure.EC.Seurat)
write.table(pure.EC.id,file = "pure.EC.id.txt",quote = F,sep = "\t",row.names = F,col.names = F)



#############################################################################################################
############################################ Mural cell subcluster ##########################################
#############################################################################################################
mural.seurat<-subset(dat.seurat.integrated.final,idents = c("Mural cell"))
save(mural.seurat,file = "Mural_all_seurat.RData")

DefaultAssay(mural.seurat)<-"RNA"

mural.seurat = NormalizeData(mural.seurat, scale.factor =10000)
mural.seurat=FindVariableFeatures(mural.seurat, selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(mural.seurat)
mural.seurat <- ScaleData(mural.seurat, features = all.genes)

mural.seurat <- RunPCA(mural.seurat, features = VariableFeatures(object = mural.seurat))

mural.seurat <- JackStraw(mural.seurat, num.replicate = 100)

ElbowPlot(mural.seurat, ndims = 50)
mural.seurat <- FindNeighbors(mural.seurat, reduction = "pca", dims = 1:30)
mural.seurat <- FindClusters(mural.seurat, resolution = 0.5)
mural.seurat <- RunUMAP(mural.seurat, reduction = "pca", dims = 1:30)

DimPlot(mural.seurat , reduction = "umap", label = T)

library(future)
library(future.apply)
plan("multiprocess", workers = 3)
options(future.globals.maxSize = 20000 * 1024^2)

mural.seurat.list <- SplitObject(mural.seurat, split.by = "genotype")
rm(list=setdiff(ls(), "mural.seurat.list"))

mural.seurat.list <- future_lapply(X = mural.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = mural.seurat.list)
mural.seurat.list <- future_lapply(X = mural.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
save(mural.seurat.list,file = "Mural.seurat.list.RData")
anchors <- FindIntegrationAnchors(object.list = mural.seurat.list, reduction = "rpca", 
                                  dims = 1:30)
save(anchors, file="Mural.anchors.RData")

rm(list=setdiff(ls(), "anchors"))
plan("multiprocess", workers = 1)
Mural.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
dim(Mural.seurat.integrated)
#[1]  2000 37223
save(Mural.seurat.integrated, file="Mural.seurat.raw.integrated.RData")
rm(list=setdiff(ls(), "Mural.seurat.integrated"))

Mural.seurat.integrated<-ScaleData(Mural.seurat.integrated,verbose=F)
Mural.seurat.integrated<-RunPCA(Mural.seurat.integrated,verbose=F)
Mural.seurat.integrated<-RunUMAP(Mural.seurat.integrated,dims=1:20)

Mural.seurat.integrated<-FindNeighbors(Mural.seurat.integrated,reduction="pca",dims=1:20)
Mural.seurat.integrated<-FindClusters(Mural.seurat.integrated,resolution=0.8)
pdf("Mural_subcluster_UMAP.pdf",width = 8,height = 6)
DimPlot(Mural.seurat.integrated,label = T)
stim.color<-c("#619CFF","#F8766D")
DimPlot(Mural.seurat.integrated,label = F,group.by = "stim",cols = stim.color)
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
DimPlot(Mural.seurat.integrated,label = F,group.by = "tumour.type",cols = tumour.type.color)
DimPlot(Mural.seurat.integrated,label = F,split.by = "genotype")

dev.off()

#Mural subcluster QC
pdf("Mural_subcluster_QC.pdf",width = 8,height = 4)
VlnPlot(Mural.seurat.integrated,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
dev.off()

#Find markers of each mural subcluster
all.markers<-FindAllMarkers(Mural.seurat.integrated,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
write.xlsx(all.markers,file = "all.Mural.markers.xlsx",rowNames=T,colNames=T)

save(Mural.seurat.integrated,file = "Mural_rpca_subcluster_pc20_r0.8_clustered.RData")

#filter low quality and contamination clusters
pure.mural.seurat<-subset(Mural.seurat.integrated,idents = c(0:3,7,9:13))
DimPlot(pure.mural.seurat)   

pure.mural.cluster.rev=as.numeric(as.vector(pure.mural.seurat@meta.data$seurat_clusters))
table(pure.mural.cluster.rev)
#pure.mural.cluster.rev
#0    1    2    3    7    9   10   11   12   13 
#6477 5735 4066 3905 1918 1457 1362 1195  821  767

#merge cluster 0 and 1 into one cluster
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(0)]="1"
pure.mural.cluster.rev<-factor(pure.mural.cluster.rev,levels=c("1","2","3","7","9","10","11","12","13"))


pure.mural.seurat@meta.data[,"cluster.new"]<-pure.mural.cluster.rev
Idents(pure.mural.seurat)<-"cluster.new"
DimPlot(pure.mural.seurat)

#celltype annotation
pure.mural.cluster.rev=as.numeric(as.vector(pure.mural.seurat@meta.data$cluster.new))
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(1)]="Pericyte 1"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(2)]="aSMC 1"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(3)]="vSMC"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(7)]="Tumour pericyte"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(9)]="aSMC 1"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(10)]="aSMC 2"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(11)]="aSMC 3"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(12)]="Pericyte 1"
pure.mural.cluster.rev[pure.mural.cluster.rev %in% c(13)]="Pericyte 2"
pure.mural.cluster.rev<-factor(pure.mural.cluster.rev,
                               levels=c("Pericyte 1","aSMC 1","vSMC","Tumour pericyte","aSMC 2","aSMC 3",
                                                                          "Pericyte 2"))
                                  
                                  
pure.mural.seurat@meta.data[,"cluster.rev"]<-pure.mural.cluster.rev
Idents(pure.mural.seurat)<-"cluster.rev"
#Fig4A
pdf("Mrual_final_UMAP.pdf",width = 6,height = 4)
DimPlot(pure.mural.seurat)
dev.off()

pdf(file = "Mural_subcluster_UMAP_by_tumour_type.pdf",width = 6,height = 6)

tumour.type<-data.frame(type=pure.mural.seurat@meta.data$tumour.type,
                        row.names = rownames(pure.mural.seurat@meta.data))
control.id<-rownames(subset(tumour.type,tumour.type$type=="Control"))
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
DimPlot(pure.mural.seurat,cells = control.id,group.by = "tumour.type",cols = "#619CFF",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Control")

Lower.grade.glioma.id<-rownames(subset(tumour.type,tumour.type$type=="Lower-grade glioma"))
DimPlot(pure.mural.seurat,cells = Lower.grade.glioma.id,group.by = "tumour.type",cols = "#00BA38",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Lower-grade Gioma")

Higher.grade.glioma.id<-rownames(subset(tumour.type,tumour.type$type=="Higher-grade glioma"))
DimPlot(pure.mural.seurat,cells = Higher.grade.glioma.id,group.by = "tumour.type",cols = "#F8766D",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Higher-grade Glioma")

dev.off()

save(pure.mural.seurat,file = "Mural_rpca_subcluster_pc20_r0.8_clustered.filtered.revised.RData")
 
#Find markers of each pure mural subtype                                 
all.filtered.markers<-FindAllMarkers(pure.mural.seurat,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)
                                  
top100.mural.filtered.markers<-all.filtered.markers %>%
                                  group_by(cluster) %>%
                                  slice_max(n = 100, order_by = avg_logFC)
write.xlsx(top100.mural.filtered.markers,file = "top100.mural.filtered.markers.xlsx",rowNames=T,colNames=T)

#Fig4B
g.select<-c("NOTCH3","PDGFRB","RGS5","AGLN","ACTA2","IGFBP5","CCL19","CCL26","MT1M","MT1X","IL6","MT1A",
            "PLN","DSTN","MYH11","MYL9","GPX3","WNT6","SLC6A12","SLC6A13","COLEC12","ABCC9","COL4A1",
            "COL18A1","CPL1A1","ANGPT2","MMP9")

Idents(pure.mural.seurat)<-factor(Idents(pure.mural.seurat),
                                  levels = c("aSMC 3","aSMC 2","aSMC 1","vSMC","Pericyte 1","Pericyte 2","Tumour pericyte"))
DotPlot(pure.mural.seurat,features = g.select,cols = "RdYlBu")+RotatedAxis()+coord_flip()

                                 
#############################################################################################################                                 
#################################### Peri-vascular fibroblast subcluster ####################################
#############################################################################################################
peri.fibo.seurat<-subset(dat.seurat.integrated.final,idents = c("Perivascular Fibroblast"))
save(peri.fibo.seurat,file = "Perivascular_fibroblast_all_seurat.RData")

DefaultAssay(peri.fibo.seurat)<-"RNA"

peri.fibo.seurat = NormalizeData(peri.fibo.seurat , scale.factor =10000)
peri.fibo.seurat =FindVariableFeatures(peri.fibo.seurat , selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(peri.fibo.seurat )
peri.fibo.seurat  <- ScaleData(peri.fibo.seurat , features = all.genes)

peri.fibo.seurat  <- RunPCA(peri.fibo.seurat , features = VariableFeatures(object = peri.fibo.seurat ))

peri.fibo.seurat  <- FindNeighbors(peri.fibo.seurat , reduction = "pca", dims = 1:20)

peri.fibo.seurat  <- RunUMAP(peri.fibo.seurat , reduction = "pca", dims = 1:20)
peri.fibo.seurat  <- FindClusters(peri.fibo.seurat , resolution = 0.3)

DimPlot(peri.fibo.seurat , reduction = "umap", label = T)

library(future)
library(future.apply)
plan("multiprocess", workers = 3)
options(future.globals.maxSize = 20000 * 1024^2)

peri.fibo.seurat.list <- SplitObject(peri.fibo.seurat, split.by = "genotype")
rm(list=setdiff(ls(), "peri.fibo.seurat.list"))

peri.fibo.seurat.list <- future_lapply(X = peri.fibo.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = peri.fibo.seurat.list)
peri.fibo.seurat.list <- future_lapply(X = peri.fibo.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
save(peri.fibo.seurat.list,file = "peri.fibo.seurat.list.RData")
anchors <- FindIntegrationAnchors(object.list = peri.fibo.seurat.list, reduction = "rpca", 
                                  dims = 1:30)
save(anchors, file="peri.fibo.anchors.RData")

peri.fibo.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
dim(peri.fibo.seurat.integrated)
#[1] 2000 2602
save(peri.fibo.seurat.integrated, file="peri.fibo.seurat.raw.integrated.RData")
rm(list=setdiff(ls(), "EC.seurat.integrated"))

peri.fibo.seurat.integrated<-ScaleData(peri.fibo.seurat.integrated,verbose=F)
peri.fibo.seurat.integrated<-RunPCA(peri.fibo.seurat.integrated,verbose=F)
peri.fibo.seurat.integrated<-RunUMAP(peri.fibo.seurat.integrated,dims=1:15)

peri.fibo.seurat.integrated<-FindNeighbors(peri.fibo.seurat.integrated,reduction="pca",dims=1:15)
peri.fibo.seurat.integrated<-FindClusters(peri.fibo.seurat.integrated,resolution=0.4)
Idents(peri.fibo.seurat.integrated)<-'integrated_snn_res.0.4'

tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
pdf("perivascular_fibroblast_subcluster_UMAP.pdf",width = 6,height = 4)
DimPlot(peri.fibo.seurat.integrated,label = T,pt.size = 0.4,)
DimPlot(peri.fibo.seurat.integrated,group.by = "stim",cols = c("#619CFF","#F8766D"),pt.size = 0.4)
DimPlot(peri.fibo.seurat.integrated,group.by = "tumour.type",pt.size = 0.4,cols = tumour.type.color)
dev.off()

#QC
pdf("peri_fibo_subcluster_QC.pdf",width = 6,height = 4)
VlnPlot(peri.fibo.seurat.integrated,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
dev.off()

save(peri.fibo.seurat.integrated,file = "perivascular_fibroblast_subcluster_pc15_r0.4.RData")

DefaultAssay(peri.fibo.seurat.integrated)<-"RNA"
all.integrated.markers<-FindAllMarkers(peri.fibo.seurat.integrated,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

top100.integrated.markers<-all.integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_logFC)

write.xlsx(top100.integrated.markers,file = "top100.peri.fibo.markers.data.xlsx",rowNames=T,colNames=T)
write.xlsx(all.integrated.markers,file = "all.peri.fibo.markers.xlsx",rowNames=T,colNames=T)

#Filter low quality and contamination clusters
pure.peri.fibo.seurat<-subset(peri.fibo.seurat.integrated,idents = c(0,2))

pure.patient.table<-data.frame(table(pure.peri.fibo.seurat@meta.data$genotype,pure.peri.fibo.seurat@meta.data$seurat_clusters))
colnames(pure.patient.table)<-c("Patient","CellType","Number")
pure.patient.table <- dcast(data = pure.patient.table,Patient~CellType)  

pure.sample.table<-data.frame(table(pure.peri.fibo.seurat@meta.data$orig.ident,pure.peri.fibo.seurat@meta.data$seurat_clusters))
colnames(pure.sample.table)<-c("Sample","CellType","Number")
pure.sample.table <- dcast(data = pure.sample.table,Sample~CellType) 

pure.tumour.type.table<-data.frame(table(pure.peri.fibo.seurat@meta.data$tumour.type,pure.peri.fibo.seurat@meta.data$seurat_clusters))
colnames(pure.tumour.type.table)<-c("TumourType","CellType","Number")
pure.tumour.type.table <- dcast(data = pure.tumour.type.table,TumourType~CellType) 

pure.peri.fibo.cell.fraction<-list("Patient"=pure.patient.table,
                                   "Sample"=pure.sample.table,
                                   "Tumour Type"=pure.tumour.type.table)
write.xlsx(pure.peri.fibo.cell.fraction,file = "pure_preivascular_fibroblast_cell_fraction.xlsx",rowNames=F,colNames=T)

save(pure.peri.fibo.seurat,file = "pure_peri_fibo_seurat.RData")

#############################################################################################################                                 
###################################### Meningeal fibroblast subcluster ######################################
#############################################################################################################
meningeal.fibo.seurat<-subset(dat.seurat.integrated.final,idents=c("Meningeal Fibroblast"))
dim(meningeal.fibo.seurat)
#[1] 26197  6189
DimPlot(meningeal.fibo.seurat)

save(meningeal.fibo.seurat,file = "meningeal_fibo_all_seurat.RData")


DefaultAssay(meningeal.fibo.seurat)<-"RNA"

meningeal.fibo.seurat = NormalizeData(meningeal.fibo.seurat , scale.factor =10000)
meningeal.fibo.seurat =FindVariableFeatures(meningeal.fibo.seurat , selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(meningeal.fibo.seurat )
meningeal.fibo.seurat  <- ScaleData(meningeal.fibo.seurat , features = all.genes)

meningeal.fibo.seurat  <- RunPCA(meningeal.fibo.seurat , features = VariableFeatures(object = meningeal.fibo.seurat ))

meningeal.fibo.seurat  <- FindNeighbors(meningeal.fibo.seurat , reduction = "pca", dims = 1:20)

meningeal.fibo.seurat  <- RunUMAP(meningeal.fibo.seurat , reduction = "pca", dims = 1:20)
meningeal.fibo.seurat  <- FindClusters(meningeal.fibo.seurat , resolution = 0.3)

DimPlot(meningeal.fibo.seurat , reduction = "umap", label = T)

library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)

meningeal.fibo.seurat.list <- SplitObject(meningeal.fibo.seurat, split.by = "genotype")
rm(list=setdiff(ls(), "meningeal.fibo.seurat.list"))

meningeal.fibo.seurat.list <- future_lapply(X = meningeal.fibo.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list =meningeal.fibo.seurat.list)
meningeal.fibo.seurat.list <- future_lapply(X = meningeal.fibo.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
save(meningeal.fibo.seurat.list,file = "meningeal.fibo.seurat.list.RData")
anchors <- FindIntegrationAnchors(object.list = meningeal.fibo.seurat.list, reduction = "rpca", 
                                  dims = 1:30)
save(anchors, file="meningeal.fibo.anchors.RData")

meningeal.fibo.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

save(meningeal.fibo.seurat.integrated, file="meningeal.fibo.seurat.raw.integrated.RData")

meningeal.fibo.seurat.integrated<-ScaleData(meningeal.fibo.seurat.integrated,verbose=F)
meningeal.fibo.seurat.integrated<-RunPCA(meningeal.fibo.seurat.integrated,verbose=F)
meningeal.fibo.seurat.integrated<-RunUMAP(meningeal.fibo.seurat.integrated,dims=1:30)

meningeal.fibo.seurat.integrated<-FindNeighbors(meningeal.fibo.seurat.integrated,reduction="pca",dims=1:30)
meningeal.fibo.seurat.integrated<-FindClusters(meningeal.fibo.seurat.integrated,resolution=0.7)

Idents(meningeal.fibo.seurat.integrated)<-'integrated_snn_res.0.7'
meningeal.fibo.seurat.integrated@meta.data$seurat_clusters<-Idents(meningeal.fibo.seurat.integrated)

tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
pdf("meningeal_fibroblast_subcluster_UMAP.pdf",width = 6,height = 6)
DimPlot(meningeal.fibo.seurat.integrated,label = T,pt.size = 0.4)
DimPlot(meningeal.fibo.seurat.integrated,group.by = "stim",cols = c("#619CFF","#F8766D"),pt.size = 0.4)
DimPlot(meningeal.fibo.seurat.integrated,group.by = "tumour.type",pt.size = 0.4,cols = tumour.type.color)
dev.off()

pdf("meningeal_fibo_subcluster_QC.pdf",width = 6,height = 4)
VlnPlot(meningeal.fibo.seurat.integrated,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
dev.off()

save(meningeal.fibo.seurat.integrated,file = "meningeal_fibroblast_subcluster_pc30_r0.7.RData")

DefaultAssay(meningeal.fibo.seurat.integrated)<-"RNA"
all.integrated.markers<-FindAllMarkers(meningeal.fibo.seurat.integrated,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

top100.integrated.markers<-all.integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_logFC)

write.xlsx(top100.integrated.markers,file = "top100.meningeal.fibo.markers.data.xlsx",rowNames=T,colNames=T)
write.xlsx(all.integrated.markers,file = "all.meningeal.fibo.markers.xlsx",rowNames=T,colNames=T)

#Filter low quality and contamination clusters
pure.meningeal.fibo.seurat<-subset(meningeal.fibo.seurat.integrated,idents = c(0,2,3,4,9,12))
save(pure.meningeal.fibo.seurat,file = "pure_meningeal_fibo_seurat.RData")
DimPlot(pure.meningeal.fibo.seurat)

UMAP.table<-data.frame(pure.meningeal.fibo.seurat@reductions$umap@cell.embeddings)
UMAP.table.f<-subset(UMAP.table,UMAP.table$UMAP_1>0)

pure.meningeal.fibo.seurat.f<-subset(pure.meningeal.fibo.seurat,cells = rownames(UMAP.table.f))
DimPlot(pure.meningeal.fibo.seurat.f)

pdf("meningeal_fibroblast_filtered_subcluster_UMAP.pdf",width = 6,height = 6)
DimPlot(pure.meningeal.fibo.seurat.f,label = F,pt.size = 0.2)+
  scale_x_continuous(limits = c(0,6.5))+scale_y_continuous(limits = c(-6,5))
stim.color<-c("#619CFF","#F8766D")
DimPlot(pure.meningeal.fibo.seurat.f,group.by = "stim",pt.size = 0.2,cols = stim.color)+
  scale_x_continuous(limits = c(0,6.5))+scale_y_continuous(limits = c(-6,5))
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
DimPlot(pure.meningeal.fibo.seurat.f,group.by = "tumour.type",pt.size=0.2,cols = tumour.type.color)+
  scale_x_continuous(limits = c(0,6.5))+scale_y_continuous(limits = c(-6,5))
dev.off()

pdf(file = "meningeal_fibroblast_filtered_subcluster_UMAP_by_tumour_type.pdf",width = 6,height = 6)
tumour.type<-data.frame(type=pure.meningeal.fibo.seurat.f@meta.data$tumour.type,
                        row.names = rownames(pure.meningeal.fibo.seurat.f@meta.data))

save(pure.meningeal.fibo.seurat.f,file = "pure_meningeal_fibo_seurat_f.RData")



#############################################################################################################                                 
################################### Oligodendrocyte fibroblast subcluster ###################################
#############################################################################################################
Oligo.seurat<-subset(dat.seurat.integrated.final,idents=c("Oligodendrocyte"))
dim(Oligo.seurat)
#[1] 26197 17287
DimPlot(Oligo.seurat)

save(Oligo.seurat,file = "Oligodendrocyte_all_seurat.RData")


DefaultAssay(Oligo.seurat)<-"RNA"

Oligo.seurat = NormalizeData(Oligo.seurat , scale.factor =10000)
Oligo.seurat =FindVariableFeatures(Oligo.seurat , selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(Oligo.seurat )
Oligo.seurat  <- ScaleData(Oligo.seurat , features = all.genes)

Oligo.seurat  <- RunPCA(Oligo.seurat , features = VariableFeatures(object = Oligo.seurat ))

Oligo.seurat  <- FindNeighbors(Oligo.seurat , reduction = "pca", dims = 1:20)

Oligo.seurat  <- RunUMAP(Oligo.seurat , reduction = "pca", dims = 1:20)
Oligo.seurat  <- FindClusters(Oligo.seurat , resolution = 0.3)

DimPlot(Oligo.seurat , reduction = "umap", label = T)

library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 20000 * 1024^2)

Oligo.seurat.list <- SplitObject(Oligo.seurat, split.by = "genotype")
rm(list=setdiff(ls(), "Oligo.seurat.list"))

Oligo.seurat.list <- future_lapply(X = Oligo.seurat.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = Oligo.seurat.list)
Oligo.seurat.list <- future_lapply(X = Oligo.seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
save(Oligo.seurat.list,file = "Oligo.seurat.list.RData")
anchors <- FindIntegrationAnchors(object.list = Oligo.seurat.list, reduction = "rpca", 
                                  dims = 1:30)
save(anchors, file="Oligo.anchors.RData")

Oligo.seurat.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

save(Oligo.seurat.integrated, file="Oligo.seurat.raw.integrated.RData")

Oligo.seurat.integrated<-ScaleData(Oligo.seurat.integrated,verbose=F)
Oligo.seurat.integrated<-RunPCA(Oligo.seurat.integrated,verbose=F)
Oligo.seurat.integrated<-RunUMAP(Oligo.seurat.integrated,dims=1:30)

Oligo.seurat.integrated<-FindNeighbors(Oligo.seurat.integrated,reduction="pca",dims=1:30)
Oligo.seurat.integrated<-FindClusters(Oligo.seurat.integrated,resolution=0.5)
Idents(Oligo.seurat.integrated)<-'integrated_snn_res.0.5'

tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
pdf("Oligodendrocyte_subcluster_UMAP.pdf",width = 6,height = 6)
DimPlot(Oligo.seurat.integrated,label = T,pt.size = 0.4,)
DimPlot(Oligo.seurat.integrated,group.by = "stim",cols = c("#619CFF","#F8766D"),pt.size = 0.4)
DimPlot(Oligo.seurat.integrated,group.by = "tumour.type",pt.size = 0.4,cols = tumour.type.color)
DimPlot(Oligo.seurat.integrated,group.by = "genotype",pt.size = 0.4)
dev.off()

pdf("Oligodendrocyte_subcluster_QC.pdf",width = 6,height = 4)
VlnPlot(Oligo.seurat.integrated,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
dev.off()

save(Oligo.seurat.integrated,file = "Oligodendrocyte_subcluster_pc30_r0.5.RData")

DefaultAssay(Oligo.seurat.integrated)<-"RNA"
all.integrated.markers<-FindAllMarkers(Oligo.seurat.integrated,only.pos = T,min.pct = 0.25,logfc.threshold = 0.25)

nk=length(unique(all.integrated.markers$cluster))

all.integrated.marker.genes=sapply(0:(nk-1), function(i){
  t.k=all.integrated.markers[all.integrated.markers$cluster==i,]
  #t.k.genes=with(t.k, gene[avg_logFC >1 & p_val_adj <0.05 ])
  t.k.genes=with(t.k, gene[p_val_adj <0.05 ])
  t.k.genes=head(t.k.genes, 100)
  paste(t.k.genes, collapse = " ")
})
all.integrated.marker.genes=data.frame(cluster=0:(nk-1), genes=all.integrated.marker.genes)

top100.integrated.markers<-all.integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_logFC)
top50.integrated.markers<-all.integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_logFC)
write.xlsx(all.integrated.marker.genes,file = "top100.oligo.markers.xlsx",rowNames=T,colNames=T)
write.xlsx(top100.integrated.markers,file = "top100.oligo.markers.data.xlsx",rowNames=T,colNames=T)
write.xlsx(all.integrated.markers,file = "all.oligo.markers.xlsx",rowNames=T,colNames=T)

pure.oligo.seurat<-subset(Oligo.seurat.integrated,idents=c(0,1,3,5,8,11))
DimPlot(pure.oligo.seurat)
save(pure.oligo.seurat,file = "pure_oligo_seurat.RData")

pdf("Oligodendrocyte_filtered_subcluster_UMAP.pdf",width = 6,height = 6)
DimPlot(pure.oligo.seurat,label = F,pt.size = 0.1)+
  scale_x_continuous(limits = c(-7,10))+scale_y_continuous(limits = c(-6,5))
stim.color<-c("#619CFF","#F8766D")
DimPlot(pure.oligo.seurat,group.by = "stim",pt.size = 0.1,cols = stim.color)+
  scale_x_continuous(limits = c(-7,10))+scale_y_continuous(limits = c(-6,5))
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")
DimPlot(pure.oligo.seurat,group.by = "tumour.type",pt.size=0.1,cols = tumour.type.color)+
  scale_x_continuous(limits = c(-7,10))+scale_y_continuous(limits = c(-6,5))
dev.off()



#############################################################################################################                                 
######################################### Proliferating subcluster ##########################################
#############################################################################################################
prolife.seurat<-subset(dat.seurat.integrated.final,idents=c("Proliferating cell"))
dim(prolife.seurat)
#[1] 26197   994
DimPlot(prolife.seurat)

save(prolife.seurat,file = "Prolife_all_seurat.RData")


DefaultAssay(prolife.seurat)<-"RNA"

prolife.seurat = NormalizeData(prolife.seurat , scale.factor =10000)
prolife.seurat =FindVariableFeatures(prolife.seurat , selection.method = "vst", nfeatures = 2000)

# Run the standard workflow for visualization and clustering
all.genes <- rownames(prolife.seurat )
prolife.seurat  <- ScaleData(prolife.seurat , features = all.genes)

prolife.seurat  <- RunPCA(prolife.seurat , features = VariableFeatures(object =prolife.seurat ))

prolife.seurat  <- FindNeighbors(prolife.seurat , reduction = "pca", dims = 1:20)

prolife.seurat  <- RunUMAP(prolife.seurat , reduction = "pca", dims = 1:20)
prolife.seurat  <- FindClusters(prolife.seurat , resolution = 0.3)

pdf("Proliferating_cell_subcluster_no_align.pdf",width = 6,height = 6)
DimPlot(prolife.seurat , reduction = "umap", label = T)
DimPlot(prolife.seurat , reduction = "umap", label = F,group.by = "stim")
DimPlot(prolife.seurat , reduction = "umap", label = F,group.by = "genotype")
DimPlot(prolife.seurat , reduction = "umap", label = F,group.by = "tumour.type")

VlnPlot(prolife.seurat,features = c("nCount_RNA","nFeature_RNA"),pt.size = 0)
dev.off()

pure.prolife.seurat<-subset(prolife.seurat,idents = c(0,1,3,4,5))
dim(pure.prolife.seurat)
#[1] 26197   764
DimPlot(pure.prolife.seurat)

save(pure.prolife.seurat,file = "pure.prolife.seurat.no.align.RData")
pure.prolife.id<-colnames(pure.prolife.seurat)
write.table(pure.prolife.id,file = "pure.prolife.id.txt",quote = F,sep = "\t",row.names = F,col.names = F)



#############################################################################################################                                 
########################################### Macrophage subcluster ###########################################
#############################################################################################################
macro.seurat<-subset(dat.seurat.integrated.final,idents=c("Macrophage"))
dim(macro.seurat)
#[1] 26197  1913
DimPlot(macro.seurat)

save(macro.seurat,file="Macrophage_all_seurat.RData")

DefaultAssay(macro.seurat)<-"RNA"

macro.seurat<-NormalizeData(macro.seurat, scale.factor = 10000)
macro.seurat<-FindVariableFeatures(macro.seurat, selection.method = "vst", nfeatures = 2000)

all.genes<-rownames(macro.seurat)
macro.seurat<-ScaleData(macro.seurat, features = all.genes)

macro.seurat<-RunPCA(macro.seurat,features = VariableFeatures(object = macro.seurat))

macro.seurat<-FindNeighbors(macro.seurat, reduction = "pca", dims = 1:20)

macro.seurat<-RunUMAP(macro.seurat, reduction = "pca", dims = 1:20)
macro.seurat<-FindClusters(macro.seurat, resolution = 0.3)

DimPlot(macro.seurat)

pure.macro.seurat<-subset(macro.seurat, idents = c(2,3,6))
dim(pure.macro.seurat)
#[1] 26197   678

DimPlot(pure.macro.seurat, label = T)

pure.macro.seurat.umap<-data.frame(pure.macro.seurat@reductions$umap@cell.embeddings)
pure.macro.seurat.umap.f<-subset(pure.macro.seurat.umap, pure.macro.seurat.umap$UMAP_1>0)

pure.macro.seurat.f<-subset(pure.macro.seurat, cells = rownames(pure.macro.seurat.umap.f))
dim(pure.macro.seurat.f)
#[1] 26197   675

save(pure.macro.seurat.f, file = "pure.macro.seurat.RData")



#############################################################################################################                                 
############################################# T cell subcluster #############################################
#############################################################################################################
T.seurat<-subset(dat.seurat.integrated.final, idents = "T cell")
dim(T.seurat)
#[1] 26197   706

save(T.seurat, file = "T_all_seurat.RData")

DefaultAssay(T.seurat)<-"RNA"

T.seurat<-NormalizeData(T.seurat, scale.factor = 10000)
T.seurat<-FindVariableFeatures(T.seurat, selection.method = "vst", nfeatures = 2000)

all.genes<-rownames(T.seurat)
T.seurat<-ScaleData(T.seurat, features = all.genes)

T.seurat<-RunPCA(T.seurat, features = VariableFeatures(T.seurat))

T.seurat<-FindNeighbors(T.seurat, reduction = "pca", dims = 1:20)

T.seurat<-RunUMAP(T.seurat, reduction = "pca", dims = 1:20)
T.seurat<-FindClusters(T.seurat, resolution = 0.3)

DimPlot(T.seurat)
VlnPlot(T.seurat, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0)

save(T.seurat, file = "Tcell.subcluster.no.aglin.RData")

pure.T.seurat<-subset(T.seurat, idents = c(0,1,4))
dim(pure.T.seurat)
#[1] 26197   415

DimPlot(pure.T.seurat)

pure.T.seurat.umap<-data.frame(pure.T.seurat@reductions$umap@cell.embeddings)
pure.T.seurat.umap.f<-subset(pure.T.seurat.umap, pure.T.seurat.umap$UMAP_1<0)

pure.T.seurat.f<-subset(pure.T.seurat, cells = rownames(pure.T.seurat.umap.f))
dim(pure.T.seurat.f)
#[1] 26197   410
DimPlot(pure.T.seurat.f)

save(pure.T.seurat.f, file = "pure.T.seurat.no.aglin.RData")



#############################################################################################################                                 
############################################# Astrocyte subcluster ##########################################
#############################################################################################################
astro.seurat<-subset(dat.seurat.integrated.final, idents = "Astrocyte")
dim(astro.seurat)
#[1] 26197  1698

DefaultAssay(astro.seurat)<-"RNA"

astro.seurat<-NormalizeData(astro.seurat, scale.factor = 10000)
astro.seurat<-FindVariableFeatures(astro.seurat,  selection.method = "vst", nfeatures = 2000)

all.genes<-rownames(astro.seurat)
astro.seurat<-ScaleData(astro.seurat, features = all.genes)

astro.seurat<-RunPCA(astro.seurat, features = VariableFeatures(astro.seurat))
astro.seurat<-FindNeighbors(astro.seurat, reduction = "pca", dims = 1:20)

astro.seurat<-RunUMAP(astro.seurat, reduction = "pca", dims = 1:20)
astro.seurat<-FindClusters(astro.seurat, resolution = 0.3)

DimPlot(astro.seurat)

pure.astro.seurat<-subset(astro.seurat, idents = 5)
dim(pure.astro.seurat)
#[1] 26197   153
DimPlot(pure.astro.seurat)

pure.astro.seurat.umap<-data.frame(pure.astro.seurat@reductions$umap@cell.embeddings)
pure.astro.seurat.umap.f<-subset(pure.astro.seurat.umap, pure.astro.seurat.umap$UMAP_1>-5)

pure.astro.seurat.f<-subset(pure.astro.seurat, cells = rownames(pure.astro.seurat.umap.f))
dim(pure.astro.seurat.f)
#[1] 26197   150
DimPlot(pure.astro.seurat.f)

save(pure.astro.seurat.f, file = "pure.Astro.no.aglin.RData")



#############################################################################################################                                 
####################################### final umap of all pure cells ########################################
#############################################################################################################
all.pure.cell.id<-c(colnames(pure.EC.Seurat), colnames(pure.mural.seurat), colnames(pure.peri.fibo.seurat),
                    colnames(pure.meningeal.fibo.seurat.f), colnames(pure.oligo.seurat), colnames(pure.prolife.seurat),
                    colnames(pure.macro.seurat.f), colnames(tumor.seurat), colnames(pure.T.seurat.f),
                    colnames(pure.astro.seurat.f))

dat.seurat.integrated.final.pure<-subset(dat.seurat.integrated.final, cells = all.pure.cell.id)
dim(dat.seurat.integrated.final.pure)
#[1] 26197 103230

save(dat.seurat.integrated.final.pure, file = "all_pure_cell_seurat.RData")

#Fig1B
pdf("All_cells_UMAP_Fig1B.pdf",width = 6,height = 6)
DimPlot(dat.seurat.integrated.final.pure)
dev.off()

#Fig1C
pdf(file = "All_cells_cluster_UMAP_by_tumour_type.pdf",width = 6,height = 6)

tumour.type<-data.frame(type=pure.all.seurat@meta.data$tumour.type,
                        row.names = rownames(pure.all.seurat@meta.data))
tumour.type.color<-c("#619CFF","#F8766D","#00BA38")

control.id<-rownames(subset(tumour.type,tumour.type$type=="Control"))
DimPlot(pure.all.seurat,cells = control.id,group.by = "tumour.type",cols = "#619CFF",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Control")

Lower.grade.glioma.id<-rownames(subset(tumour.type,tumour.type$type=="Lower-grade glioma"))
DimPlot(pure.all.seurat,cells = Lower.grade.glioma.id,group.by = "tumour.type",cols = "#00BA38",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Lower-grade Gioma")

Higher.grade.glioma.id<-rownames(subset(tumour.type,tumour.type$type=="Higher-grade glioma"))
DimPlot(pure.all.seurat,cells = Higher.grade.glioma.id,group.by = "tumour.type",cols = "#F8766D",pt.size = 0.1)+
  theme(legend.position = "none")+labs(title = "Higher-grade Glioma")

dev.off()

#data for Fig1G
table(Idents(dat.seurat.integrated.final.pure))
#Endothelial              Mural cell         Oligodendrocyte    Meningeal Fibroblast 
#57324                   27703                   10862                    3781 
#Macrophage Perivascular Fibroblast               Astrocyte      Proliferating cell 
#675                     950                     150                     764 
#T cell       High Matrix cells 
#410                     611 

table(dat.seurat.integrated.final.pure@meta.data$genotype,dat.seurat.integrated.final.pure@meta.data$cluster.rev)
#Astrocyte Endothelial High Matrix cells Macrophage Meningeal Fibroblast Mural cell
#Sample027         0       10887                 0          4                   19       2044
#Sample036         0        4642                 0          2                  104       7564
#Sample037       104       14449                 0          3                 1177       6713
#Sample038         4        7750               575        649                   71       3457
#Sample043        41        5285                 0          7                 2292       3171
#Sample045         1        8306                 1          5                   92       2127
#Sample046         0        6005                35          5                   26       2627

Oligodendrocyte Perivascular Fibroblast Proliferating cell T cell
#Sample027            1361                      77                 18     15
#Sample036             685                     196                 47     14
#Sample037            2045                     305                 27      9
#Sample038            1146                     140                248    326
#Sample043            1079                      48                 19      2
#Sample045            2172                      69                268     11
#Sample046            2374                     115                137     33

table(dat.seurat.integrated.final.pure@meta.data$tumour.type,dat.seurat.integrated.final.pure@meta.data$cluster.rev)
#Astrocyte Endothelial High Matrix cells Macrophage Meningeal Fibroblast
#Control                     3       29810                 0         10                 1946
#Higher-grade glioma         5        6292               611        653                   28
#Lower-grade glioma        142       21222                 0         12                 1807

#Mural cell Oligodendrocyte Perivascular Fibroblast Proliferating cell T cell
#Control                  18590            5666                     417                184     28
#Higher-grade glioma       2431            2731                     168                533    363
#Lower-grade glioma        6682            2465                     365                 47     19

#Find all final celltype markers after filtering
all.main.celltype.markers<-FindAllMarkers(pure.all.seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.xlsx(all.main.celltype.markers,file="all.cluster.rev.markers.xlsx",rowNames=T,colNames=T)


#############################################################################################################                                 
################################################# Session Info ##############################################
#############################################################################################################
R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X  11.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] reshape2_1.4.4     patchwork_1.1.1    ggplot2_3.3.5      Seurat_3.1.1       scRNAtoolVis_0.0.3
  [6] clustree_0.4.4     future_1.22.1      future.apply_1.8.1




