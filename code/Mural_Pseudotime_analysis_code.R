library(Seurat)
library(monocle)
library(ggplot2)
library(openxlsx)
library(patchwork)

#load Mural cell Seurat object
load("Mural_rpca_subcluster_pc20_r0.8_clustered.filtered.revised.RData")
DimPlot(pure.mural.seurat)



#############################################################################################################
######################################## Control Mural Zonation #############################################
#############################################################################################################
Idents(pure.mural.seurat)<-"tumour.type"
Control.monocle.data<-subset(pure.mural.seurat,idents=c("Control"))
Idents(Control.monocle.data)<-"cluster.rev"

dim(Control.monocle.data)
[1] 26197 18590
#filtering cells
Control.monocle.data.f<-subset(Control.monocle.data,subset=nCount_RNA>2000&nFeature_RNA>1000)
Control.monocle.data.f<-subset(Control.monocle.data.f,idents=c("aSMC 1","aSMC 2","aSMC 3","Pericyte 1","Pericyte 2",
                                                               "vSMC"))
save(Control.monocle.data.f,file = "Mural_Control_monocle_seurat.RData")

#Find ordering genes
Idents(pure.mural.seurat)<-"cluster.rev"
aSMC.vs.vSMC.DEGs<-FindMarkers(Control.monocle.data.f,ident.1 = c("aSMC 1","aSMC 2","aSMC 3"),ident.2 = "vSMC",min.pct=0.25,logfc.threshold = 0.58)
aSMC.vs.pericyte.DEGs<-FindMarkers(Control.monocle.data.f,ident.1 = c("aSMC 1","aSMC 2","aSMC 3"),ident.2 = c("Pericyte 1","Pericyte 2"),
                                   min.pct=0.25,logfc.threshold = 0.58)
vSMC.vs.pericyte.DEGs<-FindMarkers(Control.monocle.data.f,ident.1 = c("vSMC"),ident.2 = c("Pericyte 1","Pericyte 2"),
                                   min.pct=0.25,logfc.threshold = 0.58)
order.genes<-c(rownames(aSMC.vs.vSMC.DEGs),rownames(aSMC.vs.pericyte.DEGs),rownames(vSMC.vs.pericyte.DEGs))

order.genes<-unique(order.genes)
order.genes.f=order.genes[-c(grep("RPL|RPS",order.genes))]

DimPlot(Control.monocle.data.f)
dim(Control.monocle.data.f)
[1] 26197 12574

#run monocle2 pipeline
control_ann<-Control.monocle.data.f@meta.data
control_ann$celltype<-Idents(Control.monocle.data.f)

gene_ann<-data.frame(gene_short_name=rownames(Control.monocle.data.f@assays$RNA),
                     row.names = rownames(Control.monocle.data.f@assays$RNA))

head(gene_ann)

pd<-new("AnnotatedDataFrame",data=control_ann)
fd<-new("AnnotatedDataFrame",data=gene_ann)
ce=as.data.frame(Control.monocle.data.f@assays$RNA@counts)

cds<-newCellDataSet(as.matrix(ce),
                    phenoData = pd,
                    featureData = fd,
                    expressionFamily = negbinomial.size(),
                    lowerDetectionLimit = 1)

cds<-detectGenes(cds,min_expr = 1)

cds<-cds[fData(cds)$num_cells_expressed>10,]

cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)


cds<-reduceDimension(cds,max_components = 2,num_dim=6,reduction_method = "tSNE",verbose = T)
cds<-clusterCells(cds,1,2,num_clusters=5)

pData(cds)$Cluster=pData(cds)$celltype

g.select<-read.xlsx("top100.mural.filtered.markers.xlsx",rowNames = F,colNames = T)
cds<-setOrderingFilter(cds,order.genes.f)
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,method="DDRTree")
cds<-orderCells(cds)

save(cds,file="Control_Mural_monocle_result.RData")

pdf("pseudotime_test_result.pdf",width = 11,height = 6)
p1<-plot_cell_trajectory(cds,color_by="Cluster",cell_size = 0.4,cols=cols,show_branch_points = F)+
  scale_color_manual(breaks = c("Pericyte 1", "Pericyte 2", "aSMC 1","aSMC 2","aSMC 3","vSMC"), 
                     values=c("#F8766D","#FA61D7","#C49A01","#01B6EC","#A58AFF","#52B400")) 
p2<-plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.4,show_branch_points = F)
p1|p2
dev.off()


library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
plotdf=pData(cds)

ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )


diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res_1<-diff_test_res[order(diff_test_res$qval),]
write.xlsx(diff_test_res_1,file = "Control_Mural_pseudotime_DEGs.xlsx",rowNames=T,colNames=T)
sig_gene_names <- as.character(diff_test_res_1$gene_short_name[1:1000])

newdata<-data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), 
                                     max(pData(cds)$Pseudotime),length.out = 12574))
genSmoothCurves_mat<-genSmoothCurves(cds[sig_gene_names,],
                                     new_data = newdata,
                                     cores = 10)

library(pheatmap)
pdf("Mural_contrl_zonation_with_cluster.pdf",width = 6,height = 8)
p<-pheatmap(log10(genSmoothCurves_mat+1),
            scale = "row",
            cluster_rows = T,clustering_distance_rows = "correlation",clustering_method = "ward.D2",
            cluster_cols = F,
            show_rownames = F,
            show_colnames = F,cutree_rows = 5)
print(p)
dev.off()


geneClust <- p$tree_row
geneCluster2Three<- cutree(geneClust,k=5)

newOrder=genSmoothCurves_mat[p$tree_row$order,]

geneCluster2Three.table<-data.frame(geneCluster2Three)
geneCluster2Three.table$Gene<-rownames(geneCluster2Three.table)
table(geneCluster2Three.table$geneCluster2Three)
1   2   3   4   5 
249 446  98 118  89 
geneCluster2Three.table<-geneCluster2Three.table[rownames(newOrder),]

#reordering zonation genes
cluster1.genes<-data.frame(Gene=as.character(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three=="2")$Gene),Cluster=rep("1",446))
cluster2.genes<-data.frame(Gene=as.character(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three=="5")$Gene),Cluster=rep("1",89))
cluster3.genes<-data.frame(Gene=as.character(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three=="1")$Gene),Cluster=rep("2",249))
cluster4.genes<-data.frame(Gene=as.character(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three=="4")$Gene),Cluster=rep("2",118))
cluster5.genes<-data.frame(Gene=as.character(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three=="3")$Gene),Cluster=rep("2",98))


zonation.cluster.genes<-rbind(cluster1.genes,cluster2.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster3.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster4.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster5.genes)

write.xlsx(zonation.cluster.genes,file = "Mural_cell_zonation_heatmap_genes.xlsx",rowNames=F,colNames=T)

#Fig4D zonation heatmap
newOrder<-data.frame(newOrder)
newOrder<-newOrder[intersect(zonation.cluster.genes$Gene,rownames(newOrder)),]
pdf("Control_EC_zonation.pdf",width = 6,height = 8)
p1<-pheatmap(log10(newOrder+1),
             scale = "row",
             cluster_rows = F,clustering_distance_rows = "euclidean",clustering_method = "ward.D2",
             cluster_cols = F,
             show_rownames = F,
             show_colnames = F)
print(p1)
dev.off()

cells.pseudo.table<-data.frame(cell=colnames(cds),
                               pesudotime=cds@phenoData@data$Pseudotime,
                               celltype=cds@phenoData@data$Cluster)
cells.pseudo.table<-cells.pseudo.table[order(cells.pseudo.table$pesudotime),]
color<-cells.pseudo.table$celltype
values=c("#F8766D","#FA61D7","#C49A01","#01B6EC","#A58AFF","#52B400")
color<-gsub("Pericyte 1","#F8766D",color)
color<-gsub("Pericyte 2","#FA61D7",color)
color<-gsub("aSMC 1","#C49A01",color)
color<-gsub("aSMC 2","#01B6EC",color)
color<-gsub("aSMC 3","#A58AFF",color)
color<-gsub("vSMC","#52B400",color)
cells.pseudo.table$color<-color

#barplot of each genes in pseudotime
dat.seurat.normalized=exp(GetAssayData(object = Control.monocle.data.f))-1
dat.seurat.normalized.f<-dat.seurat.normalized[sig_gene_names,]
dat.seurat.normalized.f=round(dat.seurat.normalized.f, 0)
dat.seurat.normalized.f<-data.frame(dat.seurat.normalized.f[,intersect(cells.pseudo.table$cell,colnames(dat.seurat.normalized.f))])
dim(dat.seurat.normalized.f)
[1]  1000 12574


#Fig4D barplot
g.select<-c("ACTA2","PLN","SLC16A2","SLC16A3")

pdf("Control_Mural_zonation_genes_barplot_with_color.PDF",8,14)
par(mfrow=c(5,2))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
for(i in 1:length(g.select)){
  t.g=g.select[i];
  p<-barplot(as.numeric((dat.seurat.normalized.f[t.g,])), col=color,main=t.g,border = NA,ylab="counts",space = 0.1)
  box()
  curve<-loess(as.numeric((dat.seurat.normalized.f[t.g,]))~c(1:ncol(dat.seurat.normalized.f)))
  lines(p,predict(curve),col="red")
}
dev.off()



#############################################################################################################
######################################## LGG Mural Zonation #################################################
#############################################################################################################
LGG.monocle.data<-subset(pure.mural.seurat,idents=c("Lower-grade glioma"))
Idents(LGG.monocle.data)<-"cluster.rev"

LGG.monocle.data.f<-subset(LGG.monocle.data,idents=c("aSMC 1","aSMC 2","aSMC 3","Pericyte 1","Pericyte 2","vSMC"))
DimPlot(LGG.monocle.data.f)
dim(LGG.monocle.data.f)
[1] 26197  6662

LGG_ann<-LGG.monocle.data.f@meta.data
LGG_ann$celltype<-Idents(LGG.monocle.data.f)

gene_ann<-data.frame(gene_short_name=rownames(LGG.monocle.data.f@assays$RNA),
                     row.names = rownames(LGG.monocle.data.f@assays$RNA))

head(gene_ann)

LGG.pd<-new("AnnotatedDataFrame", data=LGG_ann)
LGG.fd<-new("AnnotatedDataFrame", data=gene_ann)
LGG.ce=as.data.frame(LGG.monocle.data.f@assays$RNA@counts)

LGG.cds<-newCellDataSet(as.matrix(LGG.ce),
                    phenoData = LGG.pd,
                    featureData = LGG.fd,
                    expressionFamily = negbinomial.size(),
                    lowerDetectionLimit = 1)

LGG.cds<-detectGenes(LGG.cds, min_expr = 1)

LGG.cds<-cds[fData(LGG.cds)$num_cells_expressed>10,]

LGG.cds<-estimateSizeFactors(LGG.cds)
LGG.cds<-estimateDispersions(LGG.cds)

LGG.cds<-reduceDimension(LGG.cds,max_components = 2,num_dim=6,reduction_method = "tSNE",verbose = T)
LGG.cds<-clusterCells(LGG.cds,1,2,num_clusters = 5)

pData(LGG.cds)$cluster=pData(LGG.cds)$celltype

cds<-setOrderingFilter(LGG.cds, order.genes.f)
plot_ordering_genes(LGG.cds)

cds<-reduceDimension(LGG.cds,max_components = 2,method="DDRTree")
cds<-orderCells(LGG.cds)

save(cds,file="LGG_Mural_monocle_result.RData")

#Fig4H
pdf("LGG_Mural_pseudotime_result.pdf",width=11,height=6)

p1<-plot_cell_trajectory(LGG.cds,color_by="cluster",cols=cols,cell_size = 0.4,show_branch_points = F)+
  scale_color_manual(breaks = c("Pericyte 1","Pericyte 2","aSMC 1","aSMC 2","aSMC 3","vSMC"),
                     values = c("#F8766D","#FA61D7","#C49A01","#01B6EC","#A58AFF","#52B400"))
p2<-plot_cell_trajectory(LGG.cds,color_by = "Pseudotime",cell_size = 0.4,show_branch_points = F)
p1|p2

dev.off()



#############################################################################################################
######################################## GBM Mural Zonation #################################################
#############################################################################################################
GBM.monocle.data<-subset(pure.mural.seurat,idents=c("Higher-grade glioma"))
Idents(GBM.monocle.data)<-"cluster.rev"

dim(GBM.monocle.data)
[1] 26197  2431

GBM_ann<-GBM.monocle.data@meta.data
GBM_ann$celltype<-Idents(GBM.monocle.data)

gene_ann<-data.frame(gene_short_name=rownames(GBM.monocle.data@assays$RNA),
                     row.names = rownames(GBM.monocle.data@assays$RNA))

head(gene_ann)

GBM.pd<-new("AnnotatedDataFrame", data=GBM_ann)
GBM.fd<-new("AnnotatedDataFrame", data=gene_ann)
GBM.ce=as.data.frame(GBM.monocle.data@assays$RNA@counts)

GBM.cds<-newCellDataSet(as.matrix(GBM.ce),
                    phenoData = GBM.pd,
                    featureData = GBM.fd,
                    expressionFamily = negbinomial.size(),
                    lowerDetectionLimit = 1)

GBM.cds<-detectGenes(GBM.cds, min_expr = 1)

GBM.cds<-GBM.cds[fData(GBM.cds)$num_cells_expressed>10,]

GBM.cds<-estimateSizeFactors(GBM.cds)
GBM.cds<-estimateDispersions(GBM.cds)

GBM.cds<-reduceDimension(GBM.cds,max_components = 2,num_dim=6,reduction_method = "tSNE",verbose = T)
GBM.cds<-clusterCells(GBM.cds,1,2,num_clusters = 5)

pData(GBM.cds)$cluster=pData(GBM.cds)$celltype

GBM.cds<-setOrderingFilter(GBM.cds, order.genes.f)
plot_ordering_genes(GBM.cds)

GBM.cds<-reduceDimension(GBM.cds,max_components = 2,method="DDRTree")
GBM.cds<-orderCells(GBM.cds)

save(GBM.cds,file="GBM_Mural_monocle_result.RData")

#Fig4I
pdf("GBM_Mural_pseudotime_result.pdf",width=11,height=6)

p1<-plot_cell_trajectory(GBM.cds,color_by="cluster",cols=cols,cell_size = 0.4,show_branch_points = F)+
  scale_color_manual(breaks = c("Pericyte 1","Pericyte 2","aSMC 1","aSMC 2","aSMC 3","vSMC","Tumour pericyte"),
                     values = c("#F8766D","#FA61D7","#C49A01","#01B6EC","#A58AFF","#52B400","#02C095"))
p2<-plot_cell_trajectory(GBM.cds,color_by = "Pseudotime",cell_size = 0.4,show_branch_points = F)
p1|p2

dev.off()





