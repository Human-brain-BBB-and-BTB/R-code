setwd("/Volumes/Elements/7human_sample_202111/13_sample_seurat/endothelial_subcluster/all_rpca_pc50_feature2000/Cluster_pc15_r0.8/monocle2/Control")

library(Seurat)
library(monocle)
library(ggplot2)
library(openxlsx)
library(patchwork)

#load EC seurat data object
load("EC_rpca_subcluster_pc15_r0.8_clustered_filtered_revised.RData")

#############################################################################################################
######################################## Control EC Zonation ################################################
#############################################################################################################
##identtify marker genes using for cell ordring
large.artery.vs.vein.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Large Artery"),ident.2 = c("Vein"),min.pct = 0.25,logfc.threshold = 0.25)
large.artery.vs.venule.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Large Artery"),ident.2 = c("Venule"),min.pct = 0.25,logfc.threshold = 0.25)
large.artery.vs.cap.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Large Artery"),ident.2 = c("Capillary"),min.pct = 0.25,logfc.threshold = 0.25)
large.artery.vs.artery.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Large Artery"),ident.2 = c("Artery"),min.pct = 0.25,logfc.threshold = 0.25)

artery.vs.vein.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Artery"),ident.2 = c("Vein"),min.pct = 0.25,logfc.threshold = 0.25)
artery.vs.venule.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Artery"),ident.2 = c("Venule"),min.pct = 0.25,logfc.threshold = 0.25)
artery.vs.cap.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Artery"),ident.2 = c("Capillary"),min.pct = 0.25,logfc.threshold = 0.25)

vein.vs.cap.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Vein"),ident.2 = c("Capillary"),min.pct = 0.25,logfc.threshold = 0.25)
vein.vs.venule.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Vein"),ident.2 = c("Venule"),min.pct = 0.25,logfc.threshold = 0.25)
venule.vs.cap.DEGs<-FindMarkers(pure.EC.Seurat,ident.1 = c("Venule"),ident.2 = c("Capillary"),min.pct = 0.25,logfc.threshold = 0.25)

large.artery.vs.cap.DEGs.f<-subset(large.artery.vs.cap.DEGs,abs(large.artery.vs.cap.DEGs$avg_logFC)>0.58)
large.artery.vs.artery.DEGs.f<-subset(large.artery.vs.artery.DEGs,abs(large.artery.vs.artery.DEGs$avg_logFC)>0.58)
large.artery.vs.vein.DEGs.f<-subset(large.artery.vs.vein.DEGs,abs(large.artery.vs.vein.DEGs$avg_logFC)>0.58)
large.artery.vs.venule.DEGs.f<-subset(large.artery.vs.venule.DEGs,abs(large.artery.vs.venule.DEGs$avg_logFC)>0.58)

artery.vs.cap.DEGs.f<-subset(artery.vs.cap.DEGs,abs(artery.vs.cap.DEGs$avg_logFC)>0.58)
artery.vs.vein.DEGs.f<-subset(artery.vs.vein.DEGs,abs(artery.vs.vein.DEGs$avg_logFC)>0.58)
artery.vs.venule.DEGs.f<-subset(artery.vs.venule.DEGs,abs(artery.vs.venule.DEGs$avg_logFC)>0.58)

vein.vs.venule.DEGs.f<-subset(vein.vs.venule.DEGs,abs(vein.vs.venule.DEGs$avg_logFC)>0.58)
vein.vs.cap.DEGs.f<-subset(vein.vs.cap.DEGs,abs(vein.vs.cap.DEGs$avg_logFC)>0.58)
venule.vs.cap.DEGs.f<-subset(venule.vs.cap.DEGs,abs(venule.vs.cap.DEGs$avg_logFC)>0.58)

g.select<-c(rownames(large.artery.vs.cap.DEGs.f),rownames(large.artery.vs.artery.DEGs.f),rownames(large.artery.vs.vein.DEGs.f),rownames(large.artery.vs.venule.DEGs.f),
            rownames(artery.vs.cap.DEGs.f),rownames(artery.vs.vein.DEGs.f),rownames(artery.vs.venule.DEGs.f),
            rownames(vein.vs.venule.DEGs.f),rownames(vein.vs.cap.DEGs.f),rownames(venule.vs.cap.DEGs.f))
g.select<-data.frame(gene=unique(g.select))

##Run monocle2 pipeline
#remove Tumor ECs and inflamatory EC
EC.subset<-subset(pure.EC.Seurat,idents = c("Large Artery","Artery","Capillary","Venule","Vein"))
dim(EC.subset)
[1] 26197 53235

Idents(EC.subset)<-"tumour.type"
Control.monocle.data<-subset(EC.subset,idents = "Control")
dim(Control.monocle.data)
[1] 26197 28943

##filter cells which nCount<3000 & nfeature<1000
Idents(Control.monocle.data)<-"cluster.rev"
Control.monocle.data.f<-subset(Control.monocle.data,subset=nCount_RNA>3000&nFeature_RNA>1000)
DimPlot(Control.monocle.data.f)
dim(Control.monocle.data.f)
[1] 26197 19351

#Downsampling
set.seed(1234)
downsample.id<-sample(colnames(Control.monocle.data.f),length(colnames(Control.monocle.data.f))/2,replace = F)
Control.monocle.test<-subset(Control.monocle.data.f,cells = downsample.id)
dim(Control.monocle.test)
[1] 26197  9675
table(Control.monocle.test@meta.data$cluster.rev)
Artery    Capillary Large Artery         Vein       Venule 
1387         5939          800          331         1218 
DimPlot(Control.monocle.test)
save(Control.monocle.test,file = "Control.nonocle.test.RData")

#Run standard monocle2 pipeline
control_ann<-Control.monocle.test@meta.data
control_ann$celltype<-Idents(Control.monocle.test)

gene_ann<-data.frame(gene_short_name=rownames(Control.monocle.test@assays$RNA),
                     row.names = rownames(Control.monocle.test@assays$RNA))

head(gene_ann)

pd<-new("AnnotatedDataFrame",data=control_ann)
fd<-new("AnnotatedDataFrame",data=gene_ann)
ce=as.data.frame(Control.monocle.test@assays$RNA@counts)

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

#add Seurat celltype annotation to monocle object
pData(cds)$Cluster=pData(cds)$celltype

cds<-setOrderingFilter(cds,g.select$X1)
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,method="DDRTree")
cds<-orderCells(cds)

#Fig2C visulization
pdf("pseudotime__test_result.pdf",width = 11,height = 6)
p1<-plot_cell_trajectory(cds,color_by="Cluster",cell_size = 0.8,cols=cols)+
  scale_color_manual(breaks = c("Large Artery", "Artery", "Capillary","Venule","Vein"), values=c("#E77D70","#C59733","#86AC34","#54BA6F","#55BBC2")) 
p2<-plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.8)
p1|p2
dev.off()

#Identify zonation genens
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res_1<-diff_test_res[order(diff_test_res$qval),]
write.xlsx(diff_test_res_1,file = "Control_EC_pseudotime_DEGs.xlsx",rowNames=T,colNames=T)

#choose top 1000 genes as zonation genes
sig_gene_names <- as.character(diff_test_res_1$gene_short_name[1:1000])

newdata<-data.frame(Pseudotime = seq(min(pData(cds)$Pseudotime), 
                                     max(pData(cds)$Pseudotime),length.out = 200))
genSmoothCurves_mat<-genSmoothCurves(cds[sig_gene_names,],
                                     new_data = newdata,
                                     cores = 10)

#heatmap visulization
library(pheatmap)
pdf("Control_EC_zonation_with_cluster.pdf",width = 6,height = 8)
p<-pheatmap(log10(genSmoothCurves_mat+1),
         scale = "row",
         cluster_rows = T,clustering_distance_rows = "euclidean",clustering_method = "ward.D2",
         cluster_cols = F,
         show_rownames = F,
         show_colnames = F,cutree_rows = 7)
print(p)
dev.off()

#extract heatmap cluster information
geneClust <- p$tree_row
geneCluster2Three<- cutree(geneClust,k=7)

newOrder=genSmoothCurves_mat[p$tree_row$order,]
colnames(newOrder)[ncol(newOrder)]="Cluster"

geneCluster2Three.table<-data.frame(geneCluster2Three)
geneCluster2Three.table$Gene<-rownames(geneCluster2Three.table)
table(geneCluster2Three.table$geneCluster2Three)
1   2   3   4   5   6   7 
156 189 164 145  83  54 209 
geneCluster2Three.table<-geneCluster2Three.table[rownames(newOrder),]
write.xlsx(geneCluster2Three.table,file = "Control_genes_heatmap_order.xlsx",rowNames=F,colNames=T)

#reorder heatmap genes
cluster1.genes<-data.frame(Gene=rownames(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three==6)),Cluster=rep("1",54))
cluster2.genes<-data.frame(Gene=rownames(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three==5)),Cluster=rep("2",83))
cluster3.genes<-data.frame(Gene=rownames(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three==1|geneCluster2Three.table$geneCluster2Three==2)),Cluster=rep("3",345))
cluster4.genes<-data.frame(Gene=rownames(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three==4)),Cluster=rep("4",145))
cluster5.genes<-data.frame(Gene=rownames(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three==7)),Cluster=rep("5",209))
cluster6.genes<-data.frame(Gene=rownames(subset(geneCluster2Three.table,geneCluster2Three.table$geneCluster2Three==3)),Cluster=rep("6",164))



zonation.cluster.genes<-rbind(cluster1.genes,cluster2.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster3.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster4.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster5.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster6.genes)
zonation.cluster.genes<-rbind(zonation.cluster.genes,cluster7.genes)

#Fig2D zonation heatmap
newOrder<-data.frame(newOrder)
newOrder<-newOrder[intersect(zonation.cluster.genes$Gene,rownames(newOrder)),]
pdf("Control_EC_zonation.pdf",width = 6,height = 8)
p1<-pheatmap(log10(newOrder+1),
            scale = "row",
            cluster_rows = F,clustering_distance_rows = "euclidean",clustering_method = "ward.D2",
            cluster_cols = F,
            show_rownames = F,
            show_colnames = F,gaps_row = c(54,137,482,627,836))
print(p1)
dev.off()

#extract cell order
cells.pseudo.table<-data.frame(cell=colnames(cds),
                               pesudotime=cds@phenoData@data$Pseudotime,
                               celltype=cds@phenoData@data$Cluster)
cells.pseudo.table<-cells.pseudo.table[order(cells.pseudo.table$pesudotime),]
#set colors for each celtype
color<-cells.pseudo.table$celltype
values=c("#E77D70","#C59733","#86AC34","#54BA6F","#55BBC2")
color<-gsub("Large Artery","#E77D70",color)
color<-gsub("Artery","#C59733",color)
color<-gsub("Capillary","#86AC34",color)
color<-gsub("Venule","#54BA6F",color)
color<-gsub("Vein","#55BBC2",color)
cells.pseudo.table$color<-color

#barplot of each genes in pseudotime
dat.seurat.normalized=exp(GetAssayData(object = Control.monocle.test))-1
dat.seurat.normalized.f<-dat.seurat.normalized[as.character(zonation.cluster.genes$Gene),]
dat.seurat.normalized.f=round(dat.seurat.normalized.f, 0)
dat.seurat.normalized.f<-data.frame(dat.seurat.normalized.f[,intersect(cells.pseudo.table$cell,colnames(dat.seurat.normalized.f))])
dim(dat.seurat.normalized.f)
[1] 1000 9675

#Fig2D zonation genes barplot
g.select<-c("CD74","CD99","ELN","GJA5","CAV1","SEMA3G","CAVIN2","CA4",
            "MYO1B","CA2","ACKR1","CXCL11")

pdf("Control_zonation_genes_barplot_with_color.PDF",8,16)
par(mfrow=c(6,2))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))
for(i in 1:length(g.select)){
  t.g=g.select[i];
  p<-barplot(as.numeric((dat.seurat.normalized.f[t.g,])), col=color,main=t.g,border = NA,ylab="counts",space = 0.1,byt="o")
  box()
  curve<-loess(as.numeric((dat.seurat.normalized.f[t.g,]))~c(1:ncol(dat.seurat.normalized.f)))
  lines(p,predict(curve),col="red")
}
dev.off()


#curve plot
curve.dat=GetAssayData(object = Control.monocle.test)
curve.dat.f<-curve.dat[as.character(zonation.cluster.genes$Gene),]
curve.dat.f<-data.frame(curve.dat.f[,intersect(cells.pseudo.table$cell,colnames(curve.dat.f))])
dim(curve.dat.f)
[1] 1000 9675


#Fig2D curve plot for each cluster
pdf("Control_zonation_cluster_curveplot.PDF",8,14)
par(mfrow=c(5,2))
par(mar=c(5,4,2,3))
par(oma = c(1, 1, 1, 1.1))

cluster1.curve.data<-curve.dat.f[cluster1.genes,]
cluster1.mean<-apply(cluster1.curve.data, 2, mean)
cluster1.first.curve<-loess(as.numeric(cluster1.curve.data[1,])~c(1:ncol(cluster1.curve.data)))

plot(predict(cluster1.first.curve),col="gray",type="l",ylim=c(0,5),xaxt="n",main="A-V",ylab="Sequence Counts",xlab="")
for(i in 2:length(cluster1.genes)){
    t.g=cluster1.genes[i];
    curve<-loess(as.numeric((cluster1.curve.data[t.g,]))~c(1:ncol(cluster1.curve.data)))
    box()
    lines(predict(curve),col="gray",type="l")
}
cluster1_mean_curve<-loess(cluster1.mean~c(1:ncol(cluster1.curve.data)))
lines(predict(cluster1_mean_curve),col="red",type="l",lwd=2)


cluster2.curve.data<-curve.dat.f[cluster2.genes,]
cluster2.mean<-apply(cluster2.curve.data, 2, mean)
max(cluster2.curve.data)

cluster2.first.curve<-loess(as.numeric(cluster2.curve.data[1,])~c(1:ncol(cluster2.curve.data)))
plot(predict(cluster2.first.curve),col="gray",type="l",ylim=c(0,5),xaxt="n",main="A",ylab="Sequence Counts",xlab="")
for(i in 2:length(cluster2.genes)){
  t.g=cluster2.genes[i];
  curve<-loess(as.numeric((cluster2.curve.data[t.g,]))~c(1:ncol(cluster2.curve.data)))
  box()
  lines(predict(curve),col="gray",type="l")
}
cluster2_mean_curve<-loess(cluster2.mean~c(1:ncol(cluster2.curve.data)))
lines(predict(cluster2_mean_curve),col="red",type="l",lwd=2)


cluster3.curve.data<-curve.dat.f[cluster3.genes,]
cluster3.mean<-apply(cluster3.curve.data, 2, mean)
max(cluster3.curve.data)

cluster3.first.curve<-loess(as.numeric(cluster3.curve.data[1,])~c(1:ncol(cluster3.curve.data)))
plot(predict(cluster3.first.curve),col="gray",type="l",ylim=c(0,5),xaxt="n",main="A-C",ylab="Sequence Counts",xlab="")
for(i in 2:length(cluster3.genes)){
  t.g=cluster3.genes[i];
  curve<-loess(as.numeric((cluster3.curve.data[t.g,]))~c(1:ncol(cluster3.curve.data)))
  box()
  lines(predict(curve),col="gray",type="l")
}
cluster3_mean_curve<-loess(cluster3.mean~c(1:ncol(cluster3.curve.data)))
lines(predict(cluster3_mean_curve),col="red",type="l",lwd=2)


cluster4.curve.data<-curve.dat.f[cluster4.genes,]
cluster4.mean<-apply(cluster4.curve.data, 2, mean)
max(cluster4.curve.data)

cluster4.first.curve<-loess(as.numeric(cluster4.curve.data[1,])~c(1:ncol(cluster4.curve.data)))
plot(predict(cluster4.first.curve),col="gray",type="l",ylim=c(0,5),xaxt="n",main="C",ylab="Sequence Counts",xlab="")
for(i in 2:length(cluster4.genes)){
  t.g=cluster4.genes[i];
  curve<-loess(as.numeric((cluster4.curve.data[t.g,]))~c(1:ncol(cluster4.curve.data)))
  box()
  lines(predict(curve),col="gray",type="l")
}
cluster4_mean_curve<-loess(cluster4.mean~c(1:ncol(cluster4.curve.data)))
lines(predict(cluster4_mean_curve),col="red",type="l",lwd=2)

cluster5.curve.data<-curve.dat.f[cluster5.genes,]
cluster5.mean<-apply(cluster5.curve.data, 2, mean)
max(cluster5.curve.data)

cluster5.first.curve<-loess(as.numeric(cluster5.curve.data[1,])~c(1:ncol(cluster5.curve.data)))
plot(predict(cluster5.first.curve),col="gray",type="l",ylim=c(0,5),xaxt="n",main="V-C",ylab="Sequence Counts",xlab="")
for(i in 2:length(cluster5.genes)){
  t.g=cluster5.genes[i];
  curve<-loess(as.numeric((cluster5.curve.data[t.g,]))~c(1:ncol(cluster5.curve.data)))
  box()
  lines(predict(curve),col="gray",type="l")
}
cluster5_mean_curve<-loess(cluster5.mean~c(1:ncol(cluster5.curve.data)))
lines(predict(cluster5_mean_curve),col="red",type="l",lwd=2)


cluster6.curve.data<-curve.dat.f[cluster6.genes,]
cluster6.mean<-apply(cluster6.curve.data, 2, mean)
max(cluster6.curve.data)

cluster6.first.curve<-loess(as.numeric(cluster6.curve.data[1,])~c(1:ncol(cluster6.curve.data)))
plot(predict(cluster6.first.curve),col="gray",type="l",ylim=c(0,5),xaxt="n",main="V",ylab="Sequence Counts",xlab="")
for(i in 2:length(cluster6.genes)){
  t.g=cluster6.genes[i];
  curve<-loess(as.numeric((cluster6.curve.data[t.g,]))~c(1:ncol(cluster6.curve.data)))
  box()
  lines(predict(curve),col="gray",type="l")
}
cluster6_mean_curve<-loess(cluster6.mean~c(1:ncol(cluster6.curve.data)))
lines(predict(cluster6_mean_curve),col="red",type="l",lwd=2)

dev.off()

#############################################################################################################
######################################### LGG EC Zonation ###################################################
#############################################################################################################
Idents(pure.EC.Seurat)<-"tumour.type"
lower_grade_EC_seurat<-subset(pure.EC.Seurat,ident="Lower-grade glioma")
dim(lower_grade_EC_seurat)
[1] 26197 21222
save(lower_grade_EC_seurat,file = "Lower.grade.monocle.data.RData")

load("Lower.grade.monocle.data.RData")
Idents(lower_grade_EC_seurat)<-"cluster.rev"
Lower.grade.monocle.data.f<-subset(lower_grade_EC_seurat,subset=nCount_RNA>3000&nFeature_RNA>1000)
DimPlot(Lower.grade.monocle.data.f)
dim(Lower.grade.monocle.data.f)
[1] 26197 11068

set.seed(1234)

Idents(Lower.grade.monocle.data.f)<-"cluster.rev"
Large.artery.id<-colnames(subset(Lower.grade.monocle.data.f,ident="Large Artery"))
Artery.id<-colnames(subset(Lower.grade.monocle.data.f,ident="Artery"))
Capillary.id<-colnames(subset(Lower.grade.monocle.data.f,ident="Capillary"))
Venule.id<-colnames(subset(Lower.grade.monocle.data.f,ident="Venule"))
Vein.id<-colnames(subset(Lower.grade.monocle.data.f,ident="Vein"))
TEC1.id<-colnames(subset(Lower.grade.monocle.data.f,ident="TEC1"))
TEC2.id<-colnames(subset(Lower.grade.monocle.data.f,ident="TEC2"))

Lower.grade.monocle.test<-subset(Lower.grade.monocle.data.f,cells = c(Large.artery.id,Artery.id,
                                                                      Capillary.id,Venule.id,
                                                                      Vein.id,TEC1.id,TEC2.id))
table(Lower.grade.monocle.test@meta.data$cluster.rev)
Artery    Capillary Large Artery         TEC1         TEC2         Vein       Venule 
1767         4730         1201           30          200          604         2454 
DimPlot(Lower.grade.monocle.test)

save(Lower.grade.monocle.test,file = "Lower.grade.monocle.test.with.control.zonation.genes.RData")

Lower.grade_ann<-Lower.grade.monocle.test@meta.data
Lower.grade_ann$celltype<-Idents(Lower.grade.monocle.test)

gene_ann<-data.frame(gene_short_name=rownames(Lower.grade.monocle.test@assays$RNA),
                     row.names = rownames(Lower.grade.monocle.test@assays$RNA))

head(gene_ann)

pd<-new("AnnotatedDataFrame",data=Lower.grade_ann)
fd<-new("AnnotatedDataFrame",data=gene_ann)
ce=as.data.frame(Lower.grade.monocle.test@assays$RNA@counts)

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

cds<-setOrderingFilter(cds,as.character(zonation.cluster.genes$Gene))
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,method="DDRTree")
cds<-orderCells(cds)
plot_cell_trajectory(cds, color_by = "State",cell_size = 1)
cds<-orderCells(cds,root_state = 1)

pdf("Lower_grade_pseudotime_result_with_control_zonation_genes.pdf",width = 10,height = 6)
p1<-plot_cell_trajectory(cds,color_by="Cluster",cell_size = 0.7,show_branch_points = F)+
  scale_color_manual(breaks = c("Large Artery", "Artery", "Capillary","Venule","Vein","TEC1","TEC2"),values=c("#E77D70","#C59733","#86AC34","#54BA6F","#55BBC2","#539BD3","#927CB5"))
p2<-plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.7,show_branch_points = F)
p1|p2
dev.off()

save(cds,file = "LGG_EC_monocle_result_with_control_zonation_genes.RData")




#############################################################################################################
######################################### GBM EC Zonation ###################################################
#############################################################################################################
monocle.seurat<-subset(pure.EC.Seurat,idents=c("Large Artery","Artery","Capillary","Venule","Vein","TEC1","TEC2"))
Idents(monocle.seurat)<-"tumour.type"
monocle.seurat<-subset(monocle.seurat,idents=c("Higher-grade glioma"))
dim(monocle.seurat)
Idents(monocle.seurat)<-"cluster.rev"

save(monocle.seurat,file = "monocle.seurat.RData")

higher_ann<-monocle.seurat@meta.data
higher_ann$celltype<-Idents(monocle.seurat)

gene_ann<-data.frame(gene_short_name=rownames(monocle.seurat@assays$RNA),
                     row.names = rownames(monocle.seurat@assays$RNA))

head(gene_ann)

pd<-new("AnnotatedDataFrame",data=higher_ann)
fd<-new("AnnotatedDataFrame",data=gene_ann)
ce=as.data.frame(monocle.seurat@assays$RNA@counts)

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

cds<-setOrderingFilter(cds,as.character(zonation.cluster.genes$Gene))
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2,method="DDRTree")
cds<-orderCells(cds)

pdf("pseudotime_result_with_control_zonation_genes.pdf",width = 10,height = 6)
p1<-plot_cell_trajectory(cds,color_by="Cluster",cell_size = 0.8,cols=cols,show_branch_points = F,show_backbone = F)+
  scale_color_manual(breaks = c("Large Artery", "Artery", "Capillary","Venule","Vein","TEC1","TEC2"), values=c("#E77D70","#C59733","#86AC34","#54BA6F","#55BBC2","#539BD3","#927CB5")) 
p2<-plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 0.8,show_branch_points = F,show_backbone = F)
p1|p2
dev.off()

save(cds,file = "GBM_EC_monocle_result_with_control_zonation_genes.RData")



#############################################################################################################
######################################### Session Info ######################################################
#############################################################################################################
> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X  11.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] splines   stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
[10] base     

other attached packages:
[1] pheatmap_1.0.12     patchwork_1.1.1     openxlsx_4.2.4      monocle_2.14.0     
[5]   irlba_2.3.3         VGAM_1.1-5          Biobase_2.46.0     BiocGenerics_0.32.0
[9]  Matrix_1.3-4        clustree_0.4.4      ggraph_2.1.0.9000    Seurat_3.1.1
[13] ggplot2_3.4.0        
