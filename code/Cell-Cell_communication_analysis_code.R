library(CellChat)
library(Seurat)
library(openxlsx)

#load all cells Seurat object
load("all_pure_cell_seurat.RData")

#Subset Seurat object
cellchat.seurat.all<-subset(dat.seurat.integrated.final.pure,idents=c("Endothelial","Mural cell","Perivascular Fibroblast","Oligodendrocyte","Meningeal Fibroblast"))
Idents(cellchat.seurat.all)<-"tumour.type"
cellchat.seurat.GBM<-subset(cellchat.seurat.all,idents="Higher-grade glioma")
save(cellchat.seurat.GBM,file = "CellChat_GBM_Seurat.Rdata")
cellchat.seurat.LGG<-subset(cellchat.seurat.all,idents="Lower-grade glioma")
save(cellchat.seurat.LGG,file = "CellChat_LGG_Seurat.Rdata")
cellchat.seurat.control<-subset(cellchat.seurat.all,idents="Control")
save(cellchat.seurat.control,file = "CellChat_Con_Seurat.Rdata")



#############################################################################################################
##################################### GBM Cell-Cell Communication ###########################################
#############################################################################################################
cellchat.GBM <- createCellChat(cellchat.seurat.GBM@assays$RNA@data)
meta.GBM <- data.frame(cellType = cellchat.seurat.GBM$cluster.rev, row.names =  Cells(cellchat.seurat.GBM))
cellchat.GBM <- addMeta(cellchat.GBM, meta = meta.GBM, meta.name = "cellType")
cellchat.GBM <- setIdent(cellchat.GBM, ident.use = "cellType") # set "labels" as default cell identity
groupSize.GBM <- as.numeric(table(cellchat.GBM@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human
# CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") 
#cellchat@DB <- CellChatDB.use # set the used database in the object
cellchat.GBM@DB <- CellChatDB

cellchat.GBM <- subsetData(cellchat.GBM) # subset the expression data of signaling genes for saving computation cost
cellchat.GBM <- identifyOverExpressedGenes(cellchat.GBM) 
cellchat.GBM <- identifyOverExpressedInteractions(cellchat.GBM)
cellchat.GBM <- projectData(cellchat.GBM, PPI.human)

cellchat.GBM <- computeCommunProb(cellchat.GBM,population.size = FALSE)
cellchat.GBM <- computeCommunProbPathway(cellchat.GBM)
cellchat.GBM <- aggregateNet(cellchat.GBM)

#extract the infered result of each ligand-receptor pair /pathway
GBM.df.net<-subsetCommunication(cellchat.GBM)
GBM.df.net.pathway<-subsetCommunication(cellchat.GBM,slot.name = "netP")
write.xlsx(GBM.df.net.pathway,file = "GBM_cell_interaction_pathway_result.xlsx",rowNames=F,colNames=T)
write.xlsx(GBM.df.net,file = "GBM_cell_interaction_LR_result.xlsx",rowNames=F,colNames=T)


par(mfrow = c(1,2), xpd=TRUE)
library(patchwork)
pdf("GBM_cell_interraction_in_ECM_receptor.pdf",width = 16,height = 10)

netVisual_circle(cellchat.GBM@net$count, vertex.weight = groupSize.GBM, weight.scale = T, label.edge= T, title.name = "Number of interactions in GBM",
                 color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))
netVisual_circle(cellchat.GBM@net$weight, vertex.weight = groupSize.GBM, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                 color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))

dev.off()

# Compute the network centrality scores
cellchat.GBM <- netAnalysis_computeCentrality(cellchat.GBM, slot.name = "netP")

save(cellchat.GBM,file = "GBM_cellchat_result.RData")



#############################################################################################################
##################################### LGG Cell-Cell Communication ###########################################
#############################################################################################################
cellchat.LGG <- createCellChat(cellchat.seurat.LGG@assays$RNA@data)
meta.LGG<- data.frame(cellType = cellchat.seurat.LGG$cluster.rev, row.names =  Cells(cellchat.seurat.LGG))
cellchat.LGG <- addMeta(cellchat.LGG, meta = meta.LGG, meta.name = "cellType")
cellchat.LGG <- setIdent(cellchat.LGG, ident.use = "cellType") # set "labels" as default cell identity
groupSize.LGG <- as.numeric(table(cellchat.LGG@idents)) # number of cells in each cell group

cellchat.LGG@DB <- CellChatDB

cellchat.LGG <- subsetData(cellchat.LGG) # subset the expression data of signaling genes for saving computation cost
cellchat.LGG <- identifyOverExpressedGenes(cellchat.LGG) 
cellchat.LGG <- identifyOverExpressedInteractions(cellchat.LGG)
cellchat.LGG <- projectData(cellchat.LGG, PPI.human)

cellchat.LGG <- computeCommunProb(cellchat.LGG,population.size = FALSE)
cellchat.LGG <- computeCommunProbPathway(cellchat.LGG)
cellchat.LGG <- aggregateNet(cellchat.LGG)

#extract the infered result of each ligand-receptor pair /pathway
LGG.df.net<-subsetCommunication(cellchat.LGG)
LGG.df.net.pathway<-subsetCommunication(cellchat.LGG,slot.name = "netP")
write.xlsx(LGG.df.net.pathway,file = "LGG_cell_interaction_pathway_result.xlsx",rowNames=F,colNames=T)
write.xlsx(LGG.df.net,file = "LGG_cell_interaction_LR_result.xlsx",rowNames=F,colNames=T)


par(mfrow = c(1,2), xpd=TRUE)
pdf("LGG_cell_interraction_in_ECM_receptor.pdf",width = 16,height = 10)

netVisual_circle(cellchat.LGG@net$count, vertex.weight = groupSize.LGG, weight.scale = T, label.edge= T, title.name = "Number of interactions in GBM",
                 color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))
netVisual_circle(cellchat.LGG@net$weight, vertex.weight = groupSize.LGG, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",
                 color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))

dev.off()

mat.LGG <- cellchat.LGG@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat.LGG)) {
  mat2 <- matrix(0, nrow = nrow(mat.LGG), ncol = ncol(mat.LGG), dimnames = dimnames(mat.LGG))
  mat2[i, ] <- mat.LGG[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.LGG, weight.scale = T, edge.weight.max = max(mat.LGG), title.name = rownames(mat.LGG)[i],
                   color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))
}

pathways.show.all.LGG <- cellchat.LGG@netP$pathways

# Compute the network centrality scores
cellchat.LGG <- netAnalysis_computeCentrality(cellchat.LGG, slot.name = "netP")


save(cellchat.LGG,file = "LGG_cellchat_result.RData")



#############################################################################################################
##################################### Control Cell-Cell Communication #######################################
#############################################################################################################
cellchat.Con <- createCellChat(cellchat.seurat.control@assays$RNA@data)
meta.Con <- data.frame(cellType = cellchat.seurat.control$cluster.rev, row.names =  Cells(cellchat.seurat.control))
cellchat.Con <- addMeta(cellchat.Con, meta = meta.Con, meta.name = "cellType")
cellchat.Con <- setIdent(cellchat.Con, ident.use = "cellType") # set "labels" as default cell identity
groupSize.Con <- as.numeric(table(cellchat.Con@idents)) # number of cells in each cell group

cellchat.Con@DB <- CellChatDB

cellchat.Con <- subsetData(cellchat.Con) # subset the expression data of signaling genes for saving computation cost
cellchat.Con <- identifyOverExpressedGenes(cellchat.Con) 
cellchat.Con <- identifyOverExpressedInteractions(cellchat.Con)
cellchat.Con <- projectData(cellchat.Con, PPI.human)

cellchat.Con <- computeCommunProb(cellchat.Con,population.size = FALSE)
cellchat.Con <- computeCommunProbPathway(cellchat.Con)
cellchat.Con <- aggregateNet(cellchat.Con)

#extract the infered result of each ligand-receptor pair /pathway
Con.df.net<-subsetCommunication(cellchat.Con)
Con.df.net.pathway<-subsetCommunication(cellchat.Con,slot.name = "netP")
write.xlsx(Con.df.net.pathway,file = "Control_cell_interaction_pathway_result.xlsx",rowNames=F,colNames=T)
write.xlsx(Con.df.net,file = "Control_cell_interaction_LR_result.xlsx",rowNames=F,colNames=T)

par(mfrow = c(1,2), xpd=TRUE)
pdf("Con_cell_interraction_in_ECM_receptor.pdf",width = 16,height = 10)

netVisual_circle(cellchat.Con@net$count, vertex.weight = groupSize.Con, weight.scale = T, label.edge= T, title.name = "Number of interactions in Control",
                 color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))
netVisual_circle(cellchat.Con@net$weight, vertex.weight = groupSize.Con, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

dev.off()

mat.Con <- cellchat.Con@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat.Con)) {
  mat2 <- matrix(0, nrow = nrow(mat.Con), ncol = ncol(mat.Con), dimnames = dimnames(mat.Con))
  mat2[i, ] <- mat.Con[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.Con, weight.scale = T,label.edge = T, edge.weight.max = max(mat.Con), title.name = rownames(mat.Con)[i],
                   color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))
}


# Compute the network centrality scores
cellchat.Con <- netAnalysis_computeCentrality(cellchat.Con, slot.name = "netP")

save(cellchat.Con,file = "Con_cellchat_result.RData")

#############################################################################################################
##################################### Combine all Communication  Results#####################################
#############################################################################################################
# Define the cell labels to lift up
group.new = levels(cellchat.GBM@idents)
cellchat.Con <- liftCellChat(cellchat.Con, group.new)
cellchat.LGG <- liftCellChat(cellchat.LGG, group.new)
#combine two dataset
object.list <- list(Control = cellchat.Con,LGG=cellchat.LGG, GBM = cellchat.GBM)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#Fig5B circle plot of number of interactions in each condition
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file = "Number_of_interactions_all_CirclePlot.pdf",width = 12,height = 6)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                   color.use=c("#F8766D","#38B600","#D79000","#A3A501","#01BFC4"))
}
dev.off()

#show the total nummber of cell interactions bwtween two conditions
pdf(file = "Number_of_interactions_all.pdf",width = 10,height = 6)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))+scale_fill_manual(values = c("#6787B6","#54A94E","#CE7670",""))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")+scale_fill_manual(values = c("#6787B6","#54A94E","#CE7670"))
print(gg1 + gg2)
dev.off()


#compare the spesific pathway change in certain celltype
pdf(file="signaling_changes_in_each_celltype_ScatterPlot.pdf",width = 16,height = 6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Endothelial")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mural cell")
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Perivascular Fibroblast")
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "High Matrix cells")
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Meningeal Fibroblast")
gg6 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Oligodendrocyte")
print(gg1+gg2)
print(gg3+gg4)
print(gg5+gg6)
dev.off()

#Identify dysfunctional signaling by comparing the communication probabities
pdf(file="changed_signaling_in_Endothelial_to_other.pdf",width = 12,height = 8)
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 3), max.dataset = 3, title.name = "Increased signaling in GBM", angle.x = 45, remove.isolate = T,thresh = 0.05)
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:6),  comparison = c(1, 3), max.dataset = 1, title.name = "Decreased signaling in GBM", angle.x = 45, remove.isolate = T)

gg1 + gg2
dev.off()

save(cellchat,file = "GBM_LGG_Control_combined_cellchat_result.RData")



#############################################################################################################
######################################### Session Info ######################################################
#############################################################################################################
R version 3.6.3 (2020-02-29)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: OS X  11.3

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] CellChat_1.6.0        Biobase_2.46.0        BiocGenerics_0.32.0   ggplot2_3.3.6         igraph_1.2.6          dplyr_1.0.7     Seurat_3.1.1         
  
