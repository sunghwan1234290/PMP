library(Seurat)
library(metap)
library(multtest)
library(harmony)
library(cluster)

#Data set run anlaysis workflow
#QC and selecting cells for further analysis
#read data set 
GSM3755693.data <- Read10X(data.dir = "./GSM3755693_feature_bc_matrix")
#create seurat object
GSM3755693 <- CreateSeuratObject(counts = GSM3755693.data, project = "Normal-Peritoneum",min.cells = 3, min.features = 200)
#add sample name
GSM3755693$stim <- "healthy-1"
#QC and selecting cells for further analysis
GSM3755693[["percent.mt"]] <- PercentageFeatureSet(GSM3755693, pattern = "^MT-")
GSM3755693 <- subset(GSM3755693, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#Normalizing the data
GSM3755693 <- NormalizeData(GSM3755693, verbose = FALSE)
#Identification of highly variable features (feature selection)
GSM3755693 <- FindVariableFeatures(GSM3755693, selection.method = "vst", nfeatures = 2000)

GSM3755694.data <- Read10X(data.dir = "./GSM3755694_feature_bc_matrix")
GSM3755694 <- CreateSeuratObject(counts = GSM3755694.data, project = "Normal-Peritoneum",min.cells = 3, min.features = 200)
GSM3755694$stim <- "healthy-2"
GSM3755694[["percent.mt"]] <- PercentageFeatureSet(GSM3755694, pattern = "^MT-")
GSM3755694 <- subset(GSM3755694, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
GSM3755694 <- NormalizeData(GSM3755694, verbose = FALSE)
GSM3755694 <- FindVariableFeatures(GSM3755694, selection.method = "vst", nfeatures = 2000)

GSM3755695.data <- Read10X(data.dir = "./GSM3755695_feature_bc_matrix")
GSM3755695 <- CreateSeuratObject(counts = GSM3755695.data, project = "Normal-Peritoneum",min.cells = 3, min.features = 200)
GSM3755695$stim <- "healthy-3"
GSM3755695[["percent.mt"]] <- PercentageFeatureSet(GSM3755695, pattern = "^MT-")
GSM3755695 <- subset(GSM3755695, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
GSM3755695 <- NormalizeData(GSM3755695, verbose = FALSE)
GSM3755695 <- FindVariableFeatures(GSM3755695, selection.method = "vst", nfeatures = 2000)

PMP1.data <- Read10X(data.dir = "./GSM7119899_filtered_feature_bc_matrix/")
PMP1 <- CreateSeuratObject(counts = PMP1.data, project = "PMP",min.cells = 3, min.features = 200)
PMP1$stim <- "PMP-1"
PMP1[["percent.mt"]] <- PercentageFeatureSet(PMP1, pattern = "^MT-")
PMP1 <- subset(PMP1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
PMP1 <- NormalizeData(PMP1, verbose = FALSE)
PMP1 <- FindVariableFeatures(PMP1, selection.method = "vst", nfeatures = 2000)

PMP2.data <- Read10X(data.dir = "./GSM7119900_filtered_feature_bc_matrix/")
PMP2 <- CreateSeuratObject(counts = PMP2.data, project = "PMP",min.cells = 3, min.features = 200)
PMP2$stim <- "PMP-2"
PMP2[["percent.mt"]] <- PercentageFeatureSet(PMP2, pattern = "^MT-")
PMP2 <- subset(PMP2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
PMP2 <- NormalizeData(PMP2, verbose = FALSE)
PMP2 <- FindVariableFeatures(PMP2, selection.method = "vst", nfeatures = 2000)

PMP3.data <- Read10X(data.dir = "./GSM7119901_filtered_feature_bc_matrix/")
PMP3 <- CreateSeuratObject(counts = PMP3.data, project = "PMP",min.cells = 3, min.features = 200)
PMP3$stim <- "PMP-3"
PMP3[["percent.mt"]] <- PercentageFeatureSet(PMP3, pattern = "^MT-")
PMP3 <- subset(PMP3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
PMP3 <- NormalizeData(PMP3, verbose = FALSE)
PMP3 <- FindVariableFeatures(PMP3, selection.method = "vst", nfeatures = 2000)

PMP4.data <- Read10X(data.dir = "./GSM7119902_filtered_feature_bc_matrix/")
PMP4 <- CreateSeuratObject(counts = PMP4.data, project = "PMP",min.cells = 3, min.features = 200)
PMP4$stim <- "PMP-4"
PMP4[["percent.mt"]] <- PercentageFeatureSet(PMP4, pattern = "^MT-")
PMP4 <- subset(PMP4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
PMP4 <- NormalizeData(PMP4, verbose = FALSE)
PMP4 <- FindVariableFeatures(PMP4, selection.method = "vst", nfeatures = 2000)

PMP5.data <- Read10X(data.dir = "./GSM7119903	_filtered_feature_bc_matrix/")
PMP5 <- CreateSeuratObject(counts = PMP5.data, project = "PMP",min.cells = 3, min.features = 200)
PMP5$stim <- "PMP-5"
PMP5[["percent.mt"]] <- PercentageFeatureSet(PMP5, pattern = "^MT-")
PMP5 <- subset(PMP5, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
PMP5 <- NormalizeData(PMP5, verbose = FALSE)
PMP5 <- FindVariableFeatures(PMP5, selection.method = "vst", nfeatures = 2000)

#Perform integration
PMP.anchors <- FindIntegrationAnchors(object.list = list(GSM3755693,GSM3755694,GSM3755695,PMP1,PMP2,PMP3,PMP4,PMP5), dims = 1:20)
PMP.combined <- IntegrateData(anchorset = PMP.anchors, dims = 1:20)

#Perform an integrated analysis
DefaultAssay(PMP.combined) <- "integrated"

#Scaling the data
PMP.combined <- ScaleData(PMP.combined, verbose = FALSE)
#Perform linear dimensional reduction
PMP.combined <- RunPCA(PMP.combined, npcs = 30, verbose = FALSE)
#Harmony Integration
PMP.combined <- RunHarmony(PMP.combined, group.by.vars = "orig.ident")
#Run non-linear dimensional reduction
PMP.combined <- RunUMAP(PMP.combined,reduction = "harmony", dims = 1:30) 
ElbowPlot(PMP.combined)
#Cluster the cells (Check the elbow plot to determine the optimal number of PCs)
PMP.combined <- FindNeighbors(PMP.combined, reduction = "harmony", dims = 1:6)
#Due to differences in R and Seurat versions, the number of clusters after setting the resolution may vary. Please ensure it is adjusted to 10 clusters.
PMP.combined <- FindClusters(PMP.combined, resolution = 0.25) 
#UMAP visualization
DimPlot(PMP.combined, reduction = "umap", label = TRUE,label.size = 10, repel = TRUE)
DimPlot(immune.combined, reduction = "umap", group.by = "orig.ident") + scale_colour_manual(values = alpha(c("blue", "red"), 0.1))
DimPlot(PMP.combined, reduction = "umap", split.by = "stim",ncol = 3)

#Identify conserved cell type markers
DefaultAssay(PMP.combined) <- "RNA"
#FeaturePlot,VlnPlot significant gene visualization
FeaturePlot(PMP.combined, features = c("CDH1", "MUC1", "EPCAM", "ZEB2", "SPARC", "COL3A1","CDH2"),order = TRUE)
VlnPlot(PMP.combined, features = c("CDH1", "MUC1","EPCAM","ZEB2","SPARC","COL3A1","CDH2"))
#Top marker genes by cluster 
PMP.markers <- FindAllMarkers(PMP.combined, only.pos = TRUE, min.pct = 0.2,logfc.threshold = 0.4)
write.table(PMP.markers, ./Markers.txt", quote=FALSE, sep='\t', col.names = TRUE)

#Dotplot of genes representing each cluster
DotPlot(object = PMP.combined, features = c("IL7R",
                                               "BATF",
                                               "EPCAM",
                                               "MUC1",
                                               "CD8A", 
                                               "CD8B",
                                               "COL3A1",
                                               "SPARC",
                                               "S100A8",
                                               "S100A9",
                                               "CD19",
                                               "CD37",
                                               "CD14",
                                               "CD74",
                                               "CD38",
                                               "XBP1",
                                               "IL6",
                                               "VCAM1",
                                               "CD300E",
                                               "MS4A6A"), cols = c("lightgrey", "blue"),  dot.scale = 14, scale.by = "size") 

#Violin plot of epithelial marker gene expression in healthy and PMP samples of each cluster
plots <- VlnPlot(PMP.combined, features = c("CDH1", "EPCAM", "MUC1"), split.by = "orig.ident",pt.size = 1, combine = FALSE,split.plot = TRUE) 
CombinePlots(plots = plots, ncol = 1)
#Violin plot of mesenchymal marker gene expression in healthy and PMP samples of each cluster
plots<-VlnPlot(PMP.combined, features = c("ZEB2","SPARC","COL3A1","CDH2"), split.by = "orig.ident",pt.size = 1, combine = FALSE,split.plot = TRUE)
CombinePlots(plots = plots, ncol = 1)
                                               
#Pseudo bulk using AggregateExpression
pseudo_pmp <- AggregateExpression(PMP.combined, group.by = "stim", assays = "RNA", slot = "data")
write.table(pseudo_pmp,file="./pseudo_pmp.txt",quote=FALSE,sep="\t")

#subset analysis : epithelial cell cluster
#cluster designation
cluster_of_interest <- 1
epithelial_obj <- subset(PMP.combined, idents = cluster_of_interest)
#standard workflow 
epithelial_obj <- RunPCA(epithelial_obj, npcs = 30, verbose = FALSE)  
epithelial_obj <- RunUMAP(epithelial_obj, dims = 1:10)  
ElbowPlot(epithelial_obj)
#Cluster the cells (Check the elbow plot to determine the optimal number of PCs)
epithelial_obj <- FindNeighbors(epithelial_obj, dims = 1:5)
#Due to differences in R and Seurat versions, the number of clusters after setting the resolution may vary.
epithelial_obj <- FindClusters(epithelial_obj, resolution = 0.07)
#UMAP visualization
DimPlot(epithelial_obj, reduction = "umap", label = TRUE,label.size = 10, repel = TRUE)

#epithelial cell cluster_Normal-Peritoneum analysis
healthy_cells <- subset(epithelial_obj, subset = orig.ident == "Normal-Peritoneum")
healthy_cells <- RunUMAP(healthy_cells, dims = 1:10)
ElbowPlot(healthy_cells)
healthy_cells <- FindNeighbors(healthy_cells, dims = 1:5)
healthy_cells <- FindClusters(healthy_cells, resolution = 0.15) 
DefaultAssay(healthy_cells) <- "integrated"
DimPlot(healthy_cells, reduction = "umap", label = TRUE,label.size = 10, repel = TRUE)
DefaultAssay(healthy_cells) <- "RNA"
#find Mesothelial associated genes 
FeaturePlot(healthy_cells, features = c("KRT19", "KRT8", "KRT18", "KRT7"),ncol = 2,order=TRUE)
VlnPlot(healthy_cells, features = c("KRT19", "KRT8", "KRT18", "KRT7"),ncol = 2)
#Vascular endothelial associated genes
FeaturePlot(healthy_cells, features = c("PECAM1", "CDH5","VWF", "CLDN5"),ncol = 2)
VlnPlot(healthy_cells, features = c("PECAM1", "CDH5","VWF", "CLDN5"),ncol = 2)

#epithelial cell cluster_PMP analysis
PMP_cells <- subset(epithelial_obj, subset = orig.ident == "PMP")
PMP_cells <- RunUMAP(PMP_cells, dims = 1:10)
ElbowPlot(PMP_cells)
PMP_cells <- FindNeighbors(PMP_cells, dims = 1:5)
PMP_cells <- FindClusters(PMP_cells, resolution = 0.05) 
DefaultAssay(PMP_cells) <- "integrated"
DimPlot(PMP_cells, reduction = "umap", label = TRUE,label.size = 10, repel = TRUE)
DefaultAssay(PMP_cells) <- "RNA"
#find Gastrointestinal mucinous carcinomas associated genes
FeaturePlot(PMP_cells, features = c("MUC1","MUC2","MUC5AC","MMP7"),ncol = 2,order=TRUE)
VlnPlot(PMP_cells, features = c("MUC1","MUC2","MUC5AC" ,"MMP7"),ncol = 2)
#find fibroblast associated genes
FeaturePlot(PMP_cells, features = c("COL1A1", "COL1A2","DCN", "FBLN1"),ncol = 2,order=TRUE)
VlnPlot(PMP_cells, features = c("COL1A1", "COL1A2","DCN", "FBLN1"),ncol = 2)
#find immune associated genes
FeaturePlot(PMP_cells, features = c("KIT", "CD84","IL1RL1", "TPSB2"),ncol = 2,order=TRUE)
VlnPlot(PMP_cells, features = c("KIT", "CD84","IL1RL1", "TPSB2"),ncol = 2)

#subset analysis : mesenchymal cell cluster
#cluster designation
cluster_of_interest <- 3
#standard workflow 
mesenchymal_obj <- subset(PMP.combined, idents = cluster_of_interest)
mesenchymal_obj <- RunPCA(mesenchymal_obj, npcs = 30, verbose = FALSE)  
mesenchymal_obj <- RunUMAP(mesenchymal_obj, dims = 1:10)  
ElbowPlot(mesenchymal_obj)
#Cluster the cells (Check the elbow plot to determine the optimal number of PCs)
mesenchymal_obj <- FindNeighbors(mesenchymal_obj, dims = 1:2)
#Due to differences in R and Seurat versions, the number of clusters after setting 
mesenchymal_obj <- FindClusters(mesenchymal_obj, resolution = 0.05) 
#UMAP visualization
DimPlot(mesenchymal_obj, reduction = "umap",label.size = 10, repel = TRUE)

#mesenchymal cell cluster_Normal-Peritoneum analysis
mesenchymal_healthy_cells <- subset(mesenchymal_obj, subset = orig.ident == "Normal-Peritoneum")
mesenchymal_healthy_cells <- RunUMAP(mesenchymal_healthy_cells, dims = 1:10)
ElbowPlot(mesenchymal_healthy_cells)
mesenchymal_healthy_cells <- FindNeighbors(mesenchymal_healthy_cells, dims = 1:2)
mesenchymal_healthy_cells <- FindClusters(mesenchymal_healthy_cells, resolution = 0.05) 
DefaultAssay(mesenchymal_healthy_cells) <- "integrated"
DimPlot(mesenchymal_healthy_cells, reduction = "umap", repel = TRUE)
DefaultAssay(mesenchymal_healthy_cells) <- "RNA"
#find Mesothelial associated genes
VlnPlot(mesenchymal_healthy_cells, features = c("KRT19", "KRT8", "KRT18", "KRT7"),ncol = 2)
FeaturePlot(mesenchymal_healthy_cells, features = c("KRT19", "KRT8", "KRT18", "KRT7"),ncol = 2,order=TRUE)
#find Mesothelioma associated genes
VlnPlot(mesenchymal_healthy_cells, features = c("CALB2","KRT5", "MSLN", "SEPT7"),ncol = 2)
FeaturePlot(mesenchymal_healthy_cells, features = c("CALB2","KRT5", "MSLN", "SEPT7"),ncol = 2,order=TRUE)

#mesenchymal cell cluster_PMP analysis
mesenchymal_PMP_cells <- subset(mesenchymal_obj, subset = orig.ident == "PMP")
mesenchymal_PMP_cells <- RunUMAP(mesenchymal_PMP_cells, dims = 1:10)
mesenchymal_PMP_cells <- FindNeighbors(mesenchymal_PMP_cells, dims = 1:2)
mesenchymal_PMP_cells <- FindClusters(mesenchymal_PMP_cells, resolution = 0.1) 
DefaultAssay(mesenchymal_PMP_cells) <- "integrated"
DimPlot(mesenchymal_PMP_cells, reduction = "umap", label.size = 10, repel = TRUE)
DefaultAssay(mesenchymal_PMP_cells) <- "RNA"
#find mesothelial associated genes
VlnPlot(mesenchymal_PMP_cells, features = c("KRT19", "KRT8", "KRT18", "KRT7"),ncol = 2)
FeaturePlot(mesenchymal_PMP_cells, features = c("KRT19", "KRT8", "KRT18", "KRT7"),ncol = 2)
#find mesothelioma associated genes
VlnPlot(mesenchymal_PMP_cells, features = c("CALB2","KRT5", "MSLN", "SEPT7"),ncol = 2)
FeaturePlot(mesenchymal_PMP_cells, features = c("CALB2","KRT5", "MSLN", "SEPT7"),ncol = 2)
