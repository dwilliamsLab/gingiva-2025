library(Seurat)
options(Seurat.object.assay.version = "v5")
library(tidyverse)

Epithelial <- subset(Integrated, idents=c("Epithelial"))

#too few to re-integrate

Epithelial <- ScaleData(Epithelial)
Epithelial <- FindVariableFeatures(Epithelial, nfeatures=4000)
Epithelial <- RunPCA(Epithelial)

ElbowPlot(Epithelial, ndims=50, reduction="pca")

Epithelial <- FindNeighbors(Epithelial, reduction = "integrated.rpca", dims = 1:30)
Epithelial <- FindClusters(Epithelial, resolution = 1)
Epithelial <- RunUMAP(Epithelial, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rna")

DimPlot(Epithelial, reduction = "umap.rna", pt.size = 0.2, label=TRUE) &NoLegend()


DefaultAssay(Epithelial) = "dsbADT"

Epithelial <- ScaleData(Epithelial, features = rownames(Epithelial@assays$dsbADT@features))
VariableFeatures(Epithelial) <- rownames(Epithelial@assays$dsbADT@features)
Epithelial <- RunPCA(Epithelial, reduction.name = 'apca')

Epithelial <- FindNeighbors(Epithelial, reduction = "integrated.adt.rpca", dims = 1:30)
Epithelial <- FindClusters(Epithelial, resolution = 0.5)
Epithelial <- RunUMAP(Epithelial, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="adt.umap")

DimPlot(Epithelial, reduction = "adt.umap", pt.size = 0.2, label=TRUE) &NoLegend()

# WNN integration

Epithelial <- FindMultiModalNeighbors(
  Epithelial, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")

Epithelial <- RunUMAP(Epithelial, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Epithelial <- FindClusters(Epithelial, graph.name = "wsnn", algorithm = 3, resolution = seq(from = 0.1, 
                                                                                     to = 1.0, 
                                                                                     by = 0.1))

Idents(Epithelial) <- "wsnn_res.0.5"
DimPlot(Epithelial, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.5) +theme_void() &NoLegend()


DefaultAssay(Epithelial) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Epithelial, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

DefaultAssay(Epithelial) <-"dsbADT"
Protein.markers.1 <- FindAllMarkers(Epithelial, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)


FeaturePlot(Epithelial, features=c("DCT","PMEL"), reduction="wnn.umap", min.cutoff=0) #no melanocytes

#Resolution 1: Cluster 0 express IGKC and few other genes - remove
Epithelial <- subset(Epithelial, idents=c("0"),invert=TRUE)
#Run from FindMultiModalNeighbors


Idents(Epithelial) <- "wsnn_res.0.5"
levels(Epithelial) <- c("4","1","0","3","2")
Epithelial <- RenameIdents(Epithelial, "1"="Epithelial 1","0"="Epithelial 2","3"="Epithelial 3","2"="Epithelial 4","4"="Cycling Epithelial")

Epithelial$EpithelialClusters <- Epithelial@active.ident
Idents(Epithelial) <- "EpithelialClusters"
Epithelial$FineClusters <- Epithelial@active.ident
Idents(Epithelial) <- "FineClusters"

rownames(Epithelial@assays$dsbADT@layers$counts) = adtNameData$FinalNames
rownames(Epithelial@assays$dsbADT@layers$data) = adtNameData$FinalNames
rownames(Epithelial@assays$dsbADT@features) = adtNameData$FinalNames
rownames(Epithelial@assays$ADT@counts) = adtNameData$FinalNames
rownames(Epithelial@assays$ADT@data) = adtNameData$FinalNames


DimPlot(Epithelial, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.5) +theme_void() &NoLegend()

Idents(Epithelial) <- "EpithelialClusters"
revEpithelial <- rev(levels(Epithelial))
levels(Epithelial) <-revEpithelial
Epithelial$revEpithelialClusters <-Epithelial@active.ident
Idents(Epithelial) <-"revEpithelialClusters"

revRNAFeatures <-rev(RNAFeatures)
revADTFeatures <-rev(ADTFeatures)



DefaultAssay(Epithelial) <-"RNA"
RNAFeatures <- c(
  "KIF2C","DEPDC1","MKI67", #4 Prolif
  "CXCL8", "CEACAM1","LCN2","SPRR3", "ODAPH","ODAM", #2
 "PDE7B", "MEG3","SULF1",#3
  "MMP13","MMP7","FIBIN","NNMT",  #0
  "KRT6C","DSG1","CA12", "FGFBP1"  #1
  
)

RNAdotplot <-DotPlot(Epithelial, features = RNAFeatures, col.min=0.1, scale.by="radius", dot.min=0.1, 
                     scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") 

RNAdotplot

RNAdotplot <-DotPlot(Epithelial, features = revRNAFeatures, col.min=0.1, scale.by="radius", dot.min=0.1, 
                     scale=TRUE) + 
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
RNAdotplot

DefaultAssay(Epithelial) <-"ADT"
ADTFeatures <- c(
  "CD24_p", "CD73_p","CD71_p","CD54_p","CD142_p", #2
  "CD40_p", #3
"PDPN_p","CD49f_p",#0
"NRP-1_p","IFNGR1_p","HLA-A_p" #1
)

ADTdotplot <- DotPlot(Epithelial, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot 

ADTdotplot <- DotPlot(Epithelial, features = revADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low = "white", mid = "white",high="black", midpoint = 0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

ADTdotplot

saveRDS(Epithelial, "Epithelial 24.12.14.rds")
