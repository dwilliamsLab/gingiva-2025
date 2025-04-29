library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(viridis)
library(readxl)

#Each data has been ran through the CHI CITE pipeline, which uses the following packages:
#Ambient RNA Correction – Decontx
#ADT normalization – DSB
#Identification of doublets – doubletFinder

HV181 <- readRDS("/Users/fraserdr/Desktop/CITEseq/HV181_deep_Final.rds")
HV202 <- readRDS("/Users/fraserdr/Desktop/CITEseq/HV202_deep_Final.rds")
HV214 <- readRDS("/Users/fraserdr/Desktop/CITEseq/HV214_deep_Final.rds")
HV219 <- readRDS("/Users/fraserdr/Desktop/CITEseq/HV219_Final.rds")

dataList <- c(HV181, HV202, HV214, HV219)


# Merge Data
Integrated <- merge(HV181, y = dataList[2:4])

# Integrate RNA data
DefaultAssay(Integrated) <- "RNA"

# Join layers and add cell cycle metadata
Integrated[["RNA"]] <- JoinLayers(Integrated[["RNA"]])

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Integrated <- NormalizeData(Integrated)
Integrated <- CellCycleScoring(Integrated, s.features = s.genes, g2m.features = g2m.genes,
                              set.ident = TRUE)

# run standard analysis workflow

Integrated[["RNA"]] <- split(Integrated[["RNA"]], f = Integrated$orig.ident)

Integrated <- NormalizeData(Integrated)
Integrated <- FindVariableFeatures(Integrated, nfeatures=4000)
Integrated <- ScaleData(Integrated, vars.to.regress=c("S.Score", "G2M.Score"))
Integrated <- RunPCA(Integrated)


### Integrate Datasets (Batch correction)

# Perform Integration
Integrated <- IntegrateLayers(object = Integrated, 
                             method = RPCAIntegration, 
                             orig.reduction = "pca", 
                             new.reduction = "integrated.rpca")

# re-join layers after integration
Integrated[["RNA"]] <- JoinLayers(Integrated[["RNA"]])

# Find Neighbors / Find Clusters / Run UMAP
Integrated <- FindNeighbors(Integrated, reduction = "integrated.rpca", dims = 1:30)
Integrated <- FindClusters(Integrated, resolution = 1)

# Run and Display UMAP
Integrated <- RunUMAP(Integrated, dims = 1:30, reduction = "integrated.rpca", reduction.name="rna.umap")

DimPlot(Integrated, reduction = "rna.umap", pt.size = 0.2) &NoLegend()
DimPlot(Integrated, reduction = "rna.umap", split.by = c("orig.ident"), pt.size = 0.2) &NoLegend()

DefaultAssay(Integrated) = "dsbADT"

# Load new ADT names
adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(Integrated@assays$dsbADT@layers$counts) = adtNameData$FinalNames
rownames(Integrated@assays$dsbADT@layers$data) = adtNameData$FinalNames
rownames(Integrated@assays$dsbADT@features) = adtNameData$FinalNames
rownames(Integrated@assays$ADT@counts) = adtNameData$FinalNames
rownames(Integrated@assays$ADT@data) = adtNameData$FinalNames

Integrated[["dsbADT"]] <- split(Integrated[["dsbADT"]], f = Integrated$orig.ident)

VariableFeatures(Integrated) <- rownames(Integrated@assays$dsbADT@features)
Integrated <- ScaleData(Integrated, features = rownames(Integrated@assays$dsbADT@features))
Integrated <- RunPCA(Integrated, reduction.name = 'apca')


IntegratedIntegrated <- IntegrateLayers(object = Integrated, 
                           method = RPCAIntegration, 
                           assay = "dsbADT",
                           orig.reduction = "apca",
                           new.reduction = "integrated.adt.rpca",
                           features = rownames(Integrated@assays$dsbADT@features))

# re-join layers after integration
Integrated[["dsbADT"]] <- JoinLayers(Integrated[["dsbADT"]])

# Find Neighbors / Find Clusters / Run UMAP
Integrated <- FindNeighbors(Integrated, reduction = "integrated.adt.rpca", dims = 1:30)
Integrated <- FindClusters(Integrated, resolution = 0.5)

# Run and Display UMAP
Integrated <- RunUMAP(Integrated, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="adt.umap")

DimPlot(Integrated, reduction = "adt.umap", group.by = "orig.ident", pt.size = 0.2)IntegrateLayers(object = Integrated, 
                             method = RPCAIntegration, 
                             assay = "dsbADT",
                             orig.reduction = "apca",
                             new.reduction = "integrated.adt.rpca",
                             features = rownames(Integrated@assays$dsbADT@features))

# re-join layers after integration
Integrated[["dsbADT"]] <- JoinLayers(Integrated[["dsbADT"]])

# Find Neighbors / Find Clusters / Run UMAP
Integrated <- FindNeighbors(Integrated, reduction = "integrated.adt.rpca", dims = 1:30)
Integrated <- FindClusters(Integrated, resolution = 0.5)

# Run and Display UMAP
Integrated <- RunUMAP(Integrated, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="adt.umap")
# UMAP
DimPlot(Integrated, reduction = "adt.umap", group.by = "orig.ident", pt.size = 0.2)


Integrated <- FindMultiModalNeighbors(
  Integrated, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")

Integrated <- RunUMAP(Integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Integrated <- FindClusters(Integrated, graph.name = "wsnn", algorithm = 3, resolution = seq(from = 0.1, to = 1, by = 0.1))

DimPlot(Integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2) &NoLegend()


# Find RNA and Protein Markers

Idents(Integrated) <- "wsnn_res.0.1"

DefaultAssay(Integrated) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Integrated, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)



DefaultAssay(Integrated) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(Integrated, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

Integrated <- RenameIdents(Integrated, 
             "0"="T/NK",
             "1"="Plasma","10"="Plasma",
             "2"="Fibroblast",
             "3"="Endothelial",
             "4"="Monocyte/Macrophage/mDC",
             "5"="B",
             "6"="Mast",
             "7"="SMC",
             "8" = "Neutrophil",
             "9"="Epithelial",
             "11"="Cycling",
             "12"="LEC",
             "13"="pDC")



Integrated$CellTypes <- Integrated@active.ident



Idents(Integrated) <- "wsnn_res.1"

Integrated <- RenameIdents(Integrated, 
                           "2"="Fibroblast", "10"="Fibroblast","13"="Fibroblast","16"="Fibroblast","26"="Fibroblast",
                           "0"="Endothelial", "14"="Endothelial",  "17"="Endothelial",
                           "29"="LEC", 
                           "20"="SMC","21"="SMC",
                           "23"="Epithelial", "27"="Epithelial",
                           "3"="CD4 T", "8"="CD4 T", "1"="CD8 T", "15"="CD8 T", "12"="NK",
                           "7"="B",
                           "4"="Plasma","5"="Plasma","6"="Plasma","22"="Plasma", "28"="Plasma",
                           "25"="Cycling",
                           "19"="Monocyte", "18"="Macrophage", "24"="mDC",
                           "9"="Mast",
                           "11"="Neutrophil",
                           "32"="pDC",
                           "30" = "Other", "31"="Other" #30 mixed cell populations, 31 low quality
)

Integrated$HighResCellTypes <- Integrated@active.ident
Idents(Integrated) <- "HighResCellTypes"

DefaultAssay(Integrated) <-"RNA"
RNAfeatures <- c(
  "PTCRA","CLEC4C","LILRA4", #pDC
  "FCGR3B","FFAR2","CSF3R", #Neutrophil
  "MS4A2","TPSAB1","CPA3", #Mast
  "FLT3", "CCSER1", "LGALS2", #mDC
  "MRC1", "MMP9","C1QA", #Macrophage
  "OLR1","CD300E","FCN1", #Monocyte
  "TOP2A","ASPM","MKI67", #Cycling
  "IGHA1","IGHG2","IGLC1", #Plasma
  "FCRL1","BANK1","MS4A1", #B
  "KLRD1","GNLY","NKG7", #NK
  "CCL5","CD8A","GZMK",#CD8 T
  "CTLA4","CCR6", "IL7R", #CD4 T
  "KRT5","FDCSP","S100A2", #Epi
  "NOTCH3","SMOC2","RGS5", #SMC
  "SEMA3D","MMRN1","CCL21", #LEC
  "MECOM","EMCN","VWF", #Endothelial 
  "COL1A1","COL1A2","PDGFRA" #Fibroblast
)

RNAdotplot <- DotPlot(Integrated, features = RNAfeatures, col.min=0.1, dot.min=0.1, scale=TRUE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

DefaultAssay(Integrated) <-"dsbADT"
ADTfeatures <- c(
  "CD62L","CLEC4C","CD123", #pDC
  "CD16","CD15","CD66", #Neutrophil
  "IgE", "FCER1A", "CD33", #Mast
  "CD86","CD1c","HLA-DR", #mDC
  "CD74","CD14","CD64", #Macrophage
  "CR1","CLEC12A","CD11b", #Monocyte
  "CR2","CD44","SLAMF7", #Plasma
  "CD268","CD37","CD20", #B
  "IL2RB","CD94","NCR1", #NK
  "CD57","CD314","CD8", #CD8
  "PD-1","CD69","CD4", #CD4 T"
  "CD55","CD49f","CD71", #Epithelial
  "CD29","CD146", "CD49a",#SMC
  "CD200","CD73","PDPN", #LEC
  "CD112","CD62P","CD31", #Endothelial
  "ICOSLG","CD81","CD13" #Fibroblast
)

ADTdotplot <- DotPlot(Integrated, features = ADTfeatures, col.min=0, dot.min=0.1, scale=TRUE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot 

CD8markers<- FindMarkers(Integrated, ident.1="CD8 T", ident.2="NK", min.pct=0.5)

VlnPlot(Integrated, features=c("ICOSLG","LTA"))
VlnPlot(Integrated, features=c("RNA.weight"), pt.size=0)
FeaturePlot(Integrated, features=c("CD4"), reduction="wnn.umap", min.cutoff = 5)
FeaturePlot(Integrated, features=c("percent.mt", "nFeature_RNA"), reduction="wnn.umap", min.cutoff = 5)

DimPlot(Integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, 
        pt.size = 0.2) &NoLegend()

saveRDS(Integrated, "CITEseq Integrated 24.12.22.rds")


