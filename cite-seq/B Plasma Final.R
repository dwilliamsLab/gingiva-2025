
BPlasma <- subset(Integrated, idents=c("B", "Plasma", "Cycling"))

DefaultAssay(BPlasma) = "RNA"

BPlasma[["RNA"]] = split(BPlasma[["RNA"]], f = BPlasma$orig.ident)

BPlasma = ScaleData(BPlasma)
BPlasma = FindVariableFeatures(BPlasma, nfeatures = 4000)
BPlasma = RunPCA(BPlasma, dims = 1:50)

ElbowPlot(BPlasma, ndims=50, reduction="pca")


BPlasma = IntegrateLayers(object = BPlasma, 
                          method =RPCAIntegration, 
                          orig.reduction = "pca", 
                          new.reduction = "integrated.rpca")

BPlasma = FindNeighbors(BPlasma, reduction = "integrated.rpca", dims = 1:50)
BPlasma = FindClusters(BPlasma, resolution = 1)
BPlasma = RunUMAP(BPlasma, dims = 1:50, reduction = "integrated.rpca", reduction.name="rna.umap")

DimPlot(BPlasma, reduction = 'rna.umap', group.by="RNA_snn_res.1", label = TRUE, repel = TRUE)

BPlasma[["RNA"]] = JoinLayers(BPlasma[["RNA"]])

DefaultAssay(BPlasma) = "dsbADT"

BPlasma[["dsbADT"]] = split(BPlasma[["dsbADT"]], f = BPlasma$orig.ident)

VariableFeatures(BPlasma) = rownames(BPlasma@assays$dsbADT@features)
BPlasma = ScaleData(BPlasma, features = rownames(BPlasma@assays$dsbADT@features))
BPlasma = RunPCA(BPlasma, reduction.name = 'apca')

BPlasma = IntegrateLayers(object = BPlasma, 
                          method = RPCAIntegration, 
                          orig.reduction = "apca",
                          new.reduction = "integrated.adt.rpca",
                          features = rownames(BPlasma@assays$dsbADT@features))

BPlasma[["dsbADT"]] = JoinLayers(BPlasma[["dsbADT"]])

BPlasma = FindMultiModalNeighbors(
  BPlasma, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")

BPlasma = RunUMAP(BPlasma, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
BPlasma = FindClusters(BPlasma, graph.name = "wsnn", algorithm = 3, resolution = seq(from = 0.1, 
                                                                                     to = 1.0, 
                                                                                     by = 0.1))


DimPlot(BPlasma, reduction = 'wnn.umap', label = TRUE, repel = TRUE)

Idents(BPlasma) <-"wsnn_res.0.4"
DefaultAssay(BPlasma) <-"RNA"
RNA.markers.1 <- FindAllMarkers(GC, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

DefaultAssay(BPlasma) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(BPlasma, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

#Resolution 1:
#11 Mast, 12 T/NK, 13 T/NK, 16 Endothelial, 17 Myeloid, 18 Endothelial, 19 T/NK, 20 Neutrophil doublets

BPlasma <- subset(BPlasma, idents=c("11","12","13","16","17","18","19","20"), invert=TRUE)

#Run from FindMultiModal.

#14 Fibroblast doublets, #12 low quality cells

BPlasma <- subset(BPlasma, idents=c("12","14"), invert=TRUE)

#Run from FindMultiModal



Idents(BPlasma) <-"wsnn_res.1"
BPlasma <- RenameIdents(BPlasma,"9"="8", "12"="10", "1"="0","3"="0", "2"="0","7"="0")
levels(BPlasma) <- c("11","6","8","13","10","4","0","5")
BPlasma <- RenameIdents(BPlasma, "11"="B Naive","6"="B Memory 1", "8"="B Memory 2","13"="B Atypical", 
                        "10"="Plasmablast", "4"="Plasma 1", "0" = "Plasma 2", "5"="Plasma 3")
BPlasma$BPlasmaClusters <- BPlasma@active.ident
Idents(BPlasma) <-"BPlasmaClusters"

DimPlot(BPlasma, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.5) +theme_void() +ggtitle(NULL) &NoLegend()


revBPlasmaClusters <-rev(levels(BPlasma))
levels(BPlasma) <-revBPlasmaClusters
BPlasma$revBPlasmaClusters <-BPlasma@active.ident
Idents(BPlasma) <-"revBPlasmaClusters"

revRNAfeatures <-rev(RNAFeatures)
revADTfeatures <-rev(ADTFeatures)


BPlasma$FineClusters <- BPlasma@active.ident 
Idents(BPlasma) <-"FineClusters"

DefaultAssay(BPlasma) <-"RNA"
RNAFeatures <- c(
  "IGKC","IGLC2","IGHG4", #5 Plasma 3
  "MERTK","ITGA8","FZD6","FAM13A","TBCEL","KANK1","DCC",
  "ARMH1","HSPA1A","HIST1H2BG","IGHA2","IGHA1", #4 Histone IGHA
  "RRM2","DIAPH3", "MKI67", "UBE2QL1", "KCNMA1", #10 Plasmablast
  "AICDA","CCL3","RGS13","BATF","LMO2","TRAF4","EBI3","TRAF1","MIR155HG", #Atypical GC-like 
  "LTA","IL2RB", "FCRL4","ATP8B4","TNFSF11","GSN", "HCK","IL21R", "ITGAX","DUSP4","CD86","FCRL5", #Atypical other
 "TNF","BCL2A1","EGR1",  #B Memory 2
"EPHA4", "ANK3","CR2",#B Memory 1
"CCR7","BCL6", "FCER2","TCL1A","IL4R","IGHD","IGHM"#11 Naive
  )


RNAdotplot <- DotPlot(BPlasma, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0, scale.min=1,
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

RNAdotplot <- DotPlot(BPlasma, features = revRNAfeatures, col.min=0, scale.by="size", dot.min=0, scale.min=1,
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

RNAdotplot

DefaultAssay(BPlasma) <-"dsbADT"
ADTFeatures <- c(
 "CD123_p","ITGB7_p","CR2_p","PD-L1_p","CD151_p","HVEM_p","ICOSLG_p","SLAMF7_p","CD44_p", #Plasma 2
  "CD18_p","CD86_p","CD11a_p","CD71_p", #Plasmablast
  
  "SIRPα_p","IFNGR1_p","TACI_p","CCR6_p","PD-1_p","CD52_p", "CD95_p", #13 Atypical
  "CD69_p", "CLEC12A_p", #Memory 2
  "IgM_p", "CD62L_p","BAFF-R_p","IGKC_p","CR1_p", #Memory 1
  "CD22_p","CD72_p", #Naive
  "CD27_p","CD38_p","CD19_p","HLA-DR_p","CD20_p"
)

ADTdotplot <- DotPlot(BPlasma, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

ADTdotplot <- DotPlot(BPlasma, features = revADTfeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

ADTdotplot



FeaturePlot(BPlasma, features=c("CD44","CR2","ICOSLG","CD274","CD123"), reduction="wnn.umap", min.cutoff=1)

saveRDS(BPlasma, "BPlasma with BCR 24.12.22.rds")




#Memory
DefaultAssay(BPlasma) <-"RNA"
FeaturePlot(BPlasma, feature=c("CD27","TNFRSF13B","CD38"), reduction="wnn.umap")
VlnPlot(BPlasma, feature=c("CD27","TNFRSF13B"), layer="counts")
DefaultAssay(BPlasma) <-"dsbADT"
FeaturePlot(BPlasma, feature=c("IgD","IgM","CD38","CD27","CD11c"), reduction="wnn.umap", min.cutoff=0) #CD27=TNFRSF13b
VlnPlot(BPlasma, feature=c("CD27","CD267"))


#TONSIL ATLAS
#GC
FeaturePlot(BPlasma, feature=c("BCL6","IRF8","MEF2B", "CXCR4","CD101","CD83","LMO2"), reduction="wnn.umap")
#GC Light Zone
FeaturePlot(BPlasma, feature=c("CD83","LMO2","BCL2A1","DUSP2"), reduction="wnn.umap")
#GC Dark Zone
FeaturePlot(BPlasma, feature=c("CXCR4","CD101","AICDA","FOXP1"), reduction="wnn.umap") #+CellCycle
#GC Light zone dark zone recomittment
FeaturePlot(BPlasma, feature=c("MYC","MIR155HG","NFKB1","BATF", "BCL2A1"), reduction="wnn.umap") #+CellCycle

FeaturePlot(BPlasma, features=c("CCND2", "MIR155HG","PSME2","BHLHE40", "PARVB","EBI3", "BLC2A1","LMO2", "GMDS"), min.cutoff=0, reduction="wnn.umap")


library(SeuratDisk)
BPlasma[["RNA"]] <- as(object = BPlasma[["RNA"]], Class = "Assay")
SaveH5Seurat(BPlasma, "BPlasma.h5seurat")
#Unload BPlasma and re-upload (don't save with converted RNA)
rm(BPlasma)
#Upload on Azimuth web browser. Download mappings tsv, convert to csv and rename
mappings <- read.csv("bplasma_azimuth_pred.csv", header=TRUE, row.names=1)

BPlasma <- AddMetaData(BPlasma, mappings, col.name="predicted.celltype.l1")

DimPlot(BPlasma,  pt.size = 0.2, label=TRUE, reduction="wnn.umap", group.by="predicted.celltype.l1", repel=TRUE) &NoLegend()
DimPlot(BPlasma,  pt.size = 0.2, label=TRUE, reduction="wnn.umap", group.by="Clusters", repel=TRUE) &NoLegend()





GCfeatures <- c("CXCR4","AICDA","FOXP1","MME",
                "CD83", "LMO2",
                "CD40","BCL2A1","MIR155HG","EBI3","TRAF1","NFKB1","RELB","MYC","BATF")
DotPlot(BPlasma, features = GCfeatures, col.min=0, dot.min=0.1, scale=FALSE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


library(dittoSeq)
GC.Commited <- c("CD69","EGR2","CCL3","MYC","CCND2")
Pre.GC <-c("MEF2B","RGS13")
DZ <- c("CXCR4","AICDA","FOXP1","MME")
LZ <- c("CD83", "LMO2")
S.PHASE <- c("POLA1","HIST1H4C")
G2M.PHASE <- c("MKI67","TOP2A","CDC20","CCNB1")
PROLIF <- c("STMN1","HMGB2")
REENTRY <-c("CD40","BCL2A1","MIR155HG","EBI3","TRAF1","NFKB1","NFKB2","RELB","MYC","BATF")



DefaultAssay(BPlasma) <-"RNA"
avgexp <- AggregateExpression(BPlasma, return.seurat = T, 
                              group.by = c("BPlasmaClusters"))

RNAavgexp <- avgexp$RNA

Idents(BPlasma) <- "CellTypes"
B <- subset(BPlasma, idents="B")
Idents(B) <-"Clusters"
B <- subset(B, idents="Plasma 5", invert=TRUE)
DefaultAssay(B) <-"dsbADT"

Bavgexp<- AggregateExpression(B, return.seurat = T, 
                            group.by = c("Clusters"))


dittoHeatmap(avgexp,c(GC.Commited, Pre.GC, DZ, LZ, S.PHASE, G2M.PHASE, PROLIF, REENTRY),
             order.by="BPlasmaClusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Atypical","Plasmablast","Plasma 1","Plasma 2","Plasma 3"),
             gaps_row = c(6,8,12,14,16,20,22),
             #gaps_col = c(5),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

DefaultAssay(Bavgexp) <-"RNA"
dittoHeatmap(Bavgexp,c("IGHD","IGHM","CR2","CD27","TNFRSF13B","ITGAX","MS4A1","CXCR5","BCL6"),
             order.by="BPlasmaClusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Atypical","Plasmablast","Plasma 1","Plasma 2","Plasma 3"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap


"BHLHE40", "PARVB","CCND2",     "PSME2". "HOPX"

#Tonsil Atlas Top GC Hits
DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("CCL3","RGS13","AICDA","LMO2", "MKI67","TOP2A","CCNB1","MIR155HG", "EBI3","TRAF1","TRAF4","BATF"),
             order.by="BPlasmaClusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Atypical","Plasmablast","Plasma 1","Plasma 2","Plasma 3"),
             gaps_row = c(2,4,7),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

#Tonsil Atlas Top GC Hits + King 2021 extra
DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("CCL3","RGS13","BHLHE40", "PARVB","MIR155HG", "EBI3","TRAF1","TRAF4","BATF","AICDA","LMO2","HOPX", "MKI67","TOP2A"),
             order.by="BPlasmaClusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Atypical","Plasmablast","Plasma 1","Plasma 2","Plasma 3"),
             gaps_row = c(9,12),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             data.out=FALSE,
             complex=TRUE) #ComplexHeatmap

DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("CD69", "EGR2","CCL3","MYC", "BATF","CCND2"),
             order.by="BPlasmaClusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Atypical","Plasmablast","Plasma 1","Plasma 2","Plasma 3"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

#King 2021 (https://www.science.org/doi/10.1126/sciimmunol.abe6291)
DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("CCND2","MIR155HG","PSME2","BHLHE40","PARVB","EBI3","BCL2A1","LMO2","GMDS",
                      "PRPSAP2","SERPINA9","MARCKSL1","CD27","CD38","BCL6","SUGCT","EZR","ISG20","AICDA",
                      "RASSF6","FRZB","HOPX","BTNL9","FGFR1"),
             order.by="BPlasmaClusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Atypical","Plasmablast","Plasma 1","Plasma 2","Plasma 3"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap


DefaultAssay(Bavgexp) <-"ADT"
dittoHeatmap(Bavgexp,c("IgD","IgM","CR2","CD27","CD267","CD11c","CD20","CXCR5","CD95","CD45RA","CD11b","CD127","CD38"),
             order.by="Clusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Memory 3","B Atypical"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

DefaultAssay(avgexp) <-"ADT"
dittoHeatmap(avgexp,c("CD20","IgD","IgM","CD23","CD24","CD27","CD38"),
             order.by="Clusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Memory 3","B Atypical", "Plasmablast","Plasma 1", "Plasma 2","Plasma 3","Plasma 4","Plasma 5"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmapDefaultAssay(avgexp) <-"ADT"

dittoHeatmap(avgexp,c("CD19","IgD","IgM","CD20","CD27","CD11c","CD38"),
             order.by="Clusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Memory 3","B Atypical", "Plasmablast","Plasma 1", "Plasma 2","Plasma 3","Plasma 4","Plasma 5"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

dittoHeatmap(Bavgexp,c("CD19","IgD","IgM","CD20","CD27","CD11c","CD38"),
             order.by="Clusters",
             scaled.to.max = FALSE,
             column_order=c("B Naive","B Memory 1","B Memory 2","B Memory 3","B Atypical"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap


#Manual Heatmap
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("blue","white", "#FD8A95"))

mat <- AverageExpression(BPlasma, assays = "RNA", layer="data",
                               features = c(DZ, LZ, S.PHASE, G2M.PHASE, REENTRY),
                               return.seurat = F, group.by = c("Clusters"))
mat <-as.matrix(mat$RNA)
mat <- t(apply(mat, 1, function(x) (x - min(x)) / diff(range(x))))

mat <- mat[c("CXCR4","AICDA","FOXP1","MME","CD83", "LMO2","MCM5","TYMS","CDCA7","DTL","MKI67","TOP2A","CENPF","CDC20","CD40","BCL2A1","MIR155HG","EBI3","TRAF1","NFKB1","RELB","MYC","BATF"),]

row_split = rep("DZ", 23)
row_split[5:6] = "LZ"
row_split[7:10] = "S Phase"
row_split[11:14] = "G2M Phase"
row_split[15:23] = "Reentry"

col_split = rep("B", 11)
col_split[6:11] = "Plasma"

Heatmap(mat,
        col=col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        row_split = row_split,
        column_split = col_split,
        show_heatmap_legend = TRUE,
        name="Expression",
        heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),)





B <- subset(BPlasma, idents=c("B Naive", "B Memory 1", "B Memory 2", "B Atypical"))
Plasma <- subset(BPlasma, idents=c("Plasmablast", "Plasma 1", "Plasma 2","Plasma 3"))

DefaultAssay(B) <-"RNA"
RNA.markers.1 <- FindAllMarkers(B, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

DefaultAssay(B) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(B, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)


DefaultAssay(B) <-"RNA"
RNAFeatures <- c(
  "IGKC","IGLC2","IGHG4", #5 Plasma 3
  "MERTK","ITGA8","FZD6","FAM13A","TBCEL","KANK1","DCC",
  "ARMH1","HSPA1A","HIST1H2BG","IGHA2","IGHA1", #4 Histone IGHA
  "RRM2","DIAPH3", "MKI67", "UBE2QL1", "KCNMA1", #10 Plasmablast
  
  "LMO2", "BATF", "IL2RB", "FCRL4", "EBI3","ATP8B4","TNFSF11","GSN", "HCK", "ITGAX","TRAF1", "IL21R", "DUSP4","CD86","MIR155HG","FCRL5",  #13 Atypical
  "TNF","BCL2A1","EGR1",  #B Memory 2
  "EPHA4", "ANK3","CR1", #B Memory 1
  "CCR7","BCL6", "FCER2","TCL1A","IL4R","IGHD","IGHM","BACH2" #11 Naive
)


RNAdotplot <- DotPlot(B, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

DefaultAssay(B) <-"dsbADT"
ADTFeatures <- c(
  "CD123_p","ITGB7_p","CR2_p","PD-L1_p","CD151_p","HVEM_p","ICOSLG_p","SLAMF7_p","CD44_p", #Plasma 2
  "CD18_p","CD86_p","CD11a_p","CD71_p", #Plasmablast
  
  "SIRPα_p","IFNGR1_p","TACI_p","CCR6_p","PD-1_p","CD52_p", "CD95_p", #13 Atypical
  "CD69_p", "CLEC12A_p", #Memory 2
  "IgM_p", "CD62L_p","BAFF-R_p","IGKC_p","CR1_p", #Memory 1
  "CD22_p","CD72_p", #Naive
  "CD27_p","CD38_p","CD19_p","HLA-DR_p","CD20_p"
)

ADTdotplot <- DotPlot(B, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot
