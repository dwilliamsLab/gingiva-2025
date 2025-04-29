
Stromal <- subset(Integrated, ident=c("Fibroblast", "Endothelial", "SMC","LEC")) 

Stromal[["RNA"]] <- split(Stromal[["RNA"]], f = Stromal$orig.ident)

DefaultAssay(Stromal) = "RNA"

Stromal <- ScaleData(Stromal)
Stromal <- FindVariableFeatures(Stromal, nfeatures=4000)
Stromal <- RunPCA(Stromal)

ElbowPlot(Stromal, ndims=50, reduction="pca")

Stromal <- IntegrateLayers(object = Stromal, 
                           method = RPCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.rpca")

Stromal[["RNA"]] <- JoinLayers(Stromal[["RNA"]])
Stromal <- FindNeighbors(Stromal, reduction = "integrated.rpca", dims = 1:50)
Stromal <- FindClusters(Stromal, resolution = 1)
Stromal <- RunUMAP(Stromal, dims = 1:50, reduction = "integrated.rpca", reduction.name = "umap.rna")

DimPlot(Stromal, reduction = "umap.rna", pt.size = 0.2, label=TRUE) &NoLegend()

DefaultAssay(Stromal) <- "dsbADT"
adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(Stromal@assays$dsbADT@layers$counts) = adtNameData$FinalNames_no_p
rownames(Stromal@assays$dsbADT@layers$data) = adtNameData$FinalNames_no_p
rownames(Stromal@assays$dsbADT@features) = adtNameData$FinalNames_no_p
rownames(Stromal@assays$ADT@counts) = adtNameData$FinalNames_no_p
rownames(Stromal@assays$ADT@data) = adtNameData$FinalNames_no_p


Stromal[["dsbADT"]] <- split(Stromal[["dsbADT"]], f = Stromal$orig.ident)

Stromal <- ScaleData(Stromal, features = rownames(Stromal@assays$dsbADT@features))
VariableFeatures(Stromal) <- rownames(Stromal@assays$dsbADT@features)
Stromal <- RunPCA(Stromal, reduction.name = 'apca')

Stromal <- IntegrateLayers(object = Stromal, 
                           method = RPCAIntegration, 
                           orig.reduction = "apca",
                           new.reduction = "integrated.adt.rpca",
                           features = rownames(Stromal@assays$dsbADT@features))

Stromal[["dsbADT"]] <- JoinLayers(Stromal[["dsbADT"]])

Stromal <- FindNeighbors(Stromal, reduction = "integrated.adt.rpca", dims = 1:30)
Stromal <- FindClusters(Stromal, resolution = 0.5)
Stromal <- RunUMAP(Stromal, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="adt.umap")

DimPlot(Stromal, reduction = "adt.umap", pt.size = 0.2, label=TRUE) &NoLegend()

# WNN integration
Stromal <- FindMultiModalNeighbors(
  Stromal, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:50, 1:30), modality.weight.name = "RNA.weight")

Stromal <- RunUMAP(Stromal, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Stromal <- FindClusters(Stromal, graph.name = "wsnn", algorithm = 3, resolution = seq(from = 0.1, 
                                                                                      to = 1.1, 
                                                                                      by = 0.1))

Idents(Stromal) <- "wsnn_res.2"
DimPlot(Stromal, reduction = 'wnn.umap', label = TRUE, repel = TRUE) +theme_void() &NoLegend()
DimPlot(Stromal, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2, split.by="orig.ident") &NoLegend()


DefaultAssay(Stromal) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Stromal, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)



DefaultAssay(Stromal) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(Stromal, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

#Resolution 1
#13, 14 low quality, 15 T, 17 T, 18 B, 19 Mixed Immune, 20 Mast, 21 Neutrophil, 22 Neuron, 24 Platelet, 25 B


Stromal <- subset(Stromal, idents=c("13", "14","15","17","18","19","20","21","22","24","25"), invert=TRUE)
#Run from FindMultiModal

#Resolution 1
#19 Neutrophil
Stromal <- subset(Stromal, idents=c("19"), invert=TRUE)

#Rerun from FindMultiModal up to Resolution 1.1

Idents(Stromal) <- "wsnn_res.1.1"
Stromal <- RenameIdents(Stromal,
                        "2"="1","10"="0","3"="0","17"="14")
levels(Stromal) <- c("1","15","8","13","5","6","9","7","0","11","14","4","12","16","18") 
Stromal <- RenameIdents(Stromal, "1"="Fibroblast 1", "15"="Fibroblast 2", "8" = "Fibroblast 3", "13"="Fibroblast 4", "5"="Fibroblast 5", "6"="Fibroblast 6",
                        "9" = "SMC 1", "7" = "SMC 2", "0" = "Endothelial 1", "11"= "Endothelial 2", "14"="Endothelial 3", "4" = "Endothelial 4", "12"="Endothelial 5",
                        "16" = "LEC", "18" = "Cycling")
Stromal$StromalClusters <-Stromal@active.ident
Idents(Stromal) <- "StromalClusters"
DimPlot(Stromal, reduction = 'wnn.umap', label = TRUE, repel = TRUE) +theme_void() &NoLegend()

Stromal$FineClusters <-Stromal@active.ident
Idents(Stromal) <- "FineClusters"


Idents(Stromal) <- "StromalClusters"
revStromal <- rev(levels(Stromal))
levels(Stromal) <-revStromal
Stromal$revStromalClusters <-Stromal@active.ident
Idents(Stromal) <-"revStromalClusters"

revRNAFeatures <-rev(RNAFeatures)
revADTFeatures <-rev(ADTFeatures)



DefaultAssay(Stromal) <-"RNA"
RNAFeatures <- c(
  "MKI67", "BIRC5","CENPF",#18 #Cycling
  "TNFRSF11A","IL7","MRC1","PDPN","CCL21","MMRN1",#16
  
  "TSPAN2","ATP2A3","PLLP", "GJA5","SEMA3G", #12. Large artery
  "KCNQ1","KCNK5","ADAMTSL2","RGCC", #4
  "CHRM3","RBP7","BTNL9", #14. Capillary
  "SELL","SH3TC1","TMC1", "PLA1A","LAMC2","CSF3","LAMP3","ENPP2","COTL1","CSF2RB", #11
  "CLU","PLAT","ADAMTS9", "SELE","ACKR1", #0
  
  "PAMR1","KCNJ8","ABCC9","GUCY1A2", #7
  "SOX6", "MYOCD","MUSTN1","RGS6","ACTA2","MYH11", #9
 

   "SLC19A1", "KRT17","F13A1","CP","ENPP1","BMP4","STC1",  #6
 "PRSS35","FNDC1","SFRP4","OGN", "CORIN","LRRC15", "ASPN","COL11A1", #5
 "POU3F3","HPSE2","FOXP2","LRFN5","PTGER3", "CNTN1","COMP", #13
   "SLIT2", "TMEM132C","FGL2","IGF1","ABCA8","CFD", #8
 "CXCL10","CXCL9","ADAMDEC1","IRF8","G0S2","IL6","BIRC3", "CD24","RBP5","CH25H","TNFSF13B", "DAZL","CCL19","CXCL13", #15
"CD82","IL24","CXCL8","CXCL6","LBP", "CA12","LAMC3","CXCL1" #3
)

RNAdotplot <- DotPlot(Stromal, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

RNAdotplot <- DotPlot(Stromal, features = revRNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

RNAdotplot


DefaultAssay(Stromal) <-"dsbADT"
ADTFeatures <- c(
  "CD73_p","PDPN_p", #LEC
 "CD109_p", "CD93_p", "CD151_p","CD155_p","CD58_p","CD123_p", #12
  "CD61_p", "CD49b_p", #Endothelial 4
  "CD49f_p", "NRP-1_p", #14
  "CD137_p","IFNGR1_p","CSF2RB_p",  #11
  "CD82_p","CD23_p", "DPP4_p","CD62P_p", "HLA-DR_p", #0
  
  "CD55_p", 'CD146_p',"CD29_p", "CD49a_p", #9
  
  #13
  #6
 "CR2_p","CR1_p","CD27_p","TACI_p", "CD49d_p", "CD44_p","CD9_p", #5
"ICOSLG_p", "ENTPD1_p", #8
"CD54_p","CD24_p","CD142_p", "CD81_p", #15
 "GPR56_p","CD14_p", "PD-L1_p", "CD40_p","CD95_p" #1
)

ADTdotplot <- DotPlot(Stromal, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

ADTdotplot <- DotPlot(Stromal, features = revADTFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

ADTdotplot

saveRDS(Stromal, "Stromal 25.1.7.rds")



Stromal2 <- subset(Stromal, idents="Cycling Stromal", invert=TRUE)
adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(Stromal2@assays$dsbADT@layers$counts) = adtNameData$FinalNames_no_p
rownames(Stromal2@assays$dsbADT@layers$data) = adtNameData$FinalNames_no_p
rownames(Stromal2@assays$dsbADT@features) = adtNameData$FinalNames_no_p
rownames(Stromal2@assays$ADT@counts) = adtNameData$FinalNames_no_p
rownames(Stromal2@assays$ADT@data) = adtNameData$FinalNames_no_p

avgexp <- AggregateExpression(Stromal2, return.seurat = T, 
                              group.by = c("StromalClusters"))


Stromal3 <- subset(Stromal, idents=c("Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4", "Fibroblast 5", "Fibroblast 6"))
adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(Stromal3@assays$dsbADT@layers$counts) = adtNameData$FinalNames_no_p
rownames(Stromal3@assays$dsbADT@layers$data) = adtNameData$FinalNames_no_p
rownames(Stromal3@assays$dsbADT@features) = adtNameData$FinalNames_no_p
rownames(Stromal3@assays$ADT@counts) = adtNameData$FinalNames_no_p
rownames(Stromal3@assays$ADT@data) = adtNameData$FinalNames_no_p


VlnPlot(Stromal3, features=c("CR1","CR2","CXCL12","CXCL13","CCL19","CCL21","ICAM1","VCAM1","MADCAM1","TNFRSF13B","PDPN","LTBR"))

DefaultAssay(Stromal3) <- "dsbADT"
Fibavgexp <- AggregateExpression(Stromal3, return.seurat = T, 
                                 group.by = c("StromalClusters"))

DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("CR1","CR2","CXCL12","CXCL13","CCL19","CCL21","ICAM1","VCAM1","TNFRSF13B","PDPN","LTBR"),
             order.by="StromalClusters",
             scaled.to.max = FALSE,
             column_order=c("Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4", 
                            "SMC 1", "SMC 2",
                            "Endothelial 1","Endothelial 2","Endothelial 3","Endothelial 4","Endothelial 5",
                            "LEC"),
             #gaps_row = c(4, 6, 10, 14),
             #gaps_col = c(7,13),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

DefaultAssay(Stromal) <-"RNA"
dittoHeatmap(Fibavgexp,c("CXCL13","CCL19","ICAM1","VCAM1","TNFSF13B","PDPN","LTBR","TNFRSF11A", "CR1","CR2"),
             order.by="StromalClusters",
             scaled.to.max = FALSE,
             column_order=c("Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4", "Fibroblast 5", "Fibroblast 6"),
             #gaps_row = c(4, 6, 10, 14),
             #gaps_col = c(7,13),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap



DefaultAssay(Fibavgexp) <-"ADT"
dittoHeatmap(Fibavgexp,c("CD95","CD40","PDPN","CD24","CD81","ICOSLG","CR2", "CR1","CD23","CD274"),
             order.by="StromalClusters",
             scaled.to.max = FALSE,
             column_order=c("Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4", "Fibroblast 5", "Fibroblast 6"),
             #gaps_row = c(4, 6, 10, 14),
             #gaps_col = c(7,13),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap







DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("CXCL13","CXCL12","CCL19","CR2","CLU","MADCAM1","TNFRSF11A","ACTA2","CD34","PI16","FBLN1"),
             order.by="StromalClusters",
             scaled.to.max = FALSE,
             column_order=c("Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4", 
                            "SMC 1", "SMC 2",
                            "Endothelial 1","Endothelial 2","Endothelial 3","Endothelial 4","Endothelial 5",
                            "LEC"),
             #gaps_row = c(4, 6, 10, 14),
             #gaps_col = c(7,13),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap


DefaultAssay(avgexp) <-"ADT"
dittoHeatmap(avgexp,c("CR1","CR2","CD23", "ICOSLG","PDPN"),
             order.by="StromalClusters",
             scaled.to.max = FALSE,
             column_order=c("Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4", 
                            "SMC 1", "SMC 2",
                            "Endothelial 1","Endothelial 2","Endothelial 3","Endothelial 4","Endothelial 5",
                            "LEC"),
             #gaps_row = c(4, 6, 10, 14),
             #gaps_col = c(7,13),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) 



library(readxl)
library(UCell)

hHEV<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="hHEV", col_names = FALSE)
hHEV <- hHEV$...1
aHEV<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="aHEV", col_names = FALSE)
aHEV <- aHEV$...1
VBEC<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="VBEC", col_names = FALSE)
VBEC <- VBEC$...1
TRC<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="TRC", col_names = FALSE)
TRC <- TRC$...1
MRC<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="MRC", col_names = FALSE)
MRC <- MRC$...1
FDC<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="FDC", col_names = FALSE)
FDC <- FDC$...1

AbeSignatures <- list(hHEV = hHEV, aHEV = aHEV, VBEC = VBEC, TRC = TRC, MRC = MRC, FDC = FDC)

DefaultAssay(Stromal) <- "RNA"
Stromal <- AddModuleScore_UCell(Stromal, features=AbeSignatures)
FeaturePlot(Stromal, feature=c("aHEV_UCell", "hHEV_UCell", "VBEC_UCell"), min.cutoff= 0.3, reduction="wnn.umap") 
FeaturePlot(Stromal, feature=c("TRC_UCell", "MRC_UCell", "FDC_UCell"), min.cutoff = 0.3, reduction="wnn.umap")

FeaturePlot_scCustom(Stromal, features=c("MRC_UCell", "FDC_UCell"), colors_use= viridis_plasma_dark_high,
                     na_cutoff=0.31, alpha_exp = 0.5, reduction="wnn.umap") *theme_void() *ggtitle(NULL)


FeaturePlot_scCustom(Fibro, features=c("MRC_UCell", "FDC_UCell"), colors_use= viridis_plasma_dark_high,
                     na_cutoff=0.3, alpha_exp = 0.5, reduction="wnn.umap")

VlnPlot(Stromal, feature=c("MRC_UCell", "FDC_UCell"))



FDCgenes <- c("CXCL13","CCL19","ICAM1","VCAM1","MADCAM1", "TNFSF13B","PDPN","LTBR")

HEV<-read_xlsx("/Users/fraserdr/Desktop/GeneList.xlsx", sheet="HEV", col_names = FALSE)
HEV <-HEV$...1


gene.sets <- list(FDCgenes = FDCgenes, HEV = HEV)

Stromal <- AddModuleScore_UCell(Stromal, features=gene.sets)
FeaturePlot(Stromal, feature=c("HEV_UCell"), min.cutoff = 0.25, reduction="wnn.umap") 
FeaturePlot(Stromal, feature=c("FDCgenes_UCell"), min.cutoff = 0.3, reduction="wnn.umap") 

FeaturePlot_scCustom(Stromal, features=c("HEV_UCell"), colors_use= viridis_plasma_dark_high,
                     na_cutoff=0.37, alpha_exp = 0.6, reduction="wnn.umap") +theme_void() +ggtitle(NULL)

FeaturePlot_scCustom(Stromal, features=c("FDCgenes_UCell"), colors_use= viridis_plasma_dark_high,
                     na_cutoff=0.4, alpha_exp = 0.5, reduction="wnn.umap") +theme_void() +ggtitle(NULL)

FeaturePlot_scCustom(Endo, features=c("HEV_UCell"), colors_use= viridis_plasma_dark_high,
                     na_cutoff=0.4, alpha_exp = 0.5, reduction="wnn.umap") +theme_void() +ggtitle(NULL)

df<- get_expression_data(Stromal)

umap_cor<- pbmc3k@reductions$umap@cell.embeddings  %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "cell")

df<- left_join(df, umap_cor)

head(df)


FDC2list <- list(c("CXCL13","PDPN","CCL19","ICAM1","VCAM1","IL7"))
Stromal <- AddModuleScore(
  object = Stromal,
  features = FDC2list,
  name = "FDC2list")

FeaturePlot(Stromal, feature=c("FDC2list1"), min.cutoff = 0, reduction="wnn.umap") 

Idents(Stromal) <-"StromalClusters"
Fibro <- subset(Stromal, idents=c("Fibroblast 1", "Fibroblast 2", "Fibroblast 3", "Fibroblast 4", "Fibroblast 5", "Fibroblast 6"))
Endo <- subset(Stromal, idents=c("Endothelial 1", "Endothelial 2", "Endothelial 3", "Endothelial 4", "Endothelial 5"))


Idents(Fibro) <- "StromalClusters"
DefaultAssay(Fibro) <-"RNA"

RNAFeatures <- c(
  
  "STC1", "ENPP1",   #Fib 6
  "SFRP4","COL11A1", #Fib 5
  "CNTN1","COMP",  #Fib 4
  "TMEM132C", "IGF1", #Fib 3
  "ICAM1", "LTBR","VCAM1","CCL19", "CXCL13", #Fib 2
  "CXCL8","CXCL6", "CXCL1",  #Fib 1
  "PECAM1","MCAM", "PDGFRA","THY1"
  
)

RNAdotplot <-DotPlot(Fibro, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                     scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot


DefaultAssay(Fibro) <-"dsbADT"
ADTFeatures <- c(
  "CD73_p","PDPN_p", #LEC
  "CD58_p", "CD155_p", #Endotheloal 5
  #Endothelial 4
  "CD49f_p", "CD49b_p",  #Endothelial 3
  "CD137_p", "IFNGR1_p", #Endothelial 2
  "CSF2RB_p", "CD62P_p",  #Endothelial 1
  
  
  #SMC 2
  "CD29_p", "CD49a_p", #SMC 1
  
  #Fib 6
  "CR2_p","CR1_p", "CD49d_p", #Fib 5
  #Fib 4
  "ENTPD1_p","ICOSLG_p", #Fib 3
  "CD24_p","CD81_p",  #Fibroblast 2
  "PD-L1_p", "CD40_p","CD95_p",  #Fibroblast 1
  "CD31_p", "CD146_p", "CD13_p"
)

ADTdotplot <- DotPlot(Fibro, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low = "white", mid = "white",high="black", midpoint = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot









DefaultAssay(Fibro) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Fibro, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)



DefaultAssay(Fibro) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(Fibro, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)




Idents(Stromal) <-"CellTypes"
Endo <- subset(Stromal, idents="Endothelial")

Idents(Endo)<-"StromalClusters"
DefaultAssay(Endo) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Endo, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)



DefaultAssay(Endo) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(Endo, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)





###########
#FIRST RUN

Idents(Stromal) <- "wsnn_res.1"
Stromal <- RenameIdents(Stromal, "3"="1","12"="7","15"="7","2"="0")
levels(Stromal) <- c("16","1","6","7","4","8","9","10","0","13","5","11","14") 

Stromal <- RenameIdents(Stromal, "2"="1","18"="12","17"="14","20"="0","3"="0","10"="0")

Idents(Stromal) <- "StromalClusters"

DefaultAssay(Stromal) <-"RNA"
RNAFeatures <- c(
  "MKI67", "BIRC5","CENPF",#16 #Cycling
  "TNFRSF11A","IL7","MRC1","PDPN","CCL21","MMRN1",#14
  
  "TSPAN2","ATP2A3","PLLP", "GJA5","SEMA3G", #Endothelial 5. Large artery
  "KCNQ1","KCNK5","ADAMTSL2","RGCC", #4
  "CHRM3","RBP7","BTNL9", #14. Capillary
  "HLA-DQA1", "DNASE1L3","CCL14", #3
  "SELL","SH3TC1","TMC1", "PLA1A","LAMC2","CSF3","LAMP3","ENPP2","COTL1","CSF2RB", #13
  
  
  "PAMR1","KCNJ8","ABCC9","GUCY1A2", #7
  "SOX6", "MYOCD","MUSTN1","RGS6","ACTA2","MYH11", #9
  
  
  "STXBP6", "PRSS35","FNDC1","SFRP4","OGN","LRRC15", "ASPN","COL11A1", #5 reticular?
  "SAA1","MMP3","SLC19A1", "KRT17","ENPP1","F13A1","CP","BMP4",  #6
  "LRFN5","HPSE2","CNTN1","FOXP2","COMP", #12
  "FGL2","TMEM132C","IGF1","ABCA8","CFD", #8
  "IRF8","G0S2","IL6","CD24","TNFSF13B", "DAZL","CCL19","CXCL13", #15
  "IL24","CXCL8","CXCL6","LAMC3","F3", "COL7A1","COL24A1"  #1 papillary?
  
)

RNAdotplot <- DotPlot(Stromal, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

Stromal <- RenameIdents(Stromal, "16"="Cycling Stromal", "1"="Fibroblast 1", "6"="Fibroblast 2","4"="Fibroblast 3","7"="Fibroblast 4",
                        "8"="SMC 1","9"="SMC 2","10"="Endothelial 1","0"="Endothelial 2","13"="Endothelial 3", "5"="Endothelial 4",
                        "11"="Endothelial 5", "14"="LEC")

levels(Stromal) <-c( "Fibroblast 1", "Fibroblast 2","Fibroblast 3","Fibroblast 4",
                     "SMC 1","SMC 2","Endothelial 1","Endothelial 2","Endothelial 3", "Endothelial 4",
                     "Endothelial 5", "LEC","Cycling Stromal")

Stromal$StromalClusters <-Stromal@active.ident
Idents(Stromal) <- "StromalClusters"

Stromal$FineClusters <-Stromal@active.ident
Idents(Stromal) <- "FineClusters"
DimPlot(Stromal, reduction = 'wnn.umap', label = TRUE, repel = TRUE) +theme_void() &NoLegend()

adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(Stromal@assays$dsbADT@layers$counts) = adtNameData$FinalNames
rownames(Stromal@assays$dsbADT@layers$data) = adtNameData$FinalNames
rownames(Stromal@assays$dsbADT@features) = adtNameData$FinalNames
rownames(Stromal@assays$ADT@counts) = adtNameData$FinalNames
rownames(Stromal@assays$ADT@data) = adtNameData$FinalNames


Idents(Stromal) <- "StromalClusters"

DefaultAssay(Stromal) <-"RNA"


RNAFeatures <- c(
  "MKI67", "BIRC5","CENPF",#16 #Cycling
  "TNFRSF11A","IL7", "MRC1","CCL21","MMRN1",#14
  
  "TSPAN2","ATP2A3","PLLP", "GJA5","SEMA3G", #Endothelial 5. Large artery
  "KCNQ1","KCNK5","RGCC", #Endothelial 4. Arteriole
  "CHRM3","RBP7","BTNL9", #Endothelial 3. Capillary
  "HLA-DQA1", "DNASE1L3","CCL14", #Endothelial 2. Vein
  "PLA1A","LAMC2", "SELL","PTAFR","CSF3","LAMP3", "CSF2RB", #Endothelial 1. Activated/HEV
  
  "SOX6", "MYOCD","MUSTN1","RGS6","MYH11", #SMC 2
  "PAMR1","KCNJ8","ABCC9","GUCY1A2", # SMC 1
  
  "GLDN","TNFSF13B","FGL2", "DAZL","NTRK2","ABCA8","CFD", #Fibroblast 4
  "SAA1","MMP3","SLC19A1", "KRT17","F13A1", #Fibroblast 3
  "PTGIS","STAC", "FNDC1","TENM3","SFRP4","LRRC15", "ASPN","CXCL14","COL11A1",  #Fibroblast 2
  "CCL19","CXCL3","IL6", "IL24","CXCL8","CXCL6","CA12","CHI3L2", "G0S2", "LBP","CXCL1","CXCL13"  #Fibroblast 1
  
)

RNAdotplot <- DotPlot(Stromal, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

DefaultAssay(Stromal) <-"dsbADT"
ADTFeatures <- c(
  "CD73_p","PDPN_p", #LEC
  "CD93_p", "CD112_p", "CD151_p","CD155_p","CD58_p","CD123_p", #Endotheilal 5
  "CD61_p",  "CD49b_p", #Endothelial 4
  "CD200_p", "CD49f_p",  #Endothelial 3
  "CD40L_p", "CD49d_p","CD23_p", "CD107a_p", "CD62P_p",  #Endothelial 2
  "CD137_p","PD-1_p","IFNGR1_p", "CSF2RB_p",  #Endothelial 1
  
  "CD146_p", "CD29_p","CD49a_p",  #SMC 2
  #SMC 1
  
  "HLA-DR2_p","HLA-DR_p", #Fibroblast 4
  "GPR56_p", #Fibroblast 3
  "CR2_p", "CR1_p","TACI_p","PD-L1_p","ICOSLG_p","CD44_p","CD9_p",  #Fibroblast 2
  "CD24_p", "CD40_p","CD142_p","CD95_p","CD81_p","CD54_p"  #Fibroblast 1
)

ADTdotplot <- DotPlot(Stromal, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

################
