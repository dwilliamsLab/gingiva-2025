
Idents(Integrated)<-"CellTypes"

TNK <- subset(Integrated, idents=c("T/NK", "Cycling"))

DefaultAssay(TNK) = "RNA"

#TNK[["RNA"]] = JoinLayers(TNK)

counts<- TNK[["RNA"]]$counts
TNK[["RNAsub"]] <- subset(TNK[["RNA"]], features=setdiff(rownames(TNK[["RNA"]]),"IGKC"))
DefaultAssay(TNK) <-"RNAsub"

TNK[["RNA"]] <- NULL
TNK[['RNA']] = TNK[['RNAsub']]
DefaultAssay(TNK) <-"RNA"
TNK[["RNAsub"]] <- NULL

TNK[["RNA"]] = split(TNK[["RNA"]], f = TNK$orig.ident)

TNK = NormalizeData(TNK)
TNK = FindVariableFeatures(TNK, nfeatures=4000)
TNK = ScaleData(TNK)
TNK = RunPCA(TNK, dims = 1:50)

ElbowPlot(TNK, ndims=50, reduction="pca")

TNK = IntegrateLayers(object = TNK, 
                       method =RPCAIntegration,  
                       orig.reduction = "pca", 
                       new.reduction = "integrated.rpca")

TNK[["RNA"]] = JoinLayers(TNK[["RNA"]])
TNK = FindNeighbors(TNK, reduction = "integrated.rpca", dims = 1:50)
TNK = FindClusters(TNK, resolution = 1)
TNK = RunUMAP(TNK, dims = 1:50, reduction = "integrated.rpca", reduction.name = "umap.rna")

DimPlot(TNK, reduction = "umap.rna", pt.size = 0.2, label=TRUE) &NoLegend()


DefaultAssay(TNK) = "dsbADT"

TNK[["dsbADT"]] = JoinLayers(TNK[["dsbADT"]])
TNK[["dsbADT"]] = split(TNK[["dsbADT"]], f = TNK$orig.ident)

VariableFeatures(TNK) = rownames(TNK@assays$dsbADT@features)
TNK = ScaleData(TNK, features = rownames(TNK@assays$dsbADT@features))
TNK = RunPCA(TNK, reduction.name = 'apca')

TNK = IntegrateLayers(object = TNK, 
                       method = RPCAIntegration, 
                       orig.reduction = "apca",
                       new.reduction = "integrated.adt.rpca",
                       features = rownames(TNK@assays$dsbADT@features))

TNK[["dsbADT"]] = JoinLayers(TNK[["dsbADT"]])

TNK = FindNeighbors(TNK, reduction = "integrated.adt.rpca", dims = 1:30)
TNK = FindClusters(TNK, resolution = 0.5)
TNK = RunUMAP(TNK, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="umap.adt")

DimPlot(TNK, reduction = "umap.adt", pt.size = 0.2, label=TRUE) &NoLegend()



# WNN integration
DefaultAssay(TNK) = "dsbADT"

TNK = FindMultiModalNeighbors(
  TNK, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:50, 1:30), modality.weight.name = "RNA.weight")

TNK = RunUMAP(TNK, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
TNK = FindClusters(TNK, graph.name = "wsnn", algorithm = 2, resolution = seq(from = 0.1, 
                                                                               to = 1.5, 
                                                                               by = 0.1))

Idents(TNK) <- "wsnn_res.0.8"
DimPlot(TNK, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2) &NoLegend()

VlnPlot(TNK, feature="RNA.weight")

DefaultAssay(TNK) <-"RNA"
RNA.markers.1 <- FindAllMarkers(TNK, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)



DefaultAssay(TNK) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(TNK, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)


#Resolution 1: 16 Plasma, 17 non-specific cycling, 18 myeloid,  19 neutrophil, 20 endothelial, 21 B doublets
TNK <- subset(TNK, idents=c("16","17","18","19","20","21"), invert=TRUE)

#Run from FindMultiModal

#Resolution 1.5: #5, 17, 18 low quality
TNK <- subset(TNK, idents=c("5","17","18"), invert=TRUE)
#Run from FindMultiModal


Idents(TNK) <- "wsnn_res.1.5"
TNK <- RenameIdents(TNK, "15"="CD4 T naive",
                     "0"="CD4 T CM",
                     "18"="CD4 TH17", "4"="CD4 IL17A+ FOXP3+","17"="CD4 Treg",
                     "7"="CD4 TPH", "13"="CD4 Follicular Reg",
                     
                    
                     "14"="CD8 T naive","5"="CD8 TRM",
                     
                     "1"="CD8 T 1", "2"="CD8 T 1",
                     "3"="CD8 T 2","8"="CD8 T 2", "10" ="CD8 T 2",
                     "6"="CD8 HELIOS+",
                     
                     "16"="MAIT","19"="T gd",
                     
                     "12"="CD8 TEMRA","11"="CD16+ NK", "20"="NK","9"="ILC")

TNK <- RenameIdents(TNK, "CD8 TEMRA"="T/NK mix", "CD8 HELIOS+"="CD8 T 3")
TNK <- RenameIdents(TNK, "CD8 T 1"="CD8 T", "CD8 T 2"= "CD8 T TNF+", "CD8 T 3"="CD8 T IKZF2+")
TNK <- RenameIdents(TNK, "CD4 TPH"="CD4 TFH", "CD4 Follicular Reg"="CD4 TFreg")

levels(TNK) <- c("CD4 T naive","CD4 TCM","CD4 TH17","CD4 IL17A+ FOXP3+","CD4 Treg","CD4 TFH", "CD4 TFreg", "CD8 T naive","CD8 TRM", "CD8 T", 
                  "CD8 T TNF+","CD8 T IKZF2+","MAIT","T gd","T/NK mix","CD16+ NK", "NK", "ILC")               
    
TNK$TNKClusters <- TNK@active.ident
Idents(TNK)<-"TNKClusters"

TNK$FineClusters <- TNK@active.ident
Idents(TNK)<-"FineClusters"

DimPlot(TNK, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2) +theme_void() &NoLegend()

CD8Clusters<- FindMarkers(TNK, ident.1="CD8 T 1", ident.2="CD8 T 2", min.pct=0.3)

saveRDS(TNK, "TNK 24.12.22.rds") 



DefaultAssay(TNK) <- "RNA"
RNAFeatures<-c(
  "KIT", "SPINK2","GAB1", #ILC
  "ADGRG3","KIR2DL4","KRT86", #NK
  "S1PR5","FCGR3A","KLRF1", #CD16+ NK
  "ADGRG1", "CX3CR1","FGFBP2", #CD8/NK mix (TEMRA)
  "TRDC","XCL1","KLRC2", #T gd
  "NCR3","IL4I1","SLC4A10",#MAIT
  "CADM1","PDGFD","IKZF2", #CD8 T 3 (HELIOS+)
  "TNF","IFNG","IFNG-AS1","EGR1","JUN", #CD8 T 2
   "CRTAM","CCL4","GZMK", #CD8 T
  "GZMB", "ITGAE","GZMH", #CD8 TRM
  "KLRC1","ITGA1","NELL2", #CD8 naive
  "CD200", "CDK14", "FHIT","TAFA2","PDE7B","HDAC9","RBMS3","SESN3", #CD4 T follicular reg
 "CXCR5","CXCL13","BCL6","BTBD11","ZBTB10","TNFSF8","CMTM8", "TSHZ2", #TPH
 "FANK1","HPGD","TNFRSF4","IL2RA","FOXP3", #Treg
  "IL23R","IL7",'TNFRSF18',"ZC2HC1A", "PRDM1","IL1R1","CCR6", #TH17 FOXP3
  "IL17A","IL17F","IL22","AHRR","PPARG","ADAM19",     #CD4 TH17 
   "TRAF1", "SVIL","GNAQ","ANXA1","IL7R", #CD4 CM
  "CCR7","SELL","LEF1", #CD4 naive
  "CD4","CD8A","CD3D"

)


RNAdotplot <- DotPlot(TNK, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot


revTNK <- rev(levels(TNK))
levels(TNK) <-revTNK
TNK$revTNKClusters <-TNK@active.ident
Idents(TNK) <-"revTNKClusters"
revRNAFeatures <-rev(RNAFeatures)
revADTFeatures <-rev(ADTFeatures)

RNAdotplot <- DotPlot(TNK, features = revRNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

RNAdotplot



DefaultAssay(TNK) <- "dsbADT"
ADTFeatures<-c(
 "SIRPα_p", "CD161_p","CD11c_p", #ILC
 "ITGB7_p", "CD56_p", "NCR1_p","CD244_p", #NK
  "NKp80_p","Siglec-7_p", #CD16+ NK
  "GPR56_p","SLAMF1_p", "CD16_p","CD43_p", #CD8 T/NK (TEMRA) 
   #Tgd 
    #MAIT
  "KLRG1_p", #CD8 T HELIOS
 "CD142_p","CD155_p", 'CD9_p', "CD62P_p", "CD49b_p", "CD57_p","CD84_p","CD314_p", #CD8 T
  "CD11a_p", "CD49a_p", #CD8 TRM
  "CD107a_p", "CD226_p","CD101_p", #CD8 T naive
  "CD62L_p","CD54_p","CD200_p","CD69_p",#CD4 Follicular Reg
 "CD40L_p", "CD49d_p", "CD32_p", "CD38_p","PD-1_p", #TPH
 "CTLA4_p","CCR4_p", "IL-2Rα_p","CD27_p","HLA-DR2_p","CD58_p", #Treg
 "LAIR-1_p","CD151_p",  "CD146_p", "CD7_p","CD82_p","CD5_p","CD28_p","ENTPD1_p","CD71_p", "HLA-A_p", "ICOS_p", #CD4 IL17A FOXP3
  "CD103_p", "DPP4_p",#CD4 TH17 
"CD73_p","TIGIT_p","HVEM_p","CD127_p","CD52_p", #CD4 T CM
   #CD4 T NAIVE
 "CD4_p","CD8_p","CD3_p"
)

ADTdotplot <- DotPlot(TNK, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

revADTFeatures<-rev(ADTFeatures)

ADTdotplot <- DotPlot(TNK, features = revADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


ADTdotplot

levels(CD4)<-c("CD4 Treg","CD4 Follicular Reg","CD4 T naive", "CD4 TCM", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 TPH")
DefaultAssay(CD4) <- "dsbADT"
ADTFeatures<-c("CD57", "PD-1","CD45RA","IL2RA","CD127")
ADTdotplot <- DotPlot(CD4, features = ADTFeatures, col.min=0.5, dot.min=0.3, scale=TRUE) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

#Tonsil Atlas
DefaultAssay(TNK)<- "dsbADT"
FeaturePlot(TNK, features=c("PTPRC2-PROT", "CD28-PROT"), min.cutoff=1, reduction="wnn.umap") #Naive

DefaultAssay(TNK)<- "RNA"
FeaturePlot(TNK, features=c("IL21"), reduction="wnn.umap")

VlnPlot(TNK, features=c("PRDM1","BCL6"))

FeaturePlot(TNK, features=c("CXCR5", "TCF7","IL6ST","ICOS", "TOX2","BTLA","CXCL13"), reduction="wnn.umap") #Tfh
FeaturePlot(TNK, features=c("TNFRSF4","CD200"), reduction="wnn.umap") #Tfh OXO (TNFRSF4) - T follicular regulatory -
FeaturePlot(TNK, features=c("FCRL3", "CLNK","LEF1","RBMS3", "SESN3"), reduction="wnn.umap") #T follicular regulatory --

DefaultAssay(TNK) <-"RNA"
avgexp <- AggregateExpression(TNK, return.seurat = T, 
                              group.by = c("TNKClusters"))

CD4 <- subset(TNK, idents=c("CD4 T naive", "CD4 TCM", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 Treg", "CD4 TPH","CD4 Follicular Reg"))
adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(CD4@assays$dsbADT@layers$counts) = adtNameData$FinalNames_no_p
rownames(CD4@assays$dsbADT@layers$data) = adtNameData$FinalNames_no_p
rownames(CD4@assays$dsbADT@features) = adtNameData$FinalNames_no_p
rownames(CD4@assays$ADT@counts) = adtNameData$FinalNames_no_p
rownames(CD4@assays$ADT@data) = adtNameData$FinalNames_no_p

DefaultAssay(CD4) <-"dsbADT"
CD4avgexp <- AggregateExpression(CD4, return.seurat = T, 
                              group.by = c("TNKClusters"))
CD4avgexp = ScaleData(CD4avgexp, features = rownames(CD4avgexp@assays$dsbADT@features))

DefaultAssay(avgexp) <-"RNA"
dittoHeatmap(avgexp,c("BCL6","PRDM1"),
             order.by="TNKClusters",
             scaled.to.max = FALSE,
             column_order=c("CD4 T naive", "CD4 Resting", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 Treg", "CD4 TPH",
                            "CD4 Follicular Reg", "CD8 T naive", "CD8 TRM", "CD8 T 1", "CD8 T 2", "CD8 HELIOS+", "CD8 TEMRA",
                            "MAIT", "T gd", "CD16+ NK", "NK", "ILC"),
             #gaps_row = c(4, 6, 10, 14),
             #gaps_col = c(7,13),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

DefaultAssay(CD4avgexp) <-"RNA"
dittoHeatmap(CD4avgexp,c("PRDM1","BCL6","FOXP3","IL2RA","TNFRSF4","LEF1","TCF7","SELL","RBMS3","SESN3","FCRL3"),
             order.by="TNKClusters",
             scaled.to.max = FALSE,
             column_order=c("CD4 T naive", "CD4 TCM", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 TPH","CD4 Treg","CD4 Follicular Reg"),
             gaps_row = c(2, 5, 8),
             #gaps_col = c(2,5),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap


DefaultAssay(CD4avgexp) <-"RNA"
dittoHeatmap(CD4avgexp,c("PRDM1","BCL6","IL21","CXCR5","POU2AF1","CD200", "SH2D1A","CXCL13","BTLA","ICOS","PDCD1","IL21","MAF",
                         "FOXP3","IL2RA","TNFRSF4","LEF1","TCF7","SELL","RBMS3","SESN3","FCRL3"),
             order.by="TNKClusters",
             scaled.to.max = FALSE,
             column_order=c("CD4 T naive", "CD4 TCM", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 TPH","CD4 Treg","CD4 Follicular Reg"),
             gaps_row = c(2, 12, 15),
             #gaps_col = c(2,5),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap.colors = colorRampPalette(c("white","white","lightgray", "red", "yellow"))(100),
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

DefaultAssay(CD4avgexp) <-"dsbADT"
dittoHeatmap(CD4avgexp,
             c("IL-2Rα","CD27","CD28","TIGIT","ICOS","PD-1","CD200","CD69","CD38","CD95","CD127","CD161","CD195"), #IL2RA: CD25, CD195: CCR5
             #order.by="CD4",
             scaled.to.max = FALSE,
             column_order=c("CD4 T naive", "CD4 TCM", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 TPH","CD4 Treg","CD4 Follicular Reg"),
             gaps_row = c(5,8),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap

dittoHeatmap(CD4avgexp,"CD88")

DefaultAssay(CD4avgexp) <-"ADT"
dittoHeatmap(CD4avgexp,c("CD127","IL2RA","CD45RA","PD-1","CD40","ICOS"), #IL2RA: CD25, CD195: CCR5
             #order.by="CD4",
             scaled.to.max = FALSE,
             column_order=c("CD4 Treg","CD4 Follicular Reg","CD4 T naive", "CD4 TCM", "CD4 TH17", "CD4 IL17A+ FOXP3+", "CD4 TPH"),
             gaps_row = c(4),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap



#Azimuth
library(SeuratDisk)
TNK[["RNA"]] <- as(object = TNK[["RNA"]], Class = "Assay")
SaveH5Seurat(TNK, "TNK.h5seurat")
#Upload on Azimuth web browser. Download mappings tsv, convert to csv and rename
mappings <- read.csv("tnk_azimuth_pred.csv", header=TRUE, row.names=1)

TNK <- AddMetaData(TNK, mappings, col.name="predicted.celltype.l1")

DimPlot(TNK,  pt.size = 0.2, label=TRUE, reduction="wnn.umap", group.by="predicted.celltype.l1", repel=TRUE) 
DimPlot(TNK,  pt.size = 0.2, label=TRUE, reduction="wnn.umap", group.by="TNKClusters", repel=TRUE) 



DefaultAssay(CD4) <-"RNA"
RNA.markers.1 <- FindAllMarkers(CD4, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)



DefaultAssay(CD4) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(CD4, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)
