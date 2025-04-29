
Myeloid <- subset(Integrated, idents=c("Monocyte/Macrophage/mDC","Mast","Neutrophil","Cycling", "pDC"))

DefaultAssay(Myeloid) <-"RNA"
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

Myeloid <- CellCycleScoring(Myeloid, s.features = s.genes, g2m.features = g2m.genes,
                            set.ident = TRUE)

Myeloid[["RNA"]]  <- split(Myeloid[["RNA"]], f = Myeloid$orig.ident)

Myeloid <- ScaleData(Myeloid)
Myeloid <- FindVariableFeatures(Myeloid, nfeatures=4000)
Myeloid <- RunPCA(Myeloid)

ElbowPlot(Myeloid, ndims=50, reduction="pca")

Myeloid <- IntegrateLayers(object = Myeloid, 
                           method = RPCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.rpca")

Myeloid[["RNA"]] <- JoinLayers(Myeloid[["RNA"]])
Myeloid <- FindNeighbors(Myeloid, reduction = "integrated.rpca", dims = 1:50)
Myeloid <- FindClusters(Myeloid, resolution = 1)
Myeloid <- RunUMAP(Myeloid, dims = 1:50, reduction = "integrated.rpca", reduction.name = "umap.rna")

DimPlot(Myeloid, reduction = "umap.rna", pt.size = 0.2, label=TRUE) &NoLegend()


DefaultAssay(Myeloid) = "dsbADT"

Myeloid[["dsbADT"]] <- JoinLayers(Myeloid[["dsbADT"]])
Myeloid[["dsbADT"]] <- split(Myeloid[["dsbADT"]], f = Myeloid$orig.ident)

Myeloid <- ScaleData(Myeloid, features = rownames(Myeloid@assays$dsbADT@features))
VariableFeatures(Myeloid) <- rownames(Myeloid@assays$dsbADT@features)
Myeloid <- RunPCA(Myeloid, reduction.name = 'apca')

Myeloid <- IntegrateLayers(object = Myeloid, 
                           method = RPCAIntegration, 
                           orig.reduction = "apca",
                           new.reduction = "integrated.adt.rpca",
                           features = rownames(Myeloid@assays$dsbADT@features))

Myeloid[["dsbADT"]] <- JoinLayers(Myeloid[["dsbADT"]])

Myeloid <- FindNeighbors(Myeloid, reduction = "integrated.adt.rpca", dims = 1:30)
Myeloid <- FindClusters(Myeloid, resolution = 0.5)
Myeloid <- RunUMAP(Myeloid, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="adt.umap")

DimPlot(Myeloid, reduction = "adt.umap", pt.size = 0.2, label=TRUE) &NoLegend()

#WNN integration
DefaultAssay(Myeloid) = "dsbADT"

Myeloid <- FindMultiModalNeighbors(
  Myeloid, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:50, 1:30), modality.weight.name = "RNA.weight")

Myeloid <- RunUMAP(Myeloid, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Myeloid <- FindClusters(Myeloid, graph.name = "wsnn", algorithm = 3, resolution = seq(from = 0.1, 
                                                                                      to = 1.0, 
                                                                                      by = 0.1))

Idents(Myeloid) <- "wsnn_res.1"
DimPlot(Myeloid, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2) +theme_void() &NoLegend()

DefaultAssay(Myeloid) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Myeloid, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

DefaultAssay(Myeloid) <-"dsbADT"
Protein.markers.1 <- FindAllMarkers(Myeloid, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)


#Resolution 1:7 T, 9 low quality, 10 B/Plasma, 16 Endothelial doublets

Myeloid <- subset(Myeloid, idents=c("7","9", "10","16"),invert=TRUE)
#Run from FindMultiModalNeighbors


adtNameData = read_excel("~/Desktop/CITEseq/TotalSeq A Names.xlsx")
rownames(Myeloid@assays$dsbADT@layers$counts) = adtNameData$FinalNames
rownames(Myeloid@assays$dsbADT@layers$data) = adtNameData$FinalNames
rownames(Myeloid@assays$dsbADT@features) = adtNameData$FinalNames
rownames(Myeloid@assays$ADT@counts) = adtNameData$FinalNames
rownames(Myeloid@assays$ADT@data) = adtNameData$FinalNames


Idents(Myeloid) <- "wsnn_res.1"

Myeloid <- RenameIdents(Myeloid, "13"="Neutrophil","3"="Neutrophil","1"="Neutrophil",
                        "0"="Mast","6"="Mast",
                        "9"="Monocyte","4"="Monocyte",
                        "5"="Monocyte/Macrophage", "2"="Macrophage","12"="ASDC",
                        "7"="cDC2", "14"="mregDC","8"="cDC1","10"="pDC",
                        "11"="Cycling")

levels(Myeloid)<-c("Neutrophil","Mast","Monocyte","Monocyte/Macrophage","Macrophage","ASDC","cDC2","mregDC", "cDC1" , "pDC", "Cycling")   

Myeloid$MyeloidClusters <- Myeloid@active.ident
Idents(Myeloid) <- "MyeloidClusters"

Myeloid$FineClusters <- Myeloid@active.ident
Idents(Myeloid) <- "FineClusters"

revMyeloidlevels <- rev(levels(Myeloid))
levels(Myeloid) <-revMyeloidlevels
Myeloid$revMyeloidClusters <-Myeloid@active.ident
Idents(Myeloid) <-"revMyeloidClusters"
revRNAFeatures <-rev(RNAFeatures)
revADTFeatures <-rev(ADTFeatures)

DimPlot(Myeloid, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2) +theme_void() &NoLegend()

DefaultAssay(Myeloid) <-"RNA"
RNAFeatures <- c(
  "CDC20","UBE2C","TOP2A", #Cycling
  "TCL1A","PTCRA","RGS7","COL26A1","JCHAIN", "PACSIN1","NIBAN3", #pDC
  "KCND3","DNASE1L3","BATF3", "ENPP1", "NCALD", "CLEC9A","XCR1", #cDC1
  "CCL22","CD200", "CCL19", "LAD1","CCR7","LAMP3",#mregDC
  "CCR6","CD1E","CD1C","FCER1A", "CLEC10A", #CDC2
  "CDH1","CD2","NKG7","PLAC8",  "SIGLEC6","CD5","AXL", #ASDC
  "CD209","IGF1", "PDGFC", "RNASE1","FOLR2","CD163L1","SELENOP", "C1QC",#Mac 2
  "C3","IL10","CX3CR1","KMO", #Mac 1
  "MARCO","S100A12","EREG", "CD300E","NLRP3", "VCAN", #Monocyte
  "IL1RL1","TPSAB1","CPA3",  #Mast
  "CXCR2","ADGRG3","FCGR3B"#Neutrophil
)

RNAdotplot <- DotPlot(Myeloid, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

RNAdotplot <- DotPlot(Myeloid, features = revRNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
                                   axis.title.x=element_blank(),
                                   axis.title.y=element_blank())
RNAdotplot


DefaultAssay(Myeloid) <-"dsbADT"
ADTFeatures <- c(
  #Cycling
  "CXCR3_p","CD62L_p","CD123_p",   #pDC
  "CD13_p","CD141_p", #cDC1
  "CSF2RA_p", "CD83_p","CD86_p",# mregDC
  "CD2_p","TACI_p","CD1c_p", #CDC2
  "CD5_p", "CD11a_p", #ASDC
  "CX3CR1_p", "PD-1_p", "CD74_p","CD54_p", "Siglec-1_p","CD40_p","SIRPα_p","NRP-1_p","CD32_p", "CD107a_p","CD64_p",  #Mac 2
  "TREM-1_p", "CD244_p", "CCR2_p", "CLEC12A_p", "CD101_p","Siglec-7_p", #Mac 1
  "CD55_p", "CD24_p", "CD93_p", #Monocyte
  "LAIR-1_p","DPP4_p", "CD33_p", "CD22_p","CD69_p","IgE_p","FCER1A_p",  #Mast
  "CR1_p", "LOX-1_p", "CD16_p","CD15_p","CD66_p", #Neutrophil"
  "CD11c_p","CD14_p","CD11b_p","HLA-DR_p"
)

ADTdotplot <- DotPlot(Myeloid, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

ADTdotplot <- DotPlot(Myeloid, features = revADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + 
  #coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 90, hjust=1),
                                   axis.title.x=element_blank(),
                                   axis.title.y=element_blank())

ADTdotplot

saveRDS(Myeloid, "Myeloid 24.12.22.rds")


Myeloid2 <- subset(Myeloid, idents=c("Neutrophil","Mast","Cycling"), invert=TRUE)

avgexp<- AggregateExpression(Myeloid2, return.seurat = T, 
                             group.by = c("MyeloidClusters"))

DefaultAssay(avgexp) <-"ADT"
dittoHeatmap(avgexp,c(ADTfeatures),
             order.by="MyeloidClusters",
             scaled.to.max = FALSE,
             column_order=c("Monocyte","Macrophage 1","Macrophage 2","ASDC","cDC2", "mregDC","cDC1", "pDC"),
             cluster_cols = FALSE,
             cluster_rows = FALSE,
             show_colnames = TRUE,
             annotation_legend=FALSE,
             heatmap_legend_param = list(title = "Expression", color_bar = "continuous"),
             complex=TRUE) #ComplexHeatmap






