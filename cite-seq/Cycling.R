#Subset cycling from top and run standard script
Cycling <- subset(Integrated, idents="Cycling")

Cycling[["RNA"]] = split(Cycling[["RNA"]], f = Cycling$orig.ident)

DefaultAssay(Cycling) <-"RNA"
Cycling = NormalizeData(Cycling)
Cycling = FindVariableFeatures(Cycling, nfeatures=2000)
Cycling = ScaleData(Cycling)
Cycling = RunPCA(Cycling, dims = 1:30)

Cycling = IntegrateLayers(object = Cycling, 
                      method =RPCAIntegration,  #CCA and RPCA has problems since some of orig.idents have few cells? 
                      orig.reduction = "pca", 
                      k.weight=40,
                      new.reduction = "integrated.rpca")

Cycling[["RNA"]] = JoinLayers(Cycling[["RNA"]])
Cycling = FindNeighbors(Cycling, reduction = "integrated.rpca", dims = 1:30)
Cycling = FindClusters(Cycling, resolution = 1)
Cycling = RunUMAP(Cycling, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rna")

DimPlot(Cycling, reduction = "umap.rna", pt.size = 0.2, label=TRUE) &NoLegend()


DefaultAssay(Cycling) = "dsbADT"

Cycling[["dsbADT"]] = split(Cycling[["dsbADT"]], f = Cycling$orig.ident)

VariableFeatures(Cycling) = rownames(Cycling@assays$dsbADT@features)
Cycling = ScaleData(Cycling, features = rownames(Cycling@assays$dsbADT@features))
Cycling = RunPCA(Cycling, reduction.name = 'apca')

Cycling = IntegrateLayers(object = Cycling, 
                      method = RPCAIntegration, 
                      orig.reduction = "apca",
                      new.reduction = "integrated.adt.rpca",
                      k.weight=40,
                      features = rownames(Cycling@assays$dsbADT@features))

Cycling[["dsbADT"]] = JoinLayers(Cycling[["dsbADT"]])

Cycling = FindNeighbors(Cycling, reduction = "integrated.adt.rpca", dims = 1:30)
Cycling = FindClusters(Cycling, resolution = 0.5)
Cycling = RunUMAP(Cycling, dims = 1:30, reduction = "integrated.adt.rpca", reduction.name="umap.adt")

DimPlot(Cycling, reduction = "umap.adt", pt.size = 0.2, label=TRUE) &NoLegend()



# WNN integration
Cycling = FindMultiModalNeighbors(
  Cycling, reduction.list = list("integrated.rpca", "integrated.adt.rpca"), 
  dims.list = list(1:30, 1:30), modality.weight.name = "RNA.weight")

Cycling = RunUMAP(Cycling, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Cycling = FindClusters(Cycling, graph.name = "wsnn", algorithm = 2, resolution = seq(from = 0.1, 
                                                                             to = 1, 
                                                                             by = 0.1))

Idents(Cycling) <- "wsnn_res.1"
DimPlot(Cycling, reduction = 'wnn.umap', label = TRUE, repel = TRUE, pt.size = 0.2) &NoLegend()


DefaultAssay(Cycling) <-"RNA"
RNA.markers.1 <- FindAllMarkers(Cycling, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.5)
RNA.significant.markers.1 <- RNA.markers.1[RNA.markers.1$p_val_adj < 0.2, ]

RNAmarkerTop10 <- RNA.significant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

DefaultAssay(Cycling) <-"dsbADT"

Protein.markers.1 <- FindAllMarkers(Cycling, 
                                    only.pos = TRUE, 
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.5)
Proteinsignificant.markers.1 <- Protein.markers.1[Protein.markers.1$p_val_adj < 0.2, ]

ProteinmarkerTop10 <- Proteinsignificant.markers.1 %>% 
  group_by(cluster) %>% 
  top_n(n=10, wt=avg_log2FC)

Idents(Cycling) <- "wsnn_res.1"
Cycling <- RenameIdents(Cycling, "0"="Plasma", "5"="B", "3"="CD4 T", "1"="CD8 T", "2"="Macrophage", "6"="DC","4"="Endothelial")
Cycling$Clusters <- Cycling@active.ident
Idents(Cycling)<- "Clusters"

saveRDS(Cycling, "Cycling 24.12.10.rds")

ImmuneCycling <- subset(Cycling, idents="Endothelial", invert=TRUE)
levels(ImmuneCycling)

ggplot(ImmuneCycling@meta.data, aes(x=orig.ident, fill=Clusters)) + geom_bar(position = "fill", color = "black", size=0.25) +
  theme_classic() +
  scale_fill_manual(values=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  theme(axis.text.x = element_text(size=0, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title = element_text(size=14)) +
  xlab("Patient") + ylab("Percentage") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = c('0', '25', '50', '75', '100')) 




