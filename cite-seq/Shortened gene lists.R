#All Cells

Idents(Integrated) <-"CellTypes"

DefaultAssay(Integrated) <- "RNA"
RNAFeatures<-c(
  "MKI67", #Cycling
  "KRT14", #Epithelial
  "CCL21", #LEC
  "PLVAP", #Endothelial
  "ACTA2", #SMC  
  "PDGFRA", #Fibroblast
  "CSF1R","LYZ", #Mono/Mac/DC
  "FCGR3B", #Neutrophil
  "CPA3", #Mast
  "LILRA4", #pDC
  "IGHA1", #Plasma
  "MS4A1", #B
  "GZMA", "CD3D" #T/NK
)


RNAdotplot <- DotPlot(Integrated, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot


DefaultAssay(Integrated) <- "dsbADT"
ADTFeatures<-c(
  "CD71_p", #Epithelial
  "PDPN_p", #LEC
  "CD62P_p", #Endothelial
  'CD146_p', #SMC
   'CD13_p', #Fibroblast
   "CD14_p", #Monocyte
   "CD66_p", #Neutrophil
  "FCER1A_p", #Mast
  "CLEC4C_p", #pDC
  "CD38_p", #Plasma 
  "CD20_p", #B
  "CD244_p","CD5_p" #TNK
)

#ADTdotplot <-DotPlot(Integrated, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.51, 
#                     scale=TRUE) + coord_flip() +
#  scale_color_gradient2(low="grey90", high="black") +
#  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot <- DotPlot(Integrated, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low = "white", mid = "white",high="black", midpoint = 1.2) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot




#TNK

DefaultAssay(TNK) <- "RNA"
RNAFeatures<-c(
  "KIT", #ILC
  "KLRC1","KLRF1", #NK
  "FCGR3A", #CD16+ NK
  "CX3CR1","GZMH", #CD8 TEMRA
  "TRDC", #T gd
  "IL23R",#MAIT
  "IKZF2", #HELIOS+ CD8
  "TNF", #CD8 T 2
   #CD8 T
  "ITGAE", "GZMK", #CD8 TRM
   #CD8 naive
  "CD200",#CD4 T follicular reg
   "CXCR5","BCL6",#TPH
   "FOXP3", #Treg
   #TH17 FOXP3
  "IL17A","PRDM1",#CD4 TH17 
   #CD4 CM
  "SELL","CCR7", "IL7R", #CD4 naive
  "CD4","CD8A","CD3D"
)


RNAdotplot <- DotPlot(TNK, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot



DefaultAssay(TNK) <- "dsbADT"
ADTFeatures<-c(
  "SIRPα_p", #ILC
  "CD56_p", #NK
  "NKp80_p", #CD16+ NK
  "CD16_p", #CD8 TEMRA 
   #Tgd 
"CD161_p",   #MAIT
"KLRG1_p",  #CD8 T HELIOS
"CD49b_p",'CD9_p', #CD8 T
 "CD49a_p", #CD8 TRM
 "CD101_p", #CD8 T naive
 "CD200_p","CD69_p",  #CD4 Follicular Reg
  "PD-1_p", #TPH
  "CTLA4_p","IL-2Rα_p", #Treg
  "ICOS_p", #CD4 IL17A FOXP3
  "DPP4_p", #CD4 TH17 
 "CD127_p", #CD4 TCM. CD127=IL-7ra
   #CD4 T NAIVE 
 "CD38_p","HLA-DR_p", "CD27_p","CD28_p", "CD45RA_p","CD45RO_p","CD8_p","CD4_p","CD3_p"
)

ADTdotplot <- DotPlot(TNK, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot

#ADTdotplot <- DotPlot(TNK, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
#                      scale=TRUE) + coord_flip() +
#  scale_color_gradient2(low = "white", mid = "white",high="black", midpoint = 0.5) +
#  theme(axis.text.x = element_text(angle = 45, hjust=1))







Idents(BPlasma) <-"BPlasmaClusters"

DefaultAssay(BPlasma) <-"RNA"
RNAFeatures <- c(
  "IGKC","IGHG4", #5 Plasma 3
   #2 Plasma 2
  "ITGA8", #0 Plasma 1
   "IGHA2","IGHA1",#4 Histone IGHA
  "MKI67","UBE2QL1", #Plasmablast
  "FCRL4","ITGAX","IL21R",  #13 Atypical
  "TNF", "BCL2A1","EGR1",  #B Memory 2
  "CD40","CR1",#B Memory 1
  "IL4R","IGHD","IGHM", #11 Naive
  "CD38","CD19","HLA-DRA","MS4A1"
)


RNAdotplot <-DotPlot(BPlasma, features = RNAFeatures, col.min=0.1,  dot.min=0.1, 
                     scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

DefaultAssay(BPlasma) <-"dsbADT"
ADTFeatures <- c(
  "ITGB7_p","PD-L1_p","SLAMF7_p","CD44_p", #0 CD44 Plasma
  "CD86_p","CD71_p", #Plasmablast
  "TACI_p","PD-1_p",  #Atypical
  #Memory 2
  "CR1_p","BAFF-R_p", #Memory 1
  "CD72_p",  #Naive
  "CD38_p","CD27_p","CD19_p","HLA-DR_p", "CD20_p"
)

ADTdotplot <- DotPlot(BPlasma, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot




DefaultAssay(Myeloid) <-"RNA"
RNAFeatures <- c(
  "MKI67", #Cycling
  "TCL1A", #pDC
  "XCR1", #cDC1
  "LAMP3", "CCR7",#mregDC
  "CD1C", #CDC2
  "SIGLEC6", "AXL", #ASDC
  "CD163L1",#Mac 2
  "CD14",#Mac 1
   "VCAN", #Monocyte
  "TPSAB1",  #Mast
  "CXCR2", #Neutrophil
  "HLA-DRA"
)

RNAdotplot <-DotPlot(Myeloid, features = RNAFeatures, col.min=0.1, scale.by="radius", dot.min=0.1, 
                     scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot


DefaultAssay(Myeloid) <-"dsbADT"
ADTFeatures <- c(
  #Cycling
  "CD123_p",   #pDC
  "CD141_p",  #cDC1
  "CD86_p","CD83_p",# mregDC
  "CD1c_p", #CDC2
  "CD5_p", #ASDC
  "CX3CR1_p","CD32_p", #Mac 2
  "CD101_p",#Mac 1
  "CD93_p", #Monocyte
  "FCER1A_p",  #Mast
  "CD66_p", #Neutrophil"
  "HLA-DR_p","CD11c_p","CD14_p","CD11b_p"
)

#ADTdotplot <- DotPlot(Myeloid, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
#                      scale=TRUE) + coord_flip() +
#  scale_color_gradient2(low="grey90", high="black") +
#  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot <- DotPlot(Myeloid, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.3, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low = "white", mid = "white",high="black", midpoint = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot





Idents(Stromal) <- "StromalClusters"
DefaultAssay(Stromal) <-"RNA"

RNAFeatures <- c(
  "MKI67", #Cycling
  "IL7","CCL21","MMRN1",#LEC
  
  "GJA5", #Endo 5
   "RGCC", #Endo 4
   "BTNL9", #Endo 3
  "CSF3","ENPP2", #Endo 2
  "SELE", "ACKR1", #Endo 1
  

  "KCNJ8","ABCC9", #SMC 2
  "RGS6","MYH11", #SMC 1
  
  "STC1", "ENPP1",   #Fib 6
  "SFRP4","COL11A1", #Fib 5
  "CNTN1","COMP",  #Fib 4
  "TMEM132C", "IGF1", #Fib 3
  "ICAM1", "LTBR","VCAM1","CCL19", "CXCL13", #Fib 2
  "CXCL8","CXCL6", "CXCL1",  #Fib 1
  "PECAM1","MCAM", "PDGFRA","THY1"
  
)

RNAdotplot <-DotPlot(Stromal, features = RNAFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                     scale=TRUE) + coord_flip() +
  scale_color_gradient2(low="grey90", high="black") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

RNAdotplot

DefaultAssay(Stromal) <-"dsbADT"
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

ADTdotplot <- DotPlot(Stromal, features = ADTFeatures, col.min=0, scale.by="size", dot.min=0.1, 
                      scale=TRUE) + coord_flip() +
  scale_color_gradient2(low = "white", mid = "white",high="black", midpoint = 0.3) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ADTdotplot




