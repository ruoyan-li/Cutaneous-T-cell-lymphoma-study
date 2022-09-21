library(dplyr)
library(Seurat)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
library(sctransform)
library(sceasy)
library(reticulate)

h5ad_file <- "/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL_Tcell_raw_no_nan.h5ad"
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                      outFile='/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL_Tcell_raw_no_nan.rds')

Tcell <- readRDS('/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL_Tcell_raw_no_nan.rds')
tcr <- read.csv('/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL1-8-tcr-pairs-4chain.csv')
#tcr$sample <- rep("CTCL", nrow(tcr))

inter <- intersect(rownames(Tcell@meta.data), rownames(tcr))
Tcell <- subset(Tcell, cells = inter)
tcr <- tcr[inter, ]
Tcell$TCR_pair <- tcr$TRA_pair

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

TCRgenes <- read.table("/home/jovyan/farm/tm/Tcells/rmTCRgene.list")
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize", scale.factor = 10000)
Tcell <- FindVariableFeatures(Tcell, selection.method = "vst")
Tcell <- CellCycleScoring(Tcell, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Tcell <- ScaleData(object = Tcell, vars.to.regress = c("percent_mito", "S.Score", "G2M.Score"))
Tcell <- ScaleData(object = Tcell, vars.to.regress = c("percent_mito"))
keep_gene_list <- setdiff(VariableFeatures(object = Tcell), as.vector(TCRgenes[,1]))
Tcell <- RunPCA(object = Tcell, features = keep_gene_list)

Tcell <- FindNeighbors(object = Tcell, dims = 1:30)
Tcell <- FindClusters(object = Tcell, resolution = 0.6)
Tcell <- RunUMAP(object = Tcell, dims = 1:30)

DimPlot(object = Tcell, reduction = "umap",pt.size = 0.01,group.by = 'donor_id')
DimPlot(object = Tcell, reduction = "umap",pt.size = 0.01,group.by = 'Tissue')
saveRDS(Tcell, file="/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL_Tcell.processed.rds")

Tcell <- readRDS("/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL_Tcell.processed.rds")
tcr <- read.csv('/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL1-8-tcr-pairs-4chain.csv')
#tcr$sample <- rep("CTCL", nrow(tcr))
inter <- intersect(rownames(Tcell@meta.data), rownames(tcr))
#Tcell <- subset(Tcell, cells = inter)
tcr <- tcr[inter, ]
Tcell$TCR_pair <- tcr$TCR_pair

################## 1 ######################
Tcell$clonotype <- rep("Other_TCRs",nrow(Tcell@meta.data))
Tcell$clonotype[which(Tcell$TCR_pair=="CAVRDLLTGGGNKLTF:CSARTGSNQPQHF")] <- 
  "Clonotype_donor1"

Tcell$clonotype[which(Tcell$TCR_pair=="None:CSARTGSNQPQHF" | 
                      Tcell$TCR_pair=="CAVRDLLTGGGNKLTF:None")] <- "Clonotype_donor1_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CAAGLRTRNAGNMLTF:CASRASSNSGRASYEQYF")] <- 
  "Clonotype_donor2"
Tcell$clonotype[which(Tcell$TCR_pair=="None:CASRASSNSGRASYEQYF" |
                      Tcell$TCR_pair=="CAAGLRTRNAGNMLTF:None")] <-"Clonotype_donor2_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CAVTGNQFYF:CSARTGGYGYTF")] <- 
  "Clonotype_donor3"
Tcell$clonotype[which(Tcell$TCR_pair=="None:CSARTGGYGYTF" |
                      Tcell$TCR_pair=="CAVTGNQFYF:None"  )] <- "Clonotype_donor3_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CAMVKGGSQGNLIF:CASSVSGSNTEAFF")] <- 
  "Clonotype_donor4"
Tcell$clonotype[which(Tcell$TCR_pair=="CAMVKGGSQGNLIF:None" |
                        Tcell$TCR_pair=="None:CASSVSGSNTEAFF"  )] <- "Clonotype_donor4_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CASSPSTQSNEKLFF:CAASDSWGKLQF")] <- 
  "Clonotype_donor5"
Tcell$clonotype[which(Tcell$TCR_pair=="CASSPSTQSNEKLFF:None" |
                      Tcell$TCR_pair=="None:CAASDSWGKLQF" )] <- "Clonotype_donor5_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CASSQAFGPPHGTIYF:CAAYMNSGYSTLTF")] <- 
  "Clonotype_donor6"
Tcell$clonotype[which(Tcell$TCR_pair=="CASSQAFGPPHGTIYF:None" |
                        Tcell$TCR_pair=="None:CAAYMNSGYSTLTF" )] <- "Clonotype_donor6_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CASSYRGPYNEQFF:CALSGTHQAGTALIF")] <- 
  "Clonotype_donor7"
Tcell$clonotype[which(Tcell$TCR_pair=="CASSYRGPYNEQFF:None" |
                        Tcell$TCR_pair=="None:CALSGTHQAGTALIF" )] <- "Clonotype_donor7_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="CASSLVMGVANTEAFF:CAVRGNYGQNFVF")] <- 
  "Clonotype_donor8"
Tcell$clonotype[which(Tcell$TCR_pair=="CASSLVMGVANTEAFF:None" |
                        Tcell$TCR_pair=="None:CAVRGNYGQNFVF" )] <- "Clonotype_donor8_sc"

Tcell$clonotype[which(Tcell$TCR_pair=="None:None")] <- "No_TCR"


################## 2 ######################
Tcell$clonotype <- rep("Other_TCRs",nrow(Tcell@meta.data))

cd1 <- union(grep('CAVRDLLTGGGNKLTF', Tcell$TCR_pair),
             grep('CSARTGSNQPQHF', Tcell$TCR_pair))
Tcell$clonotype[cd1] <- "Clonotype_donor1"

cd2 <- union(grep('CAAGLRTRNAGNMLTF', Tcell$TCR_pair),
             grep('CASRASSNSGRASYEQYF', Tcell$TCR_pair))
Tcell$clonotype[cd2] <- "Clonotype_donor2"

cd3 <- union(grep('CAVTGNQFYF', Tcell$TCR_pair),
             grep('CSARTGGYGYTF', Tcell$TCR_pair))
Tcell$clonotype[cd3] <- "Clonotype_donor3"

cd4 <- union(grep('CAMVKGGSQGNLIF', Tcell$TCR_pair),
             grep('CASSVSGSNTEAFF', Tcell$TCR_pair))
Tcell$clonotype[cd4] <- "Clonotype_donor4"

cd5 <- union(grep('CASSPSTQSNEKLFF', Tcell$TCR_pair),
             grep('CAASDSWGKLQF', Tcell$TCR_pair))
Tcell$clonotype[cd5] <- "Clonotype_donor5"

cd6 <- union(grep('CASSQAFGPPHGTIYF', Tcell$TCR_pair),
             grep('CAAYMNSGYSTLTF', Tcell$TCR_pair))
Tcell$clonotype[cd6] <- "Clonotype_donor6"

cd7 <- union(grep('CASSYRGPYNEQFF', Tcell$TCR_pair),
             grep('CALSGTHQAGTALIF', Tcell$TCR_pair))
Tcell$clonotype[cd7] <- "Clonotype_donor7"

cd8 <- union(grep('CASSLVMGVANTEAFF', Tcell$TCR_pair),
             grep('CAVRGNYGQNFVF', Tcell$TCR_pair))
Tcell$clonotype[cd8] <- "Clonotype_donor8"

Tcell$clonotype[which(Tcell$TCR_pair=="nan:nan:nan:nan" |
                      Tcell$TCR_pair=="None:None:None:None ")] <- "No_TCR"
##############################

#Col_3 <- c("Tomato","black","CornflowerBlue","black",
#           "Goldenrod2","black", "LightSlateBlue","black",
#           "DarkKhaki","black","Cyan2",'black',
#           "DarkSeaGreen","black","PaleVioletRed1","black",
#           "SlateGray2","Honeydew2")
Col_3 <- c("Tomato","CornflowerBlue","Goldenrod2", "LightSlateBlue",
           "DarkKhaki","Cyan2","DarkSeaGreen","PaleVioletRed1",
           "SlateGray2","Honeydew2")
DimPlot(object = Tcell, reduction = "umap",pt.size = 0.2, group.by = "clonotype",cols = Col_3)
DimPlot(object = Tcell, reduction = "umap",pt.size = 0.2, group.by = "Tissue",cols = Col_3)

### tumour cell label
#ll <- read.csv('/home/jovyan/farm/CTCL/malignant_marker/donor1-8_tumorcell_list.1.xls', sep=",")
#inter <- intersect(rownames(Tcell@meta.data), rownames(ll))
#Tcell$tu <- rep("Other", nrow(Tcell@meta.data))
#Tcell$tu[inter] <- "tumour_cell"
#DimPlot(object = Tcell, reduction = "umap",pt.size = 0.2, group.by = "tu",cols = c("grey","darkblue"))

################# Inter-tumour DEGs ####################
tu_list <- read.table("/home/jovyan/farm/CTCL/malignant_marker/donor1-8_tumorcell_list.xls",header=T, check.names=F)
Tcell_tu <- subset(x=Tcell, cells = as.vector(tu_list$cell))

Tcell_tu@active.ident <- factor(Tcell_tu$donor_id)
markers <- FindAllMarkers(object = Tcell_tu, only.pos = TRUE, 
                          min.pct = 0.1, logfc.threshold = 0.25)
x <- markers %>% group_by(cluster) %>% top_n(n =50, wt = avg_log2FC)
top100 <- markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)

cluster.averages <- AverageExpression(Tcell_tu, return.seurat = TRUE)
DoHeatmap(cluster.averages, features = top100$gene, size = 3, 
          draw.lines = FALSE)+  scale_fill_gradientn(colors = c("steelblue2", "white", "tomato2"))

################## TCR ################
tcr_data <- read.csv('/lustre/scratch117/cellgen/team205/rl20/CTCL/object-new/CTCL1-8-tcr-pairs.1.csv')
sort(table(tcr_data$TCR[which(tcr_data$donor_id=='CTCL1')]))

