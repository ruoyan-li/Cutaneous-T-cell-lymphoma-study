library(sceasy)
library(reticulate)
library(ktplots)

h5ad_file <- "/lustre/scratch117/cellgen/team205/rl20/CTCL/cellphoneDB/CTCL_sub_nomore500_anno.h5ad"
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",
                      outFile='/lustre/scratch117/cellgen/team205/rl20/CTCL/cellphoneDB/CTCL_sub_nomore500_anno.rds')

seurat_object <- readRDS("/lustre/scratch117/cellgen/team205/rl20/CTCL/cellphoneDB/CTCL_sub_nomore500_anno.rds")
meta <- read.csv('/lustre/scratch117/cellgen/team205/rl20/CTCL/cellphoneDB/CTCL_sub_nomore500_anno_meta.csv')
seurat_object@meta.data <- as.data.frame(meta)
head(seurat_object@meta.data)

#### cellphoneDB output
pvals <- read.delim("/home/jovyan/farm/CTCL/cellphoneDB_3/out_CTCL_sub500_iter1000/pvalues.txt", check.names = FALSE)
means <- read.delim("/home/jovyan/farm/CTCL/cellphoneDB_3/out_CTCL_sub500_iter1000/means.txt", check.names = FALSE)
deconv <- read.delim("/home/jovyan/farm/CTCL/cellphoneDB_3/out_CTCL_sub500_iter1000/deconvoluted.txt", check.names = FALSE)

library(SingleCellExperiment)
library(reticulate)
library(ktplots)
ad=import('anndata')

adata = ad$read_h5ad(h5ad_file)
counts <- Matrix::t(adata$X)
row.names(counts) <- row.names(adata$var)
colnames(counts) <- row.names(adata$obs)
sce <- SingleCellExperiment(list(counts = counts), colData = adata$obs, rowData = adata$var)

######## B/plasma and tumour cells
plot_cpdb(cell_type1 = 'B/plasma', cell_type2 = 'tumourcell', scdata = seurat_object,
          idents = 'anno', means = means, pvals = pvals, 
          keep_significant_only=T) + small_guide()  + small_legend(keysize=1)

p <- plot_cpdb2(cell_type1 = 'B/plasma', cell_type2 = 'tumourcell',
                scdata = sce,
                idents = 'anno', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                desiredInteractions = list(
                  c('B/plasma', 'tumourcell'),
                  c('tumourcell', 'B/plasma')),
                node_group_colors = c(
                  "B/plasma" = "red",
                  "tumourcell" = "darkgreen"),
                keep_significant_only = TRUE,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
) + scale_edge_alpha(limits = c(0, 7))
p

p <- plot_cpdb3(cell_type1 = 'tumourcell', cell_type2 = 'B/plasma',
                scdata = sce,
                idents = 'anno', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = decon, # new options from here on specific to plot_cpdb3
                keep_significant_only = TRUE,
                standard_scale = TRUE,
                remove_self = TRUE,
)
p

############### F2/F3 and tumour cell
data <- read.table('/home/jovyan/farm/CTCL/DEG-condition-psudobulk/DEG-Fibro-filterby50_100.xls')
genes <- c(rownames(data))

plot_cpdb(cell_type1 = 'F2|F3', cell_type2 = 'tumourcell', scdata = seurat_object,
          idents = 'anno', means = means, pvals = pvals, 
          keep_significant_only=F, genes=genes) + small_guide()  + small_legend(keysize=1)

p <- plot_cpdb2(cell_type1 = 'F2|F3', cell_type2 = 'tumourcell',
                scdata = sce,
                idents = 'anno', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                desiredInteractions = list(
                  c('F2', 'tumourcell'),
                  c('tumourcell', 'F2'),
                  c('F3', 'tumourcell'),
                  c('tumourcell', 'F3')),
                node_group_colors = c(
                  "F2" = "dodgerblue3",
                  "F3" = "darkgreen",
                  "tumourcell" = "red"),
                keep_significant_only = T,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
                genes=genes
) + scale_edge_alpha(limits = c(0, 7))
p

p <- plot_cpdb3(cell_type1 = 'F2|F3', cell_type2 = 'tumourcell',
                scdata = sce,
                idents = 'anno', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                keep_significant_only = T,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
                genes=genes
) + scale_edge_alpha(limits = c(0, 7))
p

############### MoDc3 and tumour cell
genes_3 <- c("ACP5","ALOX5AP","ATOX1","BCL2A1","BID","CCL22","CD40","CD58","CD74","CKB","CREB5","CSRP2",    
             "CST3","CST7","CTSZ","DNASE1L3","ECE1","EMP3", "GK5","GPR160","HLA-DPA1","HMG20B","IGSF8","ITM2C",  
             "LPXN","MAP4K4","MCOLN2","MFSD12","MIR3945HG","NAAA" ,    
             "NME4" ,"PDE4B","PLEKHA5","PPFIBP1","PTP4A3" ,"PXDC1",    
             "RAB7B","REV3L","RHOF" ,"S100A11","SDSL","SEMA4A",   
             "SESN3","SH3BP5", "TGFA","TNFRSF18","TRAF4","TSPO" ,    
             "TXNRD1")
data <- read.table('/home/jovyan/farm/CTCL/DEG-condition-psudobulk/DEG-DC-filterByExpr50_100.xls')
genes_4 <- c(rownames(data))
genes_4 <- c(genes_4, genes_3)

p <- plot_cpdb2(cell_type1 = 'MoDC3', cell_type2 = 'tumourcell',
                scdata = sce,
                idents = 'anno', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                desiredInteractions = list(
                  c('MoDC3', 'tumourcell'),
                  c('tumourcell', 'MoDC3')),
                node_group_colors = c(
                  "MoDC3" = "red",
                  "tumourcell" = "dodgerblue3"),
                keep_significant_only = F,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
                genes = genes_4
) + scale_edge_alpha(limits = c(0, 6))
p

p <- plot_cpdb3(cell_type1 = 'MoDC3', cell_type2 = 'tumourcell',
                scdata = sce,
                idents = 'anno', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                keep_significant_only = F,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
                genes = genes_3
) + scale_edge_alpha(limits = c(0, 6))
p
