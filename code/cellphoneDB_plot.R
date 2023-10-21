library(reticulate)
library(ktplots)
library(SingleCellExperiment)
library(anndata)

setwd('/Users/ruoyanli/Desktop/ALL/CTCL')

#h5ad_file <- "CTCL_sub_nomore500_anno.h5ad"
seurat_object <- readRDS("CTCL_sub_nomore_100_per_donor_celltype_anno.rds")
meta <- read.csv('CTCL_sub_nomore_100_per_donor_celltype_anno_meta.csv')
seurat_object@meta.data <- as.data.frame(meta)
head(seurat_object@meta.data)

#### cellphoneDB output
pvals <- read.delim("pvalues_cpdb4.txt", check.names = FALSE)
means <- read.delim("means_cpdb4.txt", check.names = FALSE)
deconv <- read.delim("deconvoluted_cpdb4.txt", check.names = FALSE)

#adata <- read_h5ad(h5ad_file)
counts <- seurat_object[['RNA']]@counts
#row.names(counts) <- row.names(adata$var)
#colnames(counts) <- row.names(adata$obs)
sce <- SingleCellExperiment(list(counts = counts), colData = seurat_object@meta.data, 
                            rowData = rownames(counts))

######## B/plasma and tumour cells
plot_cpdb(cell_type1 = 'B/plasma', cell_type2 = 'tumourcell', scdata = seurat_object,
          idents = 'anno', means = means, pvals = pvals, 
          keep_significant_only=T) + small_guide()  + small_legend(keysize=1)

p <- plot_cpdb2(cell_type1 = 'B/plasma', cell_type2 = 'tumourcell',
                scdata = sce,
                celltype_key = 'annotation', # column name where the cell ids are located in the metadata
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
)
p

p <- plot_cpdb3(cell_type1 = 'tumourcell', cell_type2 = 'B/plasma',
                scdata = sce,
                celltype_key = 'annotation', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = decon, # new options from here on specific to plot_cpdb3
                keep_significant_only = TRUE,
                standard_scale = TRUE,
                remove_self = TRUE,
)
p

############### F2/F3 and tumour cell
data <- read.table('DEG-Fibro-filterby50_100.xls')
genes <- c(rownames(data))

p <- plot_cpdb3(cell_type1 = 'F2|F3', cell_type2 = 'tumourcell',
                scdata = sce,
                celltype_key = 'annotation', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                keep_significant_only = T,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
                genes = genes
)
p

############### MoDc3 and tumour cell
genes_3 <- c("ACP5","ALOX5AP","ATOX1","BCL2A1","BID","CCL22","CD40","CD58","CD74","CKB","CREB5","CSRP2",    
             "CST3","CST7","CTSZ","DNASE1L3","ECE1","EMP3", "GK5","GPR160","HLA-DPA1","HMG20B","IGSF8","ITM2C",  
             "LPXN","MAP4K4","MCOLN2","MFSD12","MIR3945HG","NAAA" ,    
             "NME4" ,"PDE4B","PLEKHA5","PPFIBP1","PTP4A3" ,"PXDC1",    
             "RAB7B","REV3L","RHOF" ,"S100A11","SDSL","SEMA4A",   
             "SESN3","SH3BP5", "TGFA","TNFRSF18","TRAF4","TSPO" ,    
             "TXNRD1")
data <- read.table('DEG-DC-filterByExpr50_100.xls')
genes_4 <- c(rownames(data))
genes_4 <- c(genes_4, genes_3)

p <- plot_cpdb3(cell_type1 = 'MoDC3', cell_type2 = 'tumourcell',
                scdata = sce,
                celltype_key = 'annotation', # column name where the cell ids are located in the metadata
                means = means,
                pvals = pvals,
                deconvoluted = deconv,
                keep_significant_only = F,
                standard_scale = TRUE,
                remove_self = TRUE,
                col_option = viridis::magma(50),
                genes = genes_3
)
p
