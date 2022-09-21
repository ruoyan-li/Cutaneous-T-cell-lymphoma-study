library(infercnv)
library(Seurat)
library(Matrix)

CTCL.integrated <- readRDS("/lustre/scratch117/cellgen/team205/rl20/CTCL/object/CTCL.integrated.rmdoublets.rds")
Donor_1 <- subset(x=CTCL.integrated, cells=
		  rownames(CTCL.integrated@meta.data[which(CTCL.integrated@meta.data$donor==1), ]))
Donor_1_for_CNV <- subset(x=Donor_1, cells = 
                           rownames(Donor_1@meta.data[which(Donor_1@meta.data$broad_type=="T-cell" |
                                                                    Donor_1@meta.data$broad_type=="Pericyte" |
                                                                    Donor_1@meta.data$broad_type=="V-Endo" | 
                                                                    Donor_1@meta.data$broad_type=="Melanocyte" |  
                                                                    Donor_1@meta.data$broad_type=="Macro"), ]))

gene_list <- read.table("/home/jovyan/farm/CTCL/CNV_donor1/gencode_v21_gen_pos.complete.match.mat.order.intersect.txt")
counts_matrix <- Donor_1_for_CNV@assays$RNA@counts[as.vector(gene_list[,1]), ]


infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="/home/jovyan/farm/CTCL/CNV_donor1/dataforCNV_meta.2.txt",
                                    delim="\t",
                                    gene_order_file="/home/jovyan/farm/CTCL/CNV_donor1/gencode_v21_gen_pos.complete.match.mat.order.intersect.txt",
                                    ref_group_names=c("Pericyte","V-Endo","Melanocyte","Macro","CD8-Tcell")) 

infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/lustre/scratch117/cellgen/team205/rl20/CTCL/CNV/infercnv_out_donor1", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             num_threads=8,
			     HMM=TRUE)
