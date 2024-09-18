#!/usr/bin/env R

# ------------------------------------------------------------------------------
# title: Spatial analysis of CTCL data.
# purpose: This script creates visium figures for the CTCL project.
#          It takes the cell2location output.
# created: 2024-07-16 Tue 14:11:45 BST
# updated: 2024-08-02 Fri 13:23:37 BST
# version: 0.0.9
# status: Prototype
#
# maintainer: Ciro Ramírez-Suástegui
# author:
#   - name: Ciro Ramírez-Suástegui
#     affiliation: The Wellcome Sanger Institute
#     email: cs59@sanger.ac.uk, cramsuig@gmail.com
# contributor:
#   - name: Pasha Mazin
#     affiliation: The Wellcome Sanger Institute
#     email: pm19@sanger.ac.uk
# ------------------------------------------------------------------------------

## Environment setup ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dependencies: functions and packages
# basic ----------------------------------------------------
library(plyr)
if (requireNamespace("crayon", quietly = TRUE)) {
  red <- crayon::red
  yellow <- crayon::yellow
  cyan <- crayon::cyan
}else{
  red <- yellow <- cyan = c
}
logging::basicConfig()
section = function (i, tail_n = 60, color = cyan) {
  tail_n <- max(c(tail_n, nchar(i) + 1))
  y <- paste("##", color(i), "##", base::strrep("%", tail_n-nchar(i)), "\n")
  logging::loginfo(y)
}
# tools ----------------------------------------------------
library(Seurat)
library(NMF)
library(visutils)
# in-house/developing --------------------------------------
source(paste0(here::here(), "/code/visutils.R"))

## Global variables and paths ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir <- "/nfs/cellgeni/pasham/projects/2211.adult.skin/data.nfs"
figures_dir <- paste0(here::here(), "/figures/ctcl_visium")
# CTCL as reference, original results
c2l_file <- paste0(data_dir, "/visium/ctcl/c2l.v2.rds")
sufix <- "ref-ctcl_original"
# CTCL as reference
c2l_file <- paste0(here::here(), "/results/ref_ctcl-viss_ctcl_disease.20/predmodel")
sufix <- "ref-ctcl"
# CTCL+Healthy as reference
c2l_file <- paste0(
  here::here(), c("/results/ref_ctcl_h-viss_ctcl_all.20/predmodel",
  "/results/ref_ctcl_h-viss_bayanne_healthy.20/predmodel"))
sufix <- "ref-ctcl-h"

{ section("Loading data") ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  c2l <- unlist(lapply(c2l_file, function(x){
    if (grepl("rds$", x)){
      y <- readRDS(x)
      split.data.frame(y, splitSub(rownames(y), "|", 1))
    }else loadC2L(x)
  }), recursive = FALSE)
  logging::loginfo("Filtering to just the tissue and enhancing images")
  visium <- lapply(
    X = readRDS(paste0(data_dir, "/visium/ctcl/vs.v2.rds")),
    function(v){
      y <- v[, v$is.tissue == 1]
      y@images$slice1@image <- enhanceImage( # More contrasted images
        y@images$slice1@image, wb = TRUE, qs = c(0.1, 0.9)
      )
      return(y)
  })
  mdata <- readRDS(paste0(data_dir, "/visium/ctcl/meta.v2.rds"))
}

## Pre-processing ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{
  samples_subset <- intersect(
    grep("CTCL|WS_D", names(visium), value = TRUE),
    grep("CTCL|WS_D", names(visium), value = TRUE)
  )
  # Setting up cell type colours
  ct_color <- c(
    "Tumour cell" = "#E41A1C",
    "F3" = "#377EB8",
    "F2" = "#4DAF4A",
    "B cell" = "#FF7F00",
    "MoDC3" = "#F781BF",
    "VE3" = "#984EA3",
    "DC2" = "#be5ff5"
  )
  for (i in grep("Healthy", names(ct_color), value = TRUE, invert = TRUE)) {
    ct_color[paste("Healthy", i)] <- grDevices::adjustcolor(
      ct_color[i], alpha.f = 0.3
    )
  }

  logging::loginfo("Matching cell names")
  # tumourcell > tumor_cell
  # MoDC3 > moDC_3
  # B/plasma > B_cell
  ct_names = setNames(nm = unique(unlist(lapply(c2l, colnames))))
  ct_names["tumor_cell"] <- "Tumour cell"
  ct_names["moDC_3"] <- "MoDC3"
  ct_names["B_cell"] <- "B cell"
  for (sample_i in samples_subset) {
    cellnames <- colnames(visium[[sample_i]])
    rownames(c2l[[sample_i]]) <- gsub(".*\\|", "", rownames(c2l[[sample_i]]))
    cellint <- intersect(cellnames, rownames(c2l[[sample_i]]))
    logging::loginfo(paste0(
      sample_i, ": ", length(cellint), "/(c2l:",
      nrow(c2l[[sample_i]]), ",vis:", length(cellnames), ")"
    ))
    if (length(cellnames) > 0) {
      c2l[[sample_i]] <- c2l[[sample_i]][cellnames, ]
    }
    colnames(c2l[[sample_i]]) <- ct_names[colnames(c2l[[sample_i]])]
  }
  # trend line and NMF need the matrix with everything together
  c2lm <- as.data.frame(do.call(rbind, lapply(X = samples_subset, function(x){
    y <- c2l[[x]]; rownames(y) <- paste0(x, ".", rownames(y)); y
  })))
  # visium[["HCA_sCTCL13876503"]] <- rotateVisium(
  #   visium[["HCA_sCTCL13876503"]], n = 1, mirror = TRUE
  # )
}

## Main ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
{ section("Microenvironments (NMF analysis)") ## %%%%%%%%%%%
  c2lmt = as.matrix(c2lm)
  # we will use per-spot normalised celltype abundancies
  c2lmt = sweep(c2lmt, 1, rowSums(c2lmt), '/')
  # we will run nmf N times to assess results stability
  N = 50 # consider to increase it
  # number of factors is almost arbitrary and to be set manually depending
  # on desired granularity of microenvironments.
  for (rank in c(4, 5, 6)) {
    p_name = paste0(file.path(
      figures_dir, paste(
        "mfigure1f",
        paste(c("nmf", rank, "microenvironments"), collapse="-"),
        sufix, sep="_"
      )
    ), ".pdf")
    logging::loginfo(paste(basename(p_name)))
    start <- Sys.time()
    set.seed(1234)
    doMC::registerDoMC(2)
    nmf = plyr::llply(
      .data = 1:N,
      .fun = function(i){
        nmf(c2lmt, rank = rank)
      }, .parallel = TRUE
    )
    logging::loginfo(paste("Elapsed:", format(difftime(Sys.time(), start, unit="min"))))
    best.nmf = nmfGetBest(nmf, getNMFNormFs('max'))
    # plotting
    pdf(gsub(".pdf", "_summary.pdf", p_name), w = 7-2, h = 7+2)
    plotNMFConsSummary(best.nmf$coefn, best.nmf$cons/N, ylab.cex = 0.7, max.cex = 2)
    dev.off()
    pdf(p_name, w = 7+5, h = 7+3)
    plotNMFCons(best.nmf$coefn, best.nmf$cons/N, ylab.cex = 0.7, max.cex = 2)
    dev.off()
  }
}

# Configuring data going into each figure
fconfigs <- list(
  list(
    name = "mfigure3d",
    celltypes = c("Tumour cell", "VE3"), # F3
    samples = c("5" = "HCA_sCTCL13876505", "4" = "HCA_sCTCL13876504")
  ),
  list(
    name = "mfigure3j",
    celltypes = c("Tumour cell", "DC2"), # MoDC3
    samples = c("7" = "HCA_sCTCL13787193", "4" = "HCA_sCTCL13876504")
  ),
  list(
    name = "mfigure3j",
    celltypes = c("Tumour cell", "DC2", "MoDC3"),
    samples = c("7" = "HCA_sCTCL13787193", "4" = "HCA_sCTCL13876504")
  ),
  list(
    name = "mfigure4d",
    celltypes = c("Tumour cell", "B cell"),
    samples = c("5" = "HCA_sCTCL13876505", "8" = "HCA_sCTCL13787192")
  ),
  list(
    name = "efigure4d", # efigure4f
    celltypes = c(rep(list(c("Tumour cell", "F2")), 3), list(c("Tumour cell", "VE3"))),
    samples = c(
      "2" = "HCA_sCTCL13787190", "3" = "HCA_sCTCL13876502",
      "1" = "HCA_sCTCL13787191", "2" = "HCA_sCTCL13787190"
    )
  ),
  list(
    name = "efigure5d", # efigure5f
    celltypes = c("Tumour cell", "DC2", "MoDC3"), # +DC2
    samples = c("5" = "HCA_sCTCL13876505", "6" = "HCA_sCTCL13876503")
  ),
  list(
    name = "efigure6i",
    celltypes = c("Tumour cell", "B cell"),
    samples = c("2" = "HCA_sCTCL13787190", "3" = "HCA_sCTCL13876502")
  ),
  list(
    name = "xfigureNx",
    celltypes = c("Tumour cell", "F2", "VE3"),
    samples = c("7" = "HCA_sCTCL13787193", "4" = "HCA_sCTCL13876504")
  ),
  list(
    name = "xfigureNx", # efigure5f
    celltypes = c("Tumour cell", "F2", "VE3"),
    samples = c("5" = "HCA_sCTCL13876505", "6" = "HCA_sCTCL13876503")
  ),
  list(
    name = "efigure6i",
    celltypes = c("Tumour cell", "F2", "VE3", "B cell", "MoDC3"),
    samples = c("5" = "HCA_sCTCL13876505", "6" = "HCA_sCTCL13876503")
  ),
  list(
    name = "xfigureNx",
    celltypes = c("Tumour cell", "F2", "VE3"),
    samples = c(
      "2" = "HCA_sCTCL13787190", "3" = "HCA_sCTCL13876502",
      "1" = "HCA_sCTCL13787191", "2" = "HCA_sCTCL13787190"
    )
  )
)

{ section("Abudance plot on tissue image") ## %%%%%%%%%%%%%%
  for (i in seq_along(fconfigs)) {
    fconfigs[[i]]$he.grayscale = FALSE
    fconfigs[[i]]$img.alpha = 0.2
    fconfigs[[i]]$he.img.width = 300
  }

  for (fconfig in fconfigs) {
    p_name = paste0(file.path(
      figures_dir, paste(
        fconfig$name, "spatial-abundance",
        paste(make.names(unique(unlist(fconfig$celltypes))), collapse="-"),
        paste(names(fconfig$samples), collapse="-"),
        sufix, sep="_"
      )
    ), ".pdf")
    nrows <- ceiling(sqrt(length(fconfig$samples)))
    ncols <- ceiling(length(fconfig$samples) / nrows)
    # create folder
    logging::loginfo(paste(basename(p_name), nrows, ncols))
    pdf(p_name, w = ncols * 5, h = nrows * 3.5)
    par(mfcol = c(nrows, ncols), mar = c(0, 0, 1, 10), bty = "n")
    if (!is.list(fconfig$celltypes)) {
      fconfig$celltypes <- list(fconfig$celltypes)
    }
    if (length(fconfig$celltypes) == 1) {
      fconfig$celltypes <- rep(fconfig$celltypes, length(fconfig$samples))
    }
    for (index_i in seq_along(fconfig$samples)) {
      sample_i <- fconfig$samples[index_i]
      celltypes <- fconfig$celltypes[[index_i]]
      ftitle <- paste0(
        mdata[sample_i, "stage"], ": ",
        mdata[sample_i, "donor_id"]
      )
      ct_proportions <- c2l[[sample_i]][colnames(visium[[sample_i]]), celltypes]
      p <- plotVisiumMultyColours(
        visium[[sample_i]], ct_proportions,
        zfun = function(x) x^2, scale.per.colour = TRUE,
        col = ct_color[colnames(ct_proportions)], mode = "mean",
        he.grayscale = fconfig$he.grayscale, img.alpha = fconfig$img.alpha,
        main = ftitle, legend.ncol = 2,
        min.opacity=0
      )
    }
    dev.off()
  }
}

{ section("Abundance along dermis-epidermis axis") ## %%%%%%
  # define distance ----------------------------------------
  for (sample_i in samples_subset) {
    d <- as.matrix(dist(
      visium[[sample_i]]@images$slice1@coordinates[, c("imagerow", "imagecol")]
    ))
    spot_dist <- min(d[upper.tri(d)])
    visium[[sample_i]]$dist2junction <- -abs(
      apply(d[, visium[[sample_i]]$is.surface], 1, min) / spot_dist
    )
  }
  columns_l <- lapply(visium, function(x) colnames(x@meta.data) )
  columns_u <- unique(unlist(columns_l))
  columns_s <- sapply(columns_l, function(x) columns_u %in% x ) # shared columns
  columns_u <- columns_u[rowSums(columns_s) == length(columns_l)]
  spots <- as.data.frame(data.table::rbindlist(
    lapply(samples_subset, function(x) visium[[x]]@meta.data[, columns_u] ))
  )
  rownames(spots) <- paste0(spots$sample_id, ".", spots$barcode)
  # binarise the distance ----------------------------------
  spots$dist2junction = round(as.numeric(spots$dist2surf))
  # cor(spots$dist2junction, spots$dist2surf) # 0.9988894
  # mean(spots$dist2junction == spots$dist2surf) # 0.06986634
  # Calculating matrix with average ------------------------
  # celltype abundances for each distance bin and each sample
  dfsmtx.c2l = makeDistFeatureSampleTable(
    dist = spots$dist2junction,
    sample = spots$sample_id,
    data = c2lm,
    per.spot.norm = TRUE
  )
  # dfsmtx.c2l["0", fconfigs_trends[[1]]$celltypes[1], tail(unique(spots$sample_id), 1)]
  # mean(c2lm_condition[, fconfigs_trends[[1]]$celltypes[1]])

  
  fselect <- c("mfigure3d", "mfigure3j", "efigure5d", "mfigure4d", "xfigureNx")
  fnames <- sapply(fconfigs, "[[", "name" )
  temp <- sapply(fselect, function(x) min(which(fnames %in% x)) )
  fconfigs_trends <- fconfigs[temp] # [c(1, 2, 4, 5)]
  fconfigs_trends[[1]]$name = "mfigure3e"
  fconfigs_trends[[2]]$name = "mfigure3k"
  ftitle <- "All CTCL combined"
  if (any(grepl("WS_D", samples_subset))) {
    ftitle <- "CTCL vs Healthy"
    for (i in seq_along(fconfigs_trends)) {
      fconfigs_trends[[i]]$name <- paste0(fconfigs_trends[[i]]$name, "-healthy")
      fconfigs_trends[[i]]$celltypes <- c(
        fconfigs_trends[[i]]$celltypes,
        paste("Healthy", grep("umour", fconfigs_trends[[i]]$celltypes, value = TRUE, invert = TRUE))
      )
      fconfigs_trends[[i]]$features_facet <- sapply(
        X = fconfigs_trends[[i]]$celltypes,
        FUN = function(x) {
        if (any(grepl("Healthy", x))){
          grep("WS_D", samples_subset, value = TRUE)
        }else{
          grep("CTCL", samples_subset, value = TRUE)
        }
      }, simplify = FALSE)
    }
  }

  for (fconfig in fconfigs_trends) {
    p_name = paste0(file.path(
      figures_dir, paste(
        fconfig$name, "trend-abundance",
        paste(make.names(unique(unlist(fconfig$celltypes))), collapse="-"),
        sufix, sep="_"
      )
    ), ".pdf")
    logging::loginfo(paste(basename(p_name)))
    if (is.list(fconfig$celltypes)) next
    pdf(p_name, w = 7-2, h = 7-3)
    par(mar = c(4, 4, 2, 8), bty = "n")
    plotFeatureProfiles_fun(
      dfsmtx.c2l,
      features = fconfig$celltypes,
      features_facet = fconfig$features_facet,
      cols = ct_color[fconfig$celltypes],
      lwd = 5,
      sd.mult = 1,
      ylim = c(0, 1.3),
      xlim = range(-15:0),
      main = ftitle,
      xlab = "Distance to surface interface (spots)",
    )
    dev.off()
  }
}

## Conclusions ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Save ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
