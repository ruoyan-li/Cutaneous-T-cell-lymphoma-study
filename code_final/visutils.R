loadC2L = function(d,sd2zero=0){
  r = read.csv(paste0(d,'/q05_cell_abundance_w_sf.csv'),row.names = 1,check.names = FALSE)
  sd = read.csv(paste0(d,'/stds_cell_abundance_w_sf.csv'),row.names = 1,check.names = FALSE)
  colnames(r) = sub('q05cell_abundance_w_sf_','',colnames(r))
  colnames(sd) = sub('stdscell_abundance_w_sf_','',colnames(sd))
  sd = sd [rownames(r),colnames(r)]
  r[r<sd*sd2zero] = 0
  mm = do.call(rbind,strsplit(rownames(r),'|',T))
  r$barcode = mm[,2]
  r = split(r,mm[,1])

  for(i in 1:length(r)){
    rownames(r[[i]]) = r[[i]]$barcode
    r[[i]]$barcode = NULL
    r[[i]] = as.matrix(r[[i]])
  }
  r
}
plotNMFConsSummary <- function(
      coefs, cons,
      clcols = NULL,
      max.cex = 4.5/14*8,
      colfun = function(x) {
        visutils::num2col(
          x, c('blue','gray','orange','violet','black'), minx = 0, maxx = 1
        )
      }, ylab.cex = 1, xlab.cex = 1
){
  cons = cons[colnames(coefs),colnames(coefs)]
  cls = apply(coefs,2,which.max)
  if(is.null(clcols))
    clcols = RColorBrewer::brewer.pal(nrow(coefs),'Set3')
  o = names(cls)[order(cls)]
  parl = par(no.readonly=TRUE)
  par(mar=c(3,9,2.5,6),bty='n',oma=c(0,0,1,0),cex=1,tcl=-0.2,mgp=c(1.2,0.3,0),las=1)
  dotPlot(
    t(coefs[, o]),
    rowColours = cbind(clcols[cls[o]]),
    colColours = clcols,
    max.cex = max.cex, ylab.cex = ylab.cex,
    xlab.cex = xlab.cex,
    scaleWM = F,
    colfun = colfun,
    plot.legend = TRUE
  )
  par(parl)
}
plotFeatureProfiles_fun <- function (
  m, features, features_facet = NULL, cols = NULL, sd.mult = 2, legend. = TRUE,
  ylim = NULL, scaleY = TRUE, area.opacity = 0.2, lwd = 2,
  xlab = "Distance (spots)", ylab = "Relative abundance", main = "",
  xlim = NULL, ...
) {
  if (is.null(features_facet))
    features_facet = sapply(features, function(x) dimnames(m)[[3]], simplify = FALSE)
  if (is.null(cols))
      cols = char2col(features)
  if (is.null(names(cols)))
      names(cols) = features
  x = as.numeric(dimnames(m)[[1]])
  areas = lapply(features, function(ct) {
      do.call(rbind, apply(
        m[, gsub("Healthy ", "", ct), unlist(features_facet[ct])], 1,
        function(x) {
          x = x[!is.na(x)]
          data.frame(mean = mean(x), sd = sd(x)/sqrt(length(x)),
              n = length(x))
      }))
  })
  names(areas) = features
  if (scaleY)
      areas = lapply(areas, function(a) {
          a[, 1:2] = a[, 1:2]/max(a$mean, na.rm = TRUE)
          a
      })
  if (is.null(ylim))
      ylim = range(sapply(areas, function(a) c(min(a$mean -
          a$sd * sd.mult, na.rm = TRUE), max(a$mean + a$sd *
          sd.mult, na.rm = TRUE))))
  if (is.null(xlim))
      xlim = range(x)
  plot(1, t = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
      main = main, ...)
  for (n in names(areas)) plotArea(x, areas[[n]][, 1:2], sd.mult = sd.mult,
      col = cols[n], new = FALSE, lwd = lwd, area.transp = area.opacity)
  if (!is.null(legend.) && (!is.logical(legend.) || legend.)) {
      if (!is.list(legend.))
          legend. = list()
      if (is.null(legend.$x)) {
          legend.$x = grconvertX(1, "npc", "user")
          legend.$y = grconvertY(1, "npc", "user")
      }
      legend.$fill = cols
      legend.$legend = names(cols)
      legend.$bty = par("bty")
      legend.$border = NA
      legend.$xpd = NA
      do.call(legend, legend.)
  }
}