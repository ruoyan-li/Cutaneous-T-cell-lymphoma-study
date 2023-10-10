#### Author ---  Pasha Mazin --- Wellcome Sanger Institute

library(Seurat)
library(plyr)
library(NMF)
library(visutils)
source('/nfs/cellgeni/pasham/rcode/visutils/R/tissue.depth.R')
source('../2302.fetal.skin/src/my/visium/c2l.utils.R')

# new batch
#Spaceranger outputs are here: /lustre/scratch126/cellgen/team205/rl20/CTCL/Visium
#Cell2loc outputs are here: /lustre/scratch126/cellgen/team205/rl20/CTCL/Visium/cell2location_map_allsamples
gene.descr=readRDS('data.nfs/gene.descr.rds')
vs = readRDS('data.nfs/visium/ctcl/vs.v2.rds')
c2l = readRDS('data.nfs/visium/ctcl/c2l.v2.rds')
c2l0 = readRDS('data.nfs/visium/ctcl/c2l0.v2.rds')
mctcl = readRDS('data.nfs/visium/ctcl/meta.v2.rds')
vsf = lapply(vs,function(v)v[,v$is.tissue==1])

splitSub = function(x,del,inx,fixed=TRUE){
  sapply(strsplit(x,del,fixed=fixed),'[',inx)
}


# c2l ---------------------
# reticulate::use_condaenv('scvi')
# # healthy
# c2lh = read.csv('data.lustre/visium/ctcl/c2l_healthy/pred/ref1.20/predmodel/q05_cell_abundance_w_sf.csv',row.names = 1,check.names = FALSE)
# # ctcl
# ann = anndata::read_h5ad('/lustre/scratch126/cellgen/team205/rl20/CTCL/Visium/cell2location_map_allsamples/cell2location_map/sp.h5ad')
# c2l =  ann$obsm$q05_cell_abundance_w_sf
# rownames(c2l) = paste0(substr(rownames(c2l),28,44),'|',substr(rownames(c2l),nchar(rownames(c2l))-17,nchar(rownames(c2l))))
# table(colnames(c2l) == colnames(c2lh))
# c2l = rbind(c2l,c2lh)
# colnames(c2l) = gsub('q05cell_abundance_w_sf_','',colnames(c2l))
# 
# c2l = as.matrix(c2l)
# c2l[1:2,]
# saveRDS(c2l,'data.nfs/visium/ctcl/c2l.v2.rds')
# 
# # set <sd to 0
# sdh = read.csv('data.lustre/visium/ctcl/c2l_healthy/pred/ref1.20/predmodel/stds_cell_abundance_w_sf.csv',row.names = 1,check.names = FALSE)
# sd =  ann$obsm$stds_cell_abundance_w_sf
# rownames(sd) = paste0(substr(rownames(sd),28,44),'|',substr(rownames(sd),nchar(rownames(sd))-17,nchar(rownames(sd))))
# table(colnames(sd) == colnames(sdh))
# sd = rbind(sd,sdh)
# colnames(sd) = gsub('stdscell_abundance_w_sf_','',colnames(sd))
# sd = as.matrix(sd)
# c2l0 = c2l
# c2l0[sd>=c2l] = 0
# table(c2l0==0)
# saveRDS(c2l0,'data.nfs/visium/ctcl/c2l0.v2.rds')
# 
# 
# # prepare metadata
# mctcl = read.csv('/lustre/scratch126/cellgen/team205/rl20/CTCL/Visium/meta_visium.csv')
# mctcl$sanger_id = substr(mctcl$sample_id,28,44)
# mctcl$dataset = 'ctcl'
# mctcl = mctcl[,-1]
# mh = readRDS('data.nfs/visium/visium.sample.info.rds')
# mh = mh[mh$State == 'H' & !grepl('face|scalp',mh$Anatomical.site,ignore.case = TRUE) & mh$qc == 'ok',]
# mh = mh[,c('Donor','Age','Sex','State','Anatomical.site','File.ID')]
# mh$dataset = 'bayanne'
# colnames(mh) = colnames(mctcl)
# mctcl = rbind(mctcl,mh)
# rownames(mctcl) = mctcl$sanger_id
# 
# sids = unique(sapply(strsplit(rownames(c2l),'|',fixed=TRUE),'[',1))
# length(sids)
# sum(mctcl$sanger_id %in% sids)
# saveRDS(mctcl,'data.nfs/visium/ctcl/meta.v2.rds')
# 
# # load visiums ######
# #rsync  -urlptDR --progress --exclude '.bam*' /lustre/scratch126/cellgen/team205/rl20/CTCL/Visium/./spaceranger* /lustre/scratch127/cellgen/cellgeni/tmp/ctcl/
# sids = list.dirs('data.lustre/visium/ctcl/',recursive = F)
# vs = lapply(mctcl$sanger_id[mctcl$dataset=='ctcl'],function(s){
#   dir = list.dirs('data.lustre/visium/ctcl/spaceranger/',rec=FALSE)
#   dir = dir[grep(s,dir)]
#   myLoad10X_Spatial(dir,filter.matrix = F,ens_id = T)
#   })
# names(vs) = mctcl$sanger_id[mctcl$dataset=='ctcl']
# # add healthy
# mh = readRDS('data.nfs/visium/visium.sample.info.rds')
# vsh = readRDS('data.nfs/visium/vis.f.merged.spots.rds')
# names(vsh) = mh[names(vsh),'File.ID']
# vsh = vsh[mctcl$sanger_id[mctcl$dataset=='bayanne']]
# vs = c(vs,vsh)[mctcl$sanger_id]
# sapply(vs,dim)
# saveRDS(vs,'data.nfs/visium/ctcl/vs.v2.rds')
# visium --------------
# _qc -------------------

f = function(x){
  r = x
  r[x<100] = '<100'
  r[x>=100 & x < 500] = '100-499'
  r[x>=500 & x < 999] = '500-999'
  r[x>=1000 & x < 5000] = '1000-5000'
  r[x>5000] = '>5000'
  r
}

cf = c('<100'="#E41A1C",'100-499'="#FF7F00",'500-999'="#FFFF33",'1000-5000'="#377EB8",'>5000'="#4DAF4A")


pdf('figures/visium/ctcl/v2/qc.pdf',w=8*3.5,h=6*3)
par(mfcol=c(6,8),bty='n',mar=c(0,0,1,6))

for(n in mctcl$sanger_id[mctcl$dataset=='ctcl']){
  plotVisium(vs[[n]],vs[[n]]$nCount_Spatial,main=n,cex=0)
  plotVisium(vs[[n]],vs[[n]]$nCount_Spatial,main=n,cex=scaleTo(log1p(vs[[n]]$nCount_Spatial)))
  plotVisium(vs[[n]],vs[[n]]$nCount_Spatial,zfun=log1p,main=n,cex=scaleTo(log1p(vs[[n]]$nCount_Spatial)))
  plotVisium(vs[[n]],as.character(vs[[n]]$is.tissue),main=n,cex=0.6)
  plotVisium(vs[[n]],main=n,cex=0.6,spot.filter = vs[[n]]$is.tissue==1)
  plotVisium(vs[[n]],f(vs[[n]]$nCount_Spatial),z2col=cf,main=n,cex=1,spot.filter = vs[[n]]$is.tissue==1)
}
dev.off()

# plot percelltype ########
he.grayscale  = TRUE
img.alpha = 0.4
he.img.width = 300

c2ls = list(c2l     = split.data.frame(c2l,splitSub(rownames(c2l),'|',1)),
            c2l.sd0 = split.data.frame(c2l0,splitSub(rownames(c2l0),'|',1)))

first = c('tumourcell','B/plasma','VE2','F2','F3','MoDC3','Th','pDC','Proliferating_APC','Proliferating_T-cell','Pericyte_1','Differentiated_KC*','LC1','LC2','Mac2','Mast')
nrow=3
ncol=8


for(sid in names(vsf)){
  br=paste0(sid,'|',colnames(vsf[[sid]]))
  c2ls$c2l[[sid]] = c2ls$c2l[[sid]][br,]
  c2ls$c2l.sd0[[sid]] = c2ls$c2l.sd0[[sid]][br,]  
}

o = order(c('H'=4,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = F)
table(mctcl$stage[o])
for(a in names(c2ls)){
  print(a)
  celltypes = names(sort(-apply(do.call(rbind,c2ls[[a]]),2,sum)))
  celltypes = c(first,setdiff(celltypes,first))
  pdf(paste0('figures/visium/ctcl/v2/c2l.by.celltype.',a,'.pdf'),w=ncol*3.5,h=nrow*3)
  for(ct in celltypes){
    par(mfrow=c(nrow,ncol),mar=c(0.1,0.1,1.2,4),bty='n',oma=c(0,0,1,0))
    print(ct)
    for(sid in mctcl$sanger_id[o]){
      v = c2ls[[a]][[sid]][,ct]
      cex = scaleTo(log1p(v))
      cex[is.na(cex)] = 0
      plotVisium(vsf[[sid]],v,zfun = log1p,cex = cex,
                 main=paste0(mctcl[sid,'stage'],': ', sid),he.img.width=he.img.width,he.grayscale=he.grayscale,img.alpha=img.alpha)
    }
    mtext(ct,3,outer = TRUE)
  }
  dev.off()
}


# only ctcl #######
o = order(c('H'=4,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = F)
o = o[mctcl$stage[o]!='H']
table(mctcl$stage[o])
nrow=5
ncol=8

for(a in names(c2ls)){
  print(a)
  pdf(paste0('figures/visium/ctcl/v2/c2l.by.celltype.',a,'.only.ctcl.pdf'),w=ncol*3.5,h=nrow*3)
  par(mfrow=c(nrow,ncol),mar=c(0.1,0.1,2.2,4),bty='n',oma=c(0,0,0,0))
  for(ct in first){
    print(ct)
    for(sid in mctcl$sanger_id[o]){
      v = c2ls[[a]][[sid]][,ct]
      cex = scaleTo(log1p(v))
      cex[is.na(cex)] = 0
      plotVisium(vsf[[sid]],v,zfun = log1p,cex = cex,
                 main=paste0(mctcl[sid,'stage'],': ', sid,'\n',ct),he.img.width=he.img.width,he.grayscale=he.grayscale,img.alpha=img.alpha)
    }
  }
  dev.off()
}


cts = c('tumourcell','B/plasma','VE2','F3','MoDC3','Differentiated_KC*')#,'VE2','MoDC3','Th','pDC','Proliferating_APC','Proliferating_T-cell','Pericyte_1',,'LC1','LC2','Mac2','Mast')
a='c2l'
pdf(paste0('figures/visium/ctcl/v2/c2l.multiple_ct.',a,'.only.ctcl.v1.pdf'),w=8*3.5,h=3*3)
o = order(c('H'=4,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = F)
cols= setNames(RColorBrewer::brewer.pal(9,'Set1')[1:length(cts)],cts)
par(mfcol=c(3,8),mar=c(1,1,1,6),bty='n')
for(sid in mctcl$sanger_id[o]){
  v = c2ls[[a]][[sid]][,cts]
  cex = scaleTo(log(rowSums(v)))
  plotVisiumMultyColours(vsf[[sid]],v,zfun = function(x)x^2,scale.per.colour = TRUE,col=cols,mode = 'mean',he.grayscale=he.grayscale,img.alpha=img.alpha,main=paste0(mctcl[sid,'stage'],': ',sid))
  plotVisium(vsf[[sid]],pie.fracs=v,pie.cols=cols,he.grayscale=he.grayscale,img.alpha=img.alpha,cex=cex,z2col=cols)
  v = sweep(v,2,apply(v,2,max),'/')
  plotVisium(vsf[[sid]],pie.fracs=v,pie.cols=cols,he.grayscale=he.grayscale,img.alpha=img.alpha,cex=cex,z2col=cols,main='norm by max')
}
dev.off()

# _define border ---------------------------
# _ctcl
brd = lapply(vs[mctcl$dataset=='ctcl'],findTissueBorder)
sapply(brd,function(x)table(x$rc$tissue.piece))
# retain only largest piece, other are obviously debris
for(n in mctcl$sanger_id[mctcl$dataset=='ctcl']){
  tissue = names(sort(table(brd[[n]]$rc$tissue.piece),decreasing = TRUE))[1]
  tissue = ((!is.na(brd[[n]]$rc$tissue.piece)) & brd[[n]]$rc$tissue.piece == tissue)+0
  vs[[n]]$is.tissue = vs[[n]]@images$slice1@coordinates$tissue = tissue
  vs[[n]]@meta.data$is.border = FALSE
  vs[[n]]@meta.data$is.border[vs[[n]]$is.tissue == 1 & brd[[n]]$rc$is.border] = TRUE
  vs[[n]]@meta.data$border.inx = brd[[n]]$rc$border.inx
  vs[[n]]@meta.data$border.inx[!vs[[n]]@meta.data$is.border] = NA
}

pdf('figures/visium/ctcl/v2/define.surface.pdf',w=7*3.5,h=4*3)
par(mfcol=c(5,8),bty='n',bty='n')
for(n in mctcl$sanger_id[mctcl$dataset=='ctcl']){
  par(mar=c(0,0,1,6))
  c=plotVisium(vs[[n]],vs[[n]]$border.inx,main=n,cex=0.6,spot.filter = vs[[n]]$is.tissue==1)
  c = c[!is.na(c$z),]
  text(c$x,c$y,c$z,cex=0.65)
  o = order(vs[[n]]$border.inx)
  
  for(gn in c('KRT1','KRT10')){
    par(mar=c(0,0,1,6))
    gid = gene.descr$ens_id[gene.descr$name==gn]
    exp = vs[[n]]@assays$Spatial@counts[gid,]/vs[[n]]$nCount_Spatial*1e4
    plotVisium(vs[[n]],exp,main=n,cex=0.6,spot.filter = vs[[n]]$is.tissue==1,legend.args = list(title=paste0('cpm(',gn,')')))
    par(mar=c(4,4,1,1))
    plot(vs[[n]]$border.inx[o],exp[o]+1,log='y',pch=16,t='p',xlab='border position',ylab=paste0('cpm(',gn,')'),main=gn)#,col=ifelse(vs[[n]]$is.surface[o],'red','black'))
  }
}
dev.off()



# defined manually based on figure below 
border.spots = list(32:50,46:64,15:49,c(69:109,1:8),37:60,c(102:123,1:34),31:84,2:58)
names(border.spots) = mctcl$sanger_id[mctcl$dataset=='ctcl']
for(i in mctcl$sanger_id[mctcl$dataset=='ctcl'])
  vs[[i]]@meta.data$is.surface = vs[[i]]@meta.data$is.border & vs[[i]]@meta.data$border.inx %in% border.spots[[i]]

# define surface for healthy as border epi:
for(i in mctcl$sanger_id[mctcl$dataset=='bayanne'])
  vs[[i]]@meta.data$is.surface = !is.na(vs[[i]]@meta.data$border.inx) #& vs[[i]]@meta.data$man.ann=='epi' # border is already set to epi only



# define distance 
for(i in 1:length(vs)){
  d = as.matrix(dist(vs[[i]]@images$slice1@coordinates[,c('imagerow','imagecol')]))
  spot_dist = min(d[upper.tri(d)])
  vs[[i]]$dist2surf = apply(d[,vs[[i]]$is.surface],1,min)/spot_dist
}


for(n in names(vs)){
  vs[[n]]@meta.data$sample_id = n
  vs[[n]]@meta.data$barcode = rownames(vs[[n]]@meta.data)
  vs[[n]]@meta.data$dist2surf = -abs(vs[[n]]@meta.data$dist2surf)
}

for(i in mctcl$sanger_id[mctcl$dataset=='bayanne'])
  vs[[i]] = vs[[i]][,vs[[i]]$dist2hf>=3 & vs[[i]]$man.ann %in% c('epi','dermis') & vs[[i]]$nCount_Spatial >= 500]

#saveRDS(vs,'data.nfs/visium/ctcl/vs.v2.rds')
he.img.width=300
pdf('figures/visium/ctcl/v2/surface.pdf',w=7*3.5,h=7*3)
par(mfcol=c(7,8),bty='n',bty='n')
for(n in mctcl$sanger_id){
  par(mar=c(0,0,1,6))
  if(mctcl[n,'dataset'] == 'bayanne')
    cex = scaleTo(sqrt(vs[[n]]$nspots),1,sqrt(8),minx = 1,maxx = sqrt(7))*0.9
  else
    cex = 1
  c=plotVisium(vs[[n]],vs[[n]]$border.inx,main=n,spot.filter = vs[[n]]$is.tissue==1,he.img.width=he.img.width,cex=0.6)
  c = c[!is.na(c$z),]
  text(c$x,c$y,c$z,cex=0.65)
  o = order(vs[[n]]$border.inx)
  
  plotVisium(vs[[n]],ifelse(vs[[n]]$is.surface,'red','#FFFFFF00'),main=n,spot.filter = vs[[n]]$is.tissue==1,he.img.width=he.img.width,cex=cex)
  plotVisium(vs[[n]],vs[[n]]$dist2surf,main=n,spot.filter = vs[[n]]$is.tissue==1,legend.args =  list(title='dist2surf\nspots'),he.img.width=he.img.width,cex=cex)
  
  for(gn in c('KRT1','KRT10')){
    par(mar=c(0,0,1,6))
    gid = gene.descr$ens_id[gene.descr$name==gn]
    exp = vs[[n]]@assays$Spatial@counts[gid,]/vs[[n]]$nCount_Spatial*1e4
    plotVisium(vs[[n]],exp,main=n,spot.filter = vs[[n]]$is.tissue==1,legend.args = list(title=paste0('cpm(',gn,')')),cex=cex,he.grayscale=TRUE,img.alpha=0.5,he.img.width=he.img.width)
    par(mar=c(4,4,1,1))
    plot(vs[[n]]$border.inx[o],exp[o]+1,log='y',pch=16,t='p',xlab='border position',ylab=paste0('cpm(',gn,')'),main=gn,col=ifelse(vs[[n]]$is.surface[o],'red','black'))
  }
}
dev.off()

# c2l on depth ---------------------------
columns = intersect(colnames(vs$HCA_sCTCL13787190@meta.data),colnames(vs$WS_D_SKNsp12767489@meta.data))
spotmeta = do.call(rbind,lapply(vs,function(v){v@meta.data[,columns]}))
spotmeta = spotmeta[spotmeta$is.tissue==1,]
c2l[1:2,]
rownames(spotmeta) = paste0(spotmeta$sample_id,'|',spotmeta$barcode)
table(rownames(spotmeta) %in% rownames(c2l))

c2l  = c2l [rownames(spotmeta),]
c2l0 = c2l0[rownames(spotmeta),]

sd = FALSE # F =  use c2l as is; T = set <sd to 0
suffix = ifelse(sd,'.sd0','.asis')
if(sd){
  c2ln = sweep(c2l0,1,rowSums(c2l0),'/')
}else{
  c2ln = sweep(c2l,1,rowSums(c2l),'/')
}

hmo = rev(colnames(c2ln))
hmcols = num2col(1:100,c('white','gray','orange','red'))

dfs = makeDistFeatureSampleTable(round(spotmeta$dist2surf),spotmeta$sample_id,c2ln,f = spotmeta$dist2surf >= -15)
dfsm = apply(dfs,1:2,mean,na.rm=T)
dfsm = sweep(dfsm,2,apply(dfsm,2,max,na.rm=T),'/')
#dfsm = dfsm[,order(apply(dfsm,2,which.max))]

pdf(paste0('figures/visium/ctcl/v2/tissue.depth.hm',suffix,'.pdf'),w=20,h=8)
par(mfrow=c(1,5),mar=c(4,10,1,1),bty='n')
imageWithText(dfsm,'',xlab='Distance to surface',main='Mean sample',col=hmcols)
for(n in c('H','1A','1B','2B')){
  x = apply(dfs[,,mctcl$stage==n],1:2,mean,na.rm=T)
  x = sweep(x,2,apply(x,2,max,na.rm=T),'/')
  x = x[,hmo]
  imageWithText(x,'',xlab='Distance to surface',main=n,col=hmcols)
}
dev.off()


c2ll = lapply(names(vs),function(n)c2ln[paste0(vsf[[n]]$sample_id,'|',vsf[[n]]$barcode),])
names(c2ll) = names(vs)
table(mctcl$stage)
ct2depth = makeCTonDepthAnalysis(vsf,c2ll,colnames(c2ln),'dist2surf',mctcl$stage,min.umi = 500,min.spots = 30,
                                 min.dist2hf = 0,df = 5,norm = TRUE)
statcols = char2col(mctcl$stage)


ctorder = c('tumourcell','B/plasma','VE2','F2','F3','MoDC3','Th')
ctorder = c(ctorder,colnames(c2ln)[!(colnames(c2ln) %in% ctorder)])

# new sign
ctorder = c('pDC','Proliferating_APC','Proliferating_T-cell','Pericyte_1','Differentiated_KC*','LC1','LC2','Mac2','Mast') # LC1 is the same as LC2

pdf(paste0('figures/visium/ctcl/v2/tissue.depth.detailed',suffix,'.new.sign.pdf'),w=6*3.5,h=3.5*3)
par(mfcol=c(3,6),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,1.5,0),oma=c(0,0,0,1),bty='n')
for(ct in ctorder){
  par(mar=c(3,3,1.5,0))
  print(ct)
  xlim=c(-12,0)
  
  t = ct2depth[[ct]]
  a2ds = t$a2ds.raw
  a2d = t$a2d.raw
  
  ylim = range(0,unlist(lapply(a2ds,function(x)c(x$mean+x$sd*2)[x$x>xlim[1]])),na.rm=T)
  
  plot(1,t='n',xlim=xlim,ylim=ylim,main=ct,xlab='Distance to surface interface (spots)',ylab='cell type fraction')
  
  for(s in names(statcols))
    plotArea(a2ds[[s]]$x,a2ds[[s]][,c('mean','sd')],new = F,col=statcols[s],lwd=3)
  abline(v=0,lty=3)
  
  z = lapply(a2d,function(zz){zz$x[order(zz$y,decreasing = T)[1:max(1,nrow(zz)/10)]]})[names(statcols)]
  
  pv = do.call(rbind,lapply(names(z),function(n)data.frame(n=n,x=z[[n]])))
  pv = anova(lm(pv$x ~ pv$n))[1,5]
  pv = format(pv,scientific = TRUE,digits = 2)
  
  par(mar=c(3,10,1.5,0))
  z = lapply(z,function(x)x+runif(length(x),-0.4,0.4))
  beeswarm::beeswarm(z,cex=0.8,horizontal=TRUE,col=statcols,xlim=xlim,main=ct,xlab='Distance to surface interface (spots)',las=1,pch=16)
  abline(v=0,lty=3)
  boxplot(z,outline = F,col=statcols,horizontal = TRUE,main=paste0(ct,' (p=',pv,')'),xlab='Distance to surface interface (spots)',las=1)
  abline(v=0,lty=3)
}
dev.off()



ct2depth_sam = makeCTonDepthAnalysis(vsf,c2ll,colnames(c2ln),'dist2surf',mctcl$sanger_id,min.umi = 500,min.spots = 30,
                                     min.dist2hf = 0,df = 5,norm = TRUE)



tumor.me = c('tumourcell','B/plasma','F3','Differentiated_KC*')
ct.cols = char2col(tumor.me)[tumor.me]

o = order(c('H'=0,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = F)
o = o[mctcl$stage[o] != 'H']


pdf('figures/visium/ctcl/v2/tissue.depth.sel_celltypes.pdf',w=4*4,h=8)
par(mfrow=c(2,4),mar=c(4,4,1,1),bty='n')
for(sid in mctcl$sanger_id[o]){
  plot(1,t='n',xlim=c(-15,0),ylim=c(0,1),main=paste0(mctcl[sid,'stage'],': ',sid),xlab='Distance to surface interface (spots)',ylab='relative celltype abundance')
  for(ct in tumor.me){
    a2ds = ct2depth_sam[[ct]]$a2ds.raw
    f = !is.na(a2ds[[sid]]$mean +  2*a2ds[[sid]]$sd)
    mx = max(a2ds[[sid]]$mean +  2*a2ds[[sid]]$sd,na.rm=TRUE)
    plotArea(a2ds[[sid]]$x[f],a2ds[[sid]][f,c('mean','sd')]/mx,new = F,col=ct.cols[ct],lwd=5,area.transp=0.1)
  }
  legend('topleft',lwd=2,col=ct.cols,legend=names(ct.cols),bty='n')
}
dev.off()




mctcl$stage1 = mctcl$stage
mctcl$stage1[mctcl$stage1!='H'] = 'CTCL'
table(mctcl$stage)



pdf(paste0('figures/visium/ctcl/v2/diff.tissue.depth.hm',suffix,'.pdf'),w=11,h=8)
par(mfrow=c(1,3),mar=c(4,7,1,1),cex.axis=0.8,bty='n',oma=c(0,0,0,4))

for(cp in list(c('H','1A'),c('H','1B'),c('H','2B'),c('1A','2B'),c('1A','1B'))){
  plotDF.HM(dfs[,hmo,],mctcl$stage,cond.pair=cp,l2fc.zlim=c(-6,6),log.pseudocount=0.01,col=hmcols)
}
plotDF.HM(dfs,mctcl$stage1,cond.pair=c('H','CTCL'),l2fc.zlim=c(-6,6),log.pseudocount=0.01,col=hmcols)
dev.off()

# pairwise plots ###########
mctcl$stage1 = mctcl$stage
mctcl$stage1[mctcl$stage1!='H'] = 'CTCL'
ct2depth_diag = makeCTonDepthAnalysis(vsf,c2ll,colnames(c2ln),'dist2surf',mctcl$stage1,min.umi = 500,min.spots = 30,
                                      min.dist2hf = 0,df = 5,norm = TRUE)
xrange = as.character(-15:0)

a = 'c2l'
o = order(c('H'=0,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = F)
o = o[mctcl$stage[o] != 'H']

ctorder = c('tumourcell','B/plasma','VE2','F2','F3','MoDC3','Differentiated_KC*')
ct.cols = setNames(RColorBrewer::brewer.pal(8,'Set1')[c(1,5,4,3,2,8,7)],ctorder)
plot(1:length(ct.cols),cex=3,pch=19,col=ct.cols)
#vsf_ = vsf
for(sid in names(vsf)){
  print(sid)
  vsf[[sid]]@images$slice1@image = enhanceImage(vsf_[[sid]]@images$slice1@image,wb = TRUE,qs = c(0.1,0.9))
}
img.alpha=0.2

pdf('figures/visium/ctcl/v2/tissue.depth.sel_celltypes_pairs.pdf',w=4*9,h=8)
par(mfcol=c(2,9),mar=c(3,3,1,5),tcl=-0.2,mgp=c(1.3,0.3,0),bty='n',oma=c(0,0,1.2,1))
for(ct2 in c('B/plasma','VE2','F2','F3','MoDC3','Differentiated_KC*')){
  cts = c('tumourcell',ct2)
  plot(1,t='n',xlim=range(as.numeric(xrange)),ylim=c(0,1),main='All CTCL combined',xlab='Distance to surface interface (spots)',ylab='relative celltype abundance')
  for(ct in cts){
    a2ds = ct2depth_diag[[ct]]$a2ds.raw$CTCL[xrange,]
    f = !is.na(a2ds$mean +  2*a2ds$sd)
    mx = max(a2ds$mean +  2*a2ds$sd,na.rm=TRUE)
    plotArea(a2ds$x[f],a2ds[f,c('mean','sd')]/mx,new = F,col=ct.cols[ct],lwd=5,area.transp=0.2)
  }
  
  plot.new()
  legend('topleft',fill=ct.cols[cts],border=NA,bty='n',legend=cts,title='Cell type')
  for(sid in mctcl$sanger_id[o]){
    plot(1,t='n',xlim=c(-15,0),ylim=c(0,1),main=paste0(mctcl[sid,'stage'],': ',sid),xlab='Distance to surface interface (spots)',ylab='relative celltype abundance')
    for(ct in cts){
      a2ds = ct2depth_sam[[ct]]$a2ds.raw[[sid]][xrange,]
      f = !is.na(a2ds$mean +  2*a2ds$sd)
      mx = max(a2ds$mean +  2*a2ds$sd,na.rm=TRUE)
      plotArea(a2ds$x[f],a2ds[f,c('mean','sd')]/mx,new = F,col=ct.cols[ct],lwd=5,area.transp=0.2)
    }
    v = c2ls[[a]][[sid]][,cts]
    cex = scaleTo(log(rowSums(v)))
    zz = plotVisiumMultyColours(vsf[[sid]],v,zfun = function(x)x^2,scale.per.colour = TRUE,col=ct.cols[cts],mode = 'mean',he.grayscale=he.grayscale,img.alpha=img.alpha,main=paste0(mctcl[sid,'stage'],': ',sid),bg = '#00000000')
    zz = zz[vsf[[sid]]$is.surface,1:2]
    inx = vsf[[sid]]$border.inx[vsf[[sid]]$is.surface]
    #inx = (inx + max(inx)) %% max(inx)
    zz = zz[order(inx),]
    inx = inx[order(inx)]
    gap = which(inx[-1] - inx[-length(inx)] > 1)
    if(length(gap)==1){
      inx[1:gap] = max(inx) + inx[1:gap]
      zz = zz[order(inx),]
    }
    lines(zz,pch=19,col='black')
  }
  mtext(cts[2],side = 3,outer = TRUE)
}
dev.off()


sid='HCA_sCTCL13787193'

# NMF ctcl ###############
all(rownames(c2l)==rownames(c2l0))
f = splitSub(rownames(c2l),'|',1) %in% mctcl$sanger_id[mctcl$dataset=='ctcl']
table(f)


c2l_ = c2l[f,]
c2l0_ = c2l0[f & apply(c2l0,1,sum)>0,]

c2l_n = sweep(c2l_,1,apply(c2l_,1,sum),'/')
c2l0_n = sweep(c2l0_,1,apply(c2l0_,1,sum),'/')
c2l0_n[c2l0_==0] = 0
# paralellization
doMC::registerDoMC(3)
# repeat N times for rank microenvironments
N = 100 # set to 100

nmfres = list()
for(n in c(5,8,10)){
  set.seed(1234)
  nmfres[[paste0('c2ln.f',n)]] = llply(1:N,function(i){nmf(c2l_n, rank = n)},.parallel = T)
  set.seed(1234)
  nmfres[[paste0('c2l0n.f',n)]] = llply(1:N,function(i){nmf(c2l0_n, rank = n)},.parallel = T)
  set.seed(1234)
  nmfres[[paste0('c2l.f',n)]] = llply(1:N,function(i){nmf(c2l_, rank = n)},.parallel = T)
  set.seed(1234)
  nmfres[[paste0('c2l0.f',n)]] = llply(1:N,function(i){nmf(c2l0_, rank = n)},.parallel = T)
}

#saveRDS(nmfres,'data.nfs/visium/ctcl/nmf.ctcl.only.rds')

nmfbst = lapply(nmfres,function(x)nmfGetBest(x,nmfNormFs$max))
nmfbst


pdf('figures/visium/ctcl/v2/nmf.consensus.dotPlots.pdf',w=10,h=8)
for(i in names(nmfbst)){
  rank = nrow(nmfbst[[i]]$coefn)
  cols = RColorBrewer::brewer.pal(min(9,rank),'Set1')
  if(rank==10)
    cols = c(cols,'gray')
  cols = setNames(cols,1:rank)
  
  plotNMFCons(coefs=nmfbst[[i]]$coefn,clcols = cols,cons = nmfbst[[i]]$cons/N,max.cex = 1.7,ylab.cex = 0.6,title=i)
}
dev.off()



cln = sapply(nmfbst,function(x){
  renameClustsBySize(x$cl[tumor.me])
})
cln = cln[order(apply(cln,1,function(x)sum(x==cln[1,]))),]
imageWithText(t(cln))




# pair corr ##############
getCellTypeCor = function(t,f,pairs,norm=TRUE,...){
  if(norm)
    t = sweep(t,1,rowSums(t),'/')
  t = split.data.frame(t,f)
  sapply(t,function(tt){
    r = apply(pairs,1,function(p){
      cor(tt[,p[1]],tt[,p[2]],...)
    })
    names(r) = apply(pairs,1,paste,collapse=':')
    r
  })
}
tumor.me = c('tumourcell','B/plasma','VE2','F2','F3','MoDC3','Th','Differentiated_KC*')
pairs = expand.grid(ct2=1:length(tumor.me),ct1=1:length(tumor.me))[,2:1]
pairs = pairs[pairs$ct1!=pairs$ct2,]
#pairs = pairs[pairs$ct1<pairs$ct2,]
pairs$ct1 = tumor.me[pairs$ct1]
pairs$ct2 = tumor.me[pairs$ct2]


f = splitSub(rownames(c2l),'|',1)
ct.pair.cor = t(getCellTypeCor(c2l,f,pairs,norm=FALSE))
stage.col = setNames(RColorBrewer::brewer.pal(9,'Set1')[c(3,4,5,1)],c('H','1A','1B','2B'))

pdf('figures/visium/ctcl/v2/ct.pair.corr.only.ctcl.pdf',w=12,h=6)
par(mfrow=c(2,4),mar=c(4,10,1,1))
for(ct in tumor.me){
  o = order(c('H'=0,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = T)
  o = o[mctcl$stage[o] != 'H']
  
  v = ct.pair.cor[mctcl$sanger_id[o],startsWith(colnames(ct.pair.cor),ct)]
  colnames(v) = sub(paste0(ct,':'),'',colnames(v))
  b = barplot(v,horiz = TRUE,beside = TRUE,col=stage.col[mctcl$stage[o]],
              border = stage.col[mctcl$stage[o]],las=1,space=c(0,2),xlab='PCC',main=ct)
  x = pmin(0,apply(v,2,min))
  y = apply(b,2,mean)
  segments(-1,y,x,y,lty=2)
}
plot.new()
legend('topleft',fill=stage.col[-1],legend=names(stage.col)[-1],bty='n')

dev.off()

pdf('figures/visium/ctcl/v2/ct.pair.corr.pdf',w=12,h=6)
par(mfrow=c(2,4),mar=c(4,10,1,1))
for(ct in tumor.me){
  o = order(c('H'=0,'1A'=1,'1B'=2,'2B'=3)[mctcl$stage],decreasing = T)
  
  
  v = ct.pair.cor[mctcl$sanger_id[o],startsWith(colnames(ct.pair.cor),ct)]
  colnames(v) = sub(paste0(ct,':'),'',colnames(v))
  b = barplot(v,horiz = TRUE,beside = TRUE,col=stage.col[mctcl$stage[o]],
              border = stage.col[mctcl$stage[o]],las=1,space=c(0,2),xlab='PCC',xlim=c(-1,1),main=ct)
  x = pmin(0,apply(v,2,min))
  y = apply(b,2,mean)
  segments(-1,y,x,y,lty=2)
}
plot.new()
legend('topleft',fill=stage.col,legend=names(stage.col),bty='n')

dev.off()
