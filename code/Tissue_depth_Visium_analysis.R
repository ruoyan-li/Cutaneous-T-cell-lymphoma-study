makeDistFeatureSampleTable = function(dist,sample,data,f,FUN=mean,min.spots=5,ncores=NULL){
  require(plyr)
  require(doMC)
  doMC::registerDoMC(ncores)
  
  sample = sample[f]
  dist = dist[f]
  data = data[f,]
  features = colnames(data)
  sids = unique(sample)
  dists = sort(unique(dist))
  
  dnames = list(dists,features,sids)
  
  ldctmtx = llply(dists,function(d){
    cat(d,',',sep='')
    res = array(NA,dim=sapply(dnames[-1],length),dimnames=dnames[-1])
    flt = dist == d
    for(f in features){
      r = tapply(data[flt,f],sample[flt],function(x)ifelse(length(x)<min.spots,NA,FUN(x)))
      res[f,names(r)] = r
    }
    res
  },.parallel = T)
  dctmtx = array(NA,dim=sapply(dnames,length),dimnames=dnames)
  for(i in 1:length(ldctmtx)){
    dctmtx[i,,] = ldctmtx[[i]]
  }
  dctmtx
}


getCTonDepthSplines = function(vs,c2loc,ct,field,min.umi=0,min.spots=30,min.dist2hf=3,df=5,norm=TRUE){
  if(norm){
    c2loc = lapply(c2loc,function(x)sweep(x,1,rowSums(x),'/'))
  }
  r = lapply(names(vs),function(n){
    if(!is.null(vs[[n]]@meta.data$dist2hf))
      f = (is.na(vs[[n]]$dist2hf) | vs[[n]]$dist2hf >= min.dist2hf) & vs[[n]]$nCount_Spatial >= min.umi
    else
      f = TRUE
    x = vs[[n]]@meta.data[f,field]
    y = c2loc[[n]][f,ct]
    
    #x[vs[[n]]$man.ann[f]=='dermis'] = -x[vs[[n]]$man.ann[f]=='dermis']
    
    p = NULL
    tryCatch({
      p=predict(smooth.spline(x,y,df=df))
    },error=function(e){})
    if(sum(f)<min.spots)
      p = NULL
    if(!is.null(p)){
      p$y[p$y<0] = 0
      if(norm)
        p$y[p$y>1] = 1
    }
    list(pred=p,raw=data.frame(sid=rep(n,length(x)),x=x,y=y))
  })
  names(r) = names(vs)
  r
}


makeCTonDepthAnalysis = function(vsf,c2l,cts,field,state,min.umi=0,min.spots=30,min.dist2hf=3,df=5,norm=TRUE){
  r = lapply(cts,function(ct){
    t = getCTonDepthSplines(vsf,c2l,ct,field,min.umi=min.umi,min.spots = min.spots,min.dist2hf=min.dist2hf,df = df,norm = norm)
    a2d = do.call(rbind,lapply(t,'[[','raw'))
    a2d$State = setNames(state,names(vsf))[a2d$sid]
    a2d = split(a2d,a2d$State)
    a2ds = lapply(a2d,function(x){
      r = do.call(rbind,lapply(split(x$y,round(x$x)),function(z){
        data.frame(n=length(z),mean=mean(z),median=median(z),sd=sd(z)/sqrt(length(z)))
      }))
      r$x=as.numeric(rownames(r))
      r
    })
    t$a2d.raw = a2d
    t$a2ds.raw = a2ds
    t
  })
  names(r) = cts
  r
}

testDFConditions = function(m,f1,f2){
  pv = m[,,1]
  pv[,] = NA
  diff = pv
  for(i in 1:dim(m)[1]){
    for(j in 1:dim(m)[2]){
      x1 = m[i,j,f1]
      x2 = m[i,j,f2]
      if(sum(!is.na(x1))<2 || sum(!is.na(x2))<2){
        next
      }
      tryCatch({
        r = t.test(x1,x2)$p.value
        pv[i,j] = r
        diff[i,j] = mean(x2,na.rm=T) - mean(x1,na.rm=T)
      },error=function(e){})
    }
  }
  r = list(pv=pv,diff=diff)
  r$fdr = r$pv
  r$fdr[,] = p.adjust(r$pv,'fdr')
  r
}

plotDF.HM = function(m,cnd,cond.pair,l2fc.zlim=NULL,log.pseudocount=0,celltype.class=NULL,cols= hcl.colors(100, "YlOrRd", rev = TRUE)){
  comp = testDFConditions(log2(m+log.pseudocount),f1=cnd==cond.pair[1],f2=cnd==cond.pair[2])
  l2fc = comp$diff
  
  pvt = ifelse(comp$fdr > 0.2,'',ifelse(comp$fdr < 0.05,'*','.'))
  
  m1 = apply(m[,,cnd==cond.pair[1]],1:2,mean,na.rm=TRUE)
  m2 = apply(m[,,cnd==cond.pair[2]],1:2,mean,na.rm=TRUE)
  m1 = sweep(m1,2,apply(m1,2,max,na.rm=TRUE),'/')
  m2 = sweep(m2,2,apply(m2,2,max,na.rm=TRUE),'/')
  #m1[is.na(m1)] = 0
  #m2[is.na(m2)] = 0
  xlim = c(0.5,0.5+dim(m)[1])
  if(!is.null(celltype.class))
    xlim = c(-0.5,0.5+dim(m)[1])
  imageWithText(m1,'',main=cond.pair[1],xlab='Distance to epidermis-dermis interface (spots)',xlim=xlim,col=cols)
  plotColorAnn(1:dim(m)[2],c(-0.5,0.3),list(class=celltype.class),horis=F)
  imageWithText(m2,'',main=cond.pair[2],xlab='Distance to epidermis-dermis interface (spots)',xlim=xlim,col=cols)
  plotColorAnn(1:dim(m)[2],c(-0.5,0.3),list(class=celltype.class),horis=F)
  if(is.null(l2fc.zlim)){
    l2fc.zlim=range(l2fc,na.rm=T)
  }
  l2fc[l2fc<l2fc.zlim[1]] = l2fc.zlim[1]
  l2fc[l2fc>l2fc.zlim[2]] = l2fc.zlim[2]
  zmax=max(abs(l2fc.zlim))
  imageWithText(l2fc,pvt,zlim=c(-zmax,zmax),main=paste0('log2(',cond.pair[2],'/',cond.pair[1],') (fdr>0.1)'),
                xlab='Distance to epidermis-dermis interface (spots)',xlim=xlim)
  plotColorAnn(1:dim(m)[2],c(-0.5,0.3),list(class=celltype.class),horis=F)
  
  if(any(!is.na(l2fc))){
    zlim=range(l2fc,na.rm=T)
    plotColorLegend2(1,1.2,0.2,0.8,c(-zmax,zmax),zlim,identity,function(x)num2col(x),title='log2FC')
  }
  legend(grconvertX(1,'npc','user'),grconvertY(1,'npc','user'),xpd=NA,pch=c('*','.'),legend=c('fdr<0.05','fdr<0.2'),bty='n')
}

plotColorAnn = function(centers,width.range,labs,cols=NULL,horis = FALSE,text.cex=0.7){
  # imageWithText(matrix(1:10,ncol=2),ylim=c(-1,2.5),xlim=c(-2,5.5))
  # plotColorAnn(1:5,c(-1,0.4),list(a1=c(1,1,2,3,3),a2=c('a','b','a,','b','a')),horis=T)
  # plotColorAnn(1:2,c(-2,0.4),list(y1=c(1,2),y2=c('a','a')),horis=F)
  if(is.null(cols)){
    cols = lapply(labs,char2col)
  }
  l = abs(centers[2]-centers[1])
  w = (width.range[2] - width.range[1])/length(labs)
  centers = centers - l/2
  
  for(i in seq_along(labs)){
    if(horis){
      rect(centers,width.range[1]+(i-1)*w,centers+l,width.range[1]+i*w,border=NA,col=cols[[i]][labs[[i]]])
      text(max(centers)+l,width.range[1]+(i-0.5)*w,names(labs)[i],adj=c(0,0.5),xpd=NA,cex=text.cex)
    }else{
      rect(width.range[1]+(i-1)*w,centers,width.range[1]+i*w,centers+l,border=NA,col=cols[[i]][labs[[i]]])
      text(width.range[1]+(i-0.5)*w,min(centers),names(labs)[i],adj=c(1,0.5),srt=90,xpd=NA,cex=text.cex)
    }
  }
}