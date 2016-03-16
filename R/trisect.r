#hopkins stat
cluster.tendency = function(models, ctype='knn'){
  require(fields)
  hopkins = function(t, alt='cox'){ #distance from random points
    marked = sample(n, m); Di = mydata[marked,]; 
    Ri = sapply(1:length(low), function(i) runif(m, low[i], high[i]))
    wis = apply(rdist(Di, mydata[setdiff(1:n, marked), ]), 1, min)
    uis = apply(rdist(Ri, mydata), 1, min)
    if(alt=='hopkins'){
      H = sum(uis)/(sum(uis)+sum(wis));
    }else{#cox-lweis statistics
      H = sum(uis/wis)/m;
    }
    return(H)
  }
  coxlist = list(); t = 1000
  #fot tf we used 25 to 75 percentile data
  ld = load.weight.matrix(models); mweights = ld$mweights
  n = nrow(mweights); m = floor(n/5); 
  mydata = apply(mweights, 2, function(x) mean(x)+(x-mean(x))/sd(x) )
  low = apply(mydata, 2, quantile, probs=0.25); high = apply(mydata, 2, quantile, probs=0.75); 
  #low = apply(mydata, 2, sd); high = apply(mydata, 2, max) - apply(mydata, 2, sd)
  #centroid = apply(mydata, 2, mean); radius = apply(mydata, 2, max)/2
  #mweights = mweights[,apply(mweights, 2, quantile, probs=c(0.95))>0]
  mydata = mweights
  low = apply(mydata, 2, sd); high = apply(mydata, 2, max) - apply(mydata, 2, sd)
  cox = sapply(1:t, hopkins, alt='hopkins'); coxlist[[tf]] = cox; mean(cox)
  print(paste(tf, round(mean(mean(cox)+(cox-mean(cox))/sd(cox)), 2), sum(cox<=0.6)/t))
  #save(coxlist, file=paste0(root, project, 'hopkins-tf.Rda'))
  par(mar=c(8,4,3,3)); boxplot(coxlist, las=2, ylim=c(0.45, 1), ylab='Hopkins statistics (H)'); abline(h=0.6, col='blue')
  
  hopkinslist = list(); t = 500
  #for tf-cell we used sd to (max-sd)
  ld = load.weight.matrix(models); mweights = ld$mweights 
  start = 0; ends = cumsum(sapply(models, length)); limit = min(sapply(models, length))-1
  tissues = names(models); hopkinslist[[tf]] = list()
  for(ti in 1:length(tissues)){
    mydata <- mweights[(start+1):as.numeric(ends[ti]),]; n = nrow(mydata); m = floor(n/5); 
    tmp = apply(mydata, 2, sum); mydata = mydata[,names(mydata)[tmp>0]]
    #mydata = apply(mydata, 2, function(x) mean(x)+(x-mean(x))/sd(x) )
    #low = apply(mydata, 2, quantile, probs=0.25); high = apply(mydata, 2, quantile, probs=0.75); 
    low = apply(mydata, 2, sd); high = apply(mydata, 2, max) - apply(mydata, 2, sd)
    cox = sapply(1:t, hopkins, alt='hopkins'); hopkinslist[[tf]][[tissues[ti]]] = cox;
    print(paste(tf, tissues[ti], mean(cox), sum(cox<=0.5)/t))
    start = as.numeric(ends[ti])
  }
  
  #save(hopkinslist, file=paste0(root, project, 'hopkins-tissue.Rda'))
  boxplot(unlist(sapply(hopkinslist, function(x) sapply(x, mean))))
  
  coxlist[['all.TF-Cell.models']] = unlist(sapply(hopkinslist, function(x) sapply(x, mean)))
  p.value= unlist(sapply(hopkinslist, function(x) sapply(x, function(cox) sum(cox<=0.5)/t)))
  p.value[p.value==0] = 0.001
  boxplot(-log10(p.value), ylim=c(1, 3.2)); abline(h=1.3, col='blue')
}

decide.submodeling.dudahart = function(root, project, tfs, size=1000){
  require(fpc)
  
  tmpf = function(i, mydata){
    fit = kmeans(mydata, 2); wss = dudahart2(mydata, fit$cluster) 
    return(c(round(wss$dh,3), wss$p.value, wss$cluster1))
  }
  count1 = 0; count2 = 0; count3 = 0
  df = data.frame(matrix(NA, nrow=0, ncol=5)); colnames(df) = c('TF', 'Cell_line', 'dh', 'P.value', 'Homogenous')
  for(tf in tfs){
    print(tf)
    tmp = load.weight.matrix(root, project, tf); reduced_mw = tmp$reduced_mw; models = tmp$models; tissues = names(models)
    start = 0; ends = cumsum(sapply(models, length)); limit = min(sapply(models, length))-1
    for(ti in 1:length(tissues)){
      print(tissues[ti])
      mydata = reduced_mw[(start+1):as.numeric(ends[ti]), ]; 
      tmp = apply(mydata, 2, sum); mydata = mydata[,names(mydata)[tmp>0]]
      tmp = sapply(1:size, tmpf, mydata); 
      if(sum(tmp[2,]<0.05)!=size){ warning('sometimes tissue level is not significant') }
      df = rbind(df, data.frame(TF=tf, Cell_line=tissues[ti], dh=median(tmp[1,]), P.value=median(tmp[2,]), 
                                Homogenous=sum(tmp[3,])))
      start = as.numeric(ends[ti]); 
      count1 = count1+sum(tmp[2,]>=0.001)
      count2 = count2+sum(tmp[2,]>=0.01)
      count3 = count3+sum(tmp[2,]>=0.05)
    }
  }
  
  dhlist = list(); dhlist[['all.TF-Cell models']] = as.numeric(df$dh)
  tfdf = data.frame(matrix(NA, nrow=0, ncol=5)); colnames(tfdf) = c('TF', 'dh', 'P.value', 'Homogenous')
  for(tf in tfs){
    print(tf)
    tmp = load.weight.matrix(root, project, tf); reduced_mw = tmp$reduced_mw; mydata <- reduced_mw
    tmp = sapply(1:size, tmpf, mydata); 
    if(sum(tmp[2,]<0.05)!=size){ warning('sometimes tf level is not significant') }
    dhlist[[tf]] = tmp[1,]
    tfdf = rbind(tfdf, data.frame(TF=tf, dh=mean(tmp[1,]), P.value=mean(tmp[2,]), Homogenous=sum(tmp[3,])))
  }
  
  par(mar=c(8,4,3,3)); boxplot(dhlist, ylim=c(0,1.2), las=2, ylab='dh ratio'); abline(h=1, col='blue')
  ld = list(tissue=df, tf=tfdf, l = dhlist)
  #save(ld, file=paste0(root, project, 'dhlist.Rda') )
  print(paste('counts: ', count1, count2, count3))
  return(list(tissue=df, tf=tfdf, l = dhlist))
}

ggplot.sequence.overlap = function(){
  require(plotflow); require(RColorBrewer); require(reshape2)
  load('tissueoverlap.sameclust.rda')
  pallete = brewer.pal(12,"Paired");
  dflist2 = sapply(dflist, function(x) round(sweep(as.array(x), 2, apply(x, 2, sum), '/')*100, 2))
  tmp = list()
  for(tf in tfs){
    mdf = melt(dflist2[[tf]]) ; #mdf = cbind(mdf, lab=as.vector(sapply(apply(dflist[[tf]], 2, sum), rep, 4)))
    mdf = cbind(unique(mdf[,1:2]), lab=as.vector(sapply(apply(dflist[[tf]], 2, sum), rep, 4)))
    tmp[[tf]] = ggplot(data=melt(dflist2[[tf]]), aes(x=Var2, y=value, fill=Var1))+
      geom_bar(stat="identity")+ggtitle(tf)+
      #geom_text(data=melt(apply(dflist[[tf]], 2, sum)), aes(y=101,label=value), vjust=0)+
      geom_text(data=mdf, aes(y=101, label=lab), vjust=0)+
      scale_fill_manual(name='', values = pallete[c(8,10,7,9)])+
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
            legend.position="top",
            axis.title.x=element_blank(), axis.title.y=element_blank())
    #theme_minimal()
  }
  #merge_pdf(2, 'Fig S3 (Sequence overlap).pdf')
  size = sapply(dflist, ncol)+0.5; size[size<5] = 5
  merge_pdf(length(tfs), file = "bar.pdf", widths=size,#sapply(dflist, ncol)+0.5, 
            heights=rep(6, length(tfs)))
  print(tmp[[tfs[1]]]); 
  print(tmp[[tfs[2]]]); 
  print(tmp[[tfs[3]]]); 
  print(tmp[[tfs[4]]])
  print(tmp[[tfs[5]]]); 
  print(tmp[[tfs[6]]]); 
  print(tmp[[tfs[7]]]); 
  print(tmp[[tfs[8]]])
  print(tmp[[tfs[9]]]); 
  print(tmp[[tfs[10]]]); 
  print(tmp[[tfs[11]]]); 
  print(tmp[[tfs[12]]])
  print(tmp[[tfs[13]]]); 
  print(tmp[[tfs[14]]]); 
  print(tmp[[tfs[15]]]); 
  print(tmp[[tfs[16]]])
  print(tmp[[tfs[17]]]); 
  print(tmp[[tfs[18]]]); 
  print(tmp[[tfs[19]]]); 
  print(tmp[[tfs[20]]])
  print(tmp[[tfs[21]]]); 
  print(tmp[[tfs[22]]]); 
  print(tmp[[tfs[23]]]); 
  
}

shuffle.targets = function(root, project, tfs, ctype, clen, qi=1, xweight=1, overwrite=T, randomized=2){
  require('ChIPpeakAnno')
  require('org.Hs.eg.db')
  
  if(ctype=='som'){
    ctype = paste0(xweight*100, ctype)
  }
  if(randomized==1){
    load(paste0(root, project, 'fulltargetmats.Rda'))
  }else if(randomized==2){
    load(paste0(root, project, 'unifiedtargetmats.Rda'))
  }else{
    load(file = paste0(root,'rnaseq/A549/', 'countmat.Rda')); countmat = as.matrix(countmat$counts);
    fullsize = nrow(countmat); candidates = rownames(countmat)
  }
  for(tf in tfs){
    tfpath = paste0(root, project, tf); print(tf)
    load(paste0(tfpath, '/targetgenes_', ctype, clen, 'qi', qi, '.Rda')) 
    #currentmat is the unshuffled targetmat
    currentmat = targets$targetmat; tissues = dir(tfpath)[file.info(dir(tfpath, full.names=T))$isdir]
    map = sapply(currentmat, function(x) sapply(x, function(p) length(p@values@unlistData$feature)))
    targetmat = currentmat; ensemblmat = targets$ensemblmat; uniprotmat = targets$uniprotmat
    for(tissue in tissues){
      clabels = which(sapply(map, function(x) length(grep(tissue, names(x)))>0))
      ccount = sapply(clabels, function(x) map[[x]][grep(tissue, names(map[[x]]))]); print(ccount)
      cloc = unlist(sapply(clabels, function(x) grep(tissue, names(map[[x]]))))
      for(ci in 1:length(clabels)){
        if(randomized==1){#
          fullsize = nrow(fulltargetmats[[tf]][[tissue]]);
          targetmat[[clabels[ci]]][[cloc[ci]]] = fulltargetmats[[tf]][[tissue]][sort(sample(fullsize, ccount[ci])),]
        }else if(randomized==2){#random targets
          fullsize = nrow(unifiedtargetmats[[tf]]);
          targetmat[[clabels[ci]]][[cloc[ci]]] = unifiedtargetmats[[tf]][sort(sample(fullsize, ccount[ci])),]
        }else{
          break;
        }
        ensemblmat[[clabels[ci]]][[cloc[ci]]] = targetmat[[clabels[ci]]][[cloc[ci]]]@values@unlistData$feature
      }
    }
    entrezmat = targets$entrezmat; 
    for(ci in 1:clen){
      tissues = names(map[[ci]])
      for(tissue in tissues){
        if(randomized==3 | randomized==4){#3-random genes, 4-random genes but removing the foreground
          if(length(entrezmat[[ci]][[tissue]])>0){
            if(randomized==4){
              candidates = setdiff(rownames(countmat), entrezmat[[ci]][[tissue]]); fullsize = length(candidates); 
            }
            ccount = sum(entrezmat[[ci]][[tissue]] %in% rownames(countmat))
            entrezmat[[ci]][[tissue]] = candidates[sample(fullsize, ccount)]
          }
        }else{
          tmp = convert2EntrezID(ensemblmat[[ci]][[tissue]], orgAnn='org.Hs.eg.db', ID_type='ensembl_gene_id')
          if(length(tmp)>0){
            entrezmat[[ci]][[tissue]] = tmp[!is.na(tmp)] 
          }else{
            entrezmat[[ci]][[tissue]] = tmp
          } 
        }
      }
    }
    targets = list(ensemblmat=ensemblmat, entrezmat=entrezmat, targetmat=targetmat, uniprotmat=uniprotmat)
    if(overwrite==T){
      save(targets, file=paste0(tfpath, '/targetgenes_', ctype, clen, 'qi', qi, 'shuffled',randomized,'.Rda')) 
    }
  }
}

