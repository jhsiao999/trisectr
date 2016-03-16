
get.modelscore = function(submodels, test, best.iter, vtype='response'){
  submodels.pr <- predict(submodels, test[-1], best.iter, type=vtype);
  
  return(submodels.pr)
}

get.modelsize = function(model){
  options(warn=-1)
  if(sum(is.finite(model$oobag.improve))==0){
    best.iter = sum(model$oobag.improve>0) 
  }else{
    infval = max(model$oobag.improve[is.finite(model$oobag.improve)])
    if(infval<0){
      infval = 0
    }else{
      infval = infval*2
    }
    model$oobag.improve[is.infinite(model$oobag.improve)&model$oobag.improve>0] = infval
    model$oobag.improve[is.infinite(model$oobag.improve)&model$oobag.improve<0] = -infval
    best.iter <- gbm.perf(model, method="OOB", plot.it=F)
  }
  options(warn=0)
  return(best.iter)
}

get.performance = function(submodels.pred, digits=3, verbose=T){
  auc = performance(submodels.pred,"auc")@y.values[[1]]
  if(verbose==T){
    print(auc) 
  }
  return(round(auc, digits=digits))
}

#how do I allow the existing users to access parameters from the gbm function; using ...
cellSpecClsfr = function(dataset, distb='huberized', pwm.ids=NULL, verbose=T, interaction.depth=15, 
                         bag.fraction=0.5, shrinkage=0.05, n.minobsinnode=30, n.trees = 100, 
                         ...){
  require('ROCR')
  require('gbm')
  
  if(!is.null(pwm.ids)){
    train = train[ ,c('label', pwm.ids)]; test = test[ ,c('label', pwm.ids)]
  }
  submodels = gbm.fit(dataset$train[-1], dataset$train$label, distribution=distb, n.trees = 100, 
                   interaction.depth=interaction.depth, n.minobsinnode=30, shrinkage=shrinkage, 
                   bag.fraction=bag.fraction, verbose=F, keep.data=F)
  
  best.iter = get.modelsize(submodels); 
  submodels$n.trees = best.iter
  submodels$trees = submodels$trees[1:best.iter] 
  submodels$train.error = submodels$train.error[1:best.iter] 
  submodels$valid.error = submodels$valid.error[1:best.iter]
  submodels$oobag.improve = submodels$oobag.improve[1:best.iter]
  
  modelscores = get.modelscore(submodels, dataset$test, best.iter); 
  model.pred = prediction(modelscores, dataset$test$label)
  auc = get.performance(model.pred)
  
  if(verbose==T){
    print(paste(best.iter, auc)) 
  }
  
  return(list(model=submodels, auc=auc))
}

feature.importance = function (object, n.trees=object$n.trees, single.tree=F, normalize=T) {
  if (n.trees < 1) {
    stop("n.trees must be greater than 0.")
  }
  if (n.trees > object$n.trees) {
    warning("Exceeded total number of GBM terms. Results use n.trees=", object$n.trees, " terms.\n")
    n.trees <- object$n.trees
  }
  rel.inf <- relative.inf(object, n.trees, single.tree=single.tree) #modified for single tree
  rel.inf[rel.inf < 0] <- 0
  if (normalize) 
    rel.inf <- 100 * rel.inf/sum(rel.inf)
  
  return(data.frame(var = object$var.names, rel.inf = rel.inf)) 
}

relative.inf = function (object, n.trees, single.tree=F) {
  if (missing(n.trees)) {
    if (object$train.fraction < 1) {
      n.trees <- gbm.perf(object, method = "test", plot.it = FALSE)
    }
    else if (!is.null(object$cv.error)) {
      n.trees <- gbm.perf(object, method = "cv", plot.it = FALSE)
    }
    else {
      n.trees <- length(object$train.error)
    }
    cat(paste("n.trees not given. Using", n.trees, "trees.\n"))
  }
  get.rel.inf <- function(obj) {
    lapply(split(obj[[6]], obj[[1]]), sum)
  }
  if(single.tree==T){
    temp <- unlist(lapply(object$trees[n.trees:n.trees], get.rel.inf))
  }else{
    temp <- unlist(lapply(object$trees[1:n.trees], get.rel.inf))
    if(length(grep('[.]', names(temp)))>0){
      names(temp) = unlist(strsplit(names(temp), '[.]'))[seq(2,(2*length(temp)),2)]
    }
  }
  rel.inf.compact <- unlist(lapply(split(temp, names(temp)), sum))
  rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) != "-1"]
  rel.inf <- rep(0, length(object$var.names))
  i <- as.numeric(names(rel.inf.compact)) + 1
  rel.inf[i] <- rel.inf.compact
  names(rel.inf) <- object$var.names
  return(rel.inf = rel.inf)
}

load.weight.matrix = function(models){
  require('gbm')
  
  cells = names(models)
  col.cells = unlist(sapply(cells, function(x) rep(x, models[[x]]$n.trees)))
  #oobag.improve = sapply(models, function(x) x$oobag.improve)
  mweights = sapply(models, function(x) sapply(1:x$n.trees, 
                                               function(y) feature.importance(x, n.trees=y, single.tree=T)$rel.inf))
  
  mweights = t(Reduce(cbind, mweights)); colnames(mweights) = models[[cells[1]]]$var.names
  fvals = apply(mweights, 2, sum); mweights = mweights[,colnames(mweights)[fvals>0]]
  return(list(mweights=mweights, col.cells=col.cells))#, oobag.improve=oobag.improve))
}

get.cluster.membership = function(models, clen, ctype='knn', plotit=T){
  require('kohonen')
  require('gplots')
  
  ld = load.weight.matrix(models); mweights = ld$mweights
  if(ctype=='knn'){
    fit = kmeans(mweights, clen) 
    cluster.membership = as.data.frame.matrix(table(ld$col.cells, fit$cluster))
  }else{
    fit = xyf(mweights, classvec2classmat(ld$col.cells), xweight=xweight,
              somgrid(4, clen/4, "hexagonal"), rlen=500)
    cluster.membership = as.data.frame.matrix(table(ld$col.cells, fit$unit.classif))
  }
  
  if(plotit==T){
    plot.cluster.membership(cluster.membership, ctype=ctype)
  }
  return(fit)
}

plot.cluster.membership = function(cluster.membership, ctype='knn'){
  require('RColorBrewer')
  require('gplots')
  
  colbreaks = c(seq(0,0.2,length=100), seq(0.2001,1,length=100), 
                seq(1.0001,max(cluster.membership, na.rm=T),length=100))
  pallete = colorRampPalette(brewer.pal(9,"YlGn")[c(1,5,8)])(n = 299)
  #pdf(paste0(getwd(), '/cluster.membership.pdf'))
  heatmap.2(t(as.matrix(cluster.membership)),  main = ctype, 
            cellnote=t(as.matrix(cluster.membership)), cexRow=1, notecol="black", 
            density.info="none", trace="none", margins =c(6,7), col=pallete,  dendrogram='col',        
            cexCol=1.2,  Rowv=F, breaks=colbreaks, keysize=0.8)
  #dev.off()
}

make.cluster.ensembles = function(models, fit, clen, ctype='knn'){
  require('gbm')
  
  ld = load.weight.matrix(models); mweights = ld$mweights
  if(is.null(fit)){
    fit = get.cluster.membership(models, clen=clen, ctype=ctype, plotit=F)
  }
  if(ctype=='knn'){
    cluster_indices = fit$cluster
  }else{
    cluster_indices = fit$unit.classif
  }
  cells = names(models)
  allmodels = unlist(models, recursive=F)
  clusters = list();
  for(ci in 1:clen){
    flen = sum(cluster_indices==ci); clust.clsfr = list(); 
    if(flen==0){
      clusters[[ci]] = clust.clsfr; next
    }
    weight = rep(1, length(models)); names(weight) = cells
    tmp = table(ld$col.cells[cluster_indices==ci]); weight[names(tmp)] = tmp
    clust.clsfr$initF = sum(sapply(cells, function(x) models[[x]]$initF*weight[x]))/sum(weight)
    clust.clsfr$fit = runif(floor(sum(sapply(models, function(x) x$nTrain))/length(models)), min=-1); 
    clust.clsfr$train.error = unlist(sapply(models, function(x) x$train.error))[cluster_indices==ci];
    clust.clsfr$valid.error = rep(NaN, flen); 
    clust.clsfr$oobag.improve = unlist(sapply(models, function(x) x$oobag.improve))[cluster_indices==ci]
    clust.clsfr$trees = list(); clust.clsfr$c.splits = list(); clust.clsfr$bag.fraction = 0.5; 
    clust.clsfr$distribution$name = models[[cells[1]]]$distribution$name; 
    clust.clsfr$interaction.depth = models[[cells[1]]]$interaction.depth; 
    clust.clsfr$n.minobsinnode = models[[cells[1]]]$n.minobsinnode; 
    clust.clsfr$num.classes = models[[cells[1]]]$num.classes; clust.clsfr$n.trees = flen; 
    clust.clsfr$nTrain = as.integer(floor(sum(sapply(models, function(x) x$nTrain))/length(models)));
    clust.clsfr$train.fraction = models[[cells[1]]]$train.fraction; 
    clust.clsfr$response.name = models[[cells[1]]]$response.name; 
    clust.clsfr$shrinkage = models[[cells[1]]]$shrinkage; 
    clust.clsfr$var.levels = models[[cells[1]]]$var.levels; 
    #for(i in 1:ncol(reduced[-1])){
    #clust.clsfr$var.levels[[i]] <- quantile(reduced[, i+1], prob = (0:10)/10, na.rm = TRUE) 
    #}
    clust.clsfr$var.monotone = models[[cells[1]]]$var.monotone; 
    clust.clsfr$var.names = models[[cells[1]]]$var.names; 
    clust.clsfr$var.type = models[[cells[1]]]$var.type; 
    clust.clsfr$verbose = models[[cells[1]]]$verbose;
    clust.clsfr$trees = unlist(sapply(models, function(x) x$trees), recursive=F)[cluster_indices==ci]
    class(clust.clsfr) = 'gbm'; clusters[[ci]] = clust.clsfr
  }#end of cluster
  
  return(clusters)
}

get.targets = function(datalist, clusters, clen){
  require('ChIPpeakAnno')
  require('org.Hs.eg.db')
  require('gbm')
  require('biomaRt')
  
  sequences = lapply(datalist, function(x) rbind(x$train[x$train$label==1,], x$test[x$test$label==1,]))
  scoretab = sapply(1:clen, function(x){ 
    if(class(clusters[[x]])!='gbm'){
      rep(NA, sum(sapply(sequences, nrow)))
    }else{
      unlist(sapply(cells, function(y){
        if(sum(grepl(y, names(clusters[[x]]$trees)))==0){
          return(rep(NA, nrow(sequences[[y]])))
        }else{
          return(get.modelscore(clusters[[x]], sequences[[y]], clusters[[x]]$n.trees))
        }
      }))
    }
  })
  status = scoretab>=1; #positive sequences deemed by the cluster
  sequences = lapply(datalist, function(x) rbind(x$trainseq[x$trainseq$label==1,],
                                                 x$testseq[x$testseq$label==1,]))
  col.cells = unlist(sapply(cells, function(x) rep(x, nrow(sequences[[x]]))))
  sequences = Reduce(rbind, sequences)
  mart = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  cells = names(datalist); data(TSS.human.GRCh37)
  ensemblmat = list(); entrezmat = list(); targetmat = list(); 
  for(ci in 1:clen){
    ensemblmat[[ci]] = list(); entrezmat[[ci]] = list(); targetmat[[ci]] = list();
    for(cell in cells){
      cell.seq = sequences[col.cells==cell & rowSums(status, na.rm=T)==1 & status[,ci]==T & 
                             !is.na(status[,ci]), ]
      if(nrow(cell.seq)==0){
        next
      }
      peaks = RangedData(IRanges(start = cell.seq$start, end = cell.seq$end),
                         space = sapply(strsplit(as.character(cell.seq$seqnames), 'chr'), '[[', 2));
      annotatedPeaks = annotatePeakInBatch(peaks, AnnotationData = TSS.human.GRCh37, output='nearestStart',
                                           PeakLocForDistance='middle', FeatureLocForDistance = "TSS"); head(annotatedPeaks)
      targetmat[[ci]][[cell]] = annotatedPeaks
      ensemblmat[[ci]][[cell]] = annotatedPeaks@elementMetadata@listData$feature;
      entrezmat[[ci]][[cell]] = convert2EntrezID(ensemblmat[[ci]][[cell]], orgAnn='org.Hs.eg.db', ID_type='ensembl_gene_id')
      if(length(entrezmat[[ci]][[cell]])>0)
        entrezmat[[ci]][[cell]] = entrezmat[[ci]][[cell]][!is.na(entrezmat[[ci]][[cell]])]
    }
  }
  targets = list(ensemblmat=ensemblmat, entrezmat=entrezmat, targetmat=targetmat)
  return(targets)
}

setup.exp.coherence = function(targets, exprsn, exptheK=1, verbose=T){
  entrezmat = targets$entrezmat; clen = length(targets$entrezmat);
  cells = intersect(unique(unlist(sapply(entrezmat, names))), colnames(exprsn))
  coexpressed = function(gi, gj, t1, t2){
    return((exprsn[gi,t1]>=exptheK & exprsn[gj,t2]>=exptheK))
  }
  
  expression.cluster = sapply(1:clen, function(z){ if(verbose==T) print(paste('cluster#', z))
    unlist(sapply(cells, function(x) unlist(sapply(cells, function(y){
      if(x<y){
        entrezmat[[z]][[x]] = entrezmat[[z]][[x]][entrezmat[[z]][[x]] %in% rownames(exprsn)]
        entrezmat[[z]][[y]] = entrezmat[[z]][[y]][entrezmat[[z]][[y]] %in% rownames(exprsn)]
        from.cluster = expand.grid(entrezmat[[z]][[x]], entrezmat[[z]][[y]])
        from.cluster$Var1 = as.character(from.cluster$Var1); 
        from.cluster$Var2 = as.character(from.cluster$Var2)
        from.cluster = from.cluster[from.cluster$Var1!=from.cluster$Var2,]
        if(nrow(from.cluster)!=0){
          tmp = sapply(1:nrow(from.cluster), function(w) 
            coexpressed(from.cluster$Var1[w], from.cluster$Var2[w], x, y))
          names(tmp) = apply(from.cluster, 1, paste, collapse=';'); tmp
        }
      }
    }))))
  })
  
  #cells = intersect(unique(unlist(sapply(entrezmat, names))), colnames(exprsn))
  expression.cell = sapply(cells, function(z){ if(verbose==T) print(z)
    unlist(sapply(1:clen, function(x) unlist(sapply(1:clen, function(y){
      if(x<y){
        entrezmat[[x]][[z]] = entrezmat[[x]][[z]][entrezmat[[x]][[z]] %in% rownames(exprsn)]
        entrezmat[[y]][[z]] = entrezmat[[y]][[z]][entrezmat[[y]][[z]] %in% rownames(exprsn)]
        from.cell = expand.grid(entrezmat[[x]][[z]], entrezmat[[y]][[z]])
        from.cell$Var1 = as.character(from.cell$Var1); 
        from.cell$Var2 = as.character(from.cell$Var2)
        from.cell = from.cell[from.cell$Var1!=from.cell$Var2,]
        if(nrow(from.cell)!=0){
          tmp = sapply(1:nrow(from.cell), function(w) 
            coexpressed(from.cell$Var1[w], from.cell$Var2[w], z, z))
          names(tmp) = apply(from.cell, 1, paste, collapse=';'); tmp
        }
      }
    }))))
  })
  
  return(list(cluster=expression.cluster, cell=expression.cell))
}

exp.coherence = function(coherence){#do i need target
  expression.cell = coherence$cell; expression.cluster = coherence$cluster
  from.cell = unlist(expression.cell); c = sum(from.cell); d = length(from.cell)-c;
  df = matrix(NA, nrow = length(expression.cluster), ncol = 2); colnames(df) = c('odds', 'p.value')
  for(ci in 1:length(expression.cluster)){
    if(length(expression.cluster[[ci]])==0 | length(unlist(expression.cluster[[ci]]))<=1 ) next# 
    from.cluster = unlist(expression.cluster[[ci]]); a = sum(from.cluster); b = length(from.cluster)-a; 
    ftab = matrix(c(a, b, c, d),nrow=2,dimnames=list(paired = c("yes","no"), from.acc = c("cluster","cell")));
    ftest = fisher.test(ftab, alternative = "two.sided")
    df[ci, 'odds'] = round(ftest$estimate, 3); df[ci, 'p.value'] = round(ftest$p.value, 100)
  }
  
  return(data.frame(df))
}

setup.pathway.coherence = function(genemat, pathwaymat, verbose=T){
  clen = length(genemat); cells = sort(unique(unlist(sapply(genemat, names))))
  
  interacting = function(df, i){ 
    cond = df$Var1[i] %in% rownames(pathwaymat) & df$Var2[i] %in% colnames(pathwaymat); 
    if(cond){
      return(pathwaymat[df$Var1[i], df$Var2[i]])
    }else{
      return(0)
    }
  }
  
  pathway.cluster = sapply(1:clen, function(z){ if(verbose==T) print(paste('cluster#', z))
    unlist(sapply(cells, function(x) unlist(sapply(cells, function(y){ 
      if(x<y){
        #if(length(unique(genemat[[z]][[x]]))>=length(genemat[[z]][[x]]) & 
         #  length(unique(genemat[[z]][[y]]))>=length(genemat[[z]][[y]])){ print(paste(z,x,y));
          from.cluster = expand.grid(genemat[[z]][[x]], genemat[[z]][[y]])
          from.cluster$Var1 = as.character(from.cluster$Var1); 
          from.cluster$Var2 = as.character(from.cluster$Var2)
          from.cluster = from.cluster[from.cluster$Var1!=from.cluster$Var2,]
          if(nrow(from.cluster)!=0){
            tmp = sapply(1:nrow(from.cluster), function(w) interacting(from.cluster, w)) 
            names(tmp) = apply(from.cluster, 1, paste, collapse=';'); tmp
          }
        #}
      }
    }))))
  }) 
  
  pathway.cell = sapply(cells, function(z){ if(verbose==T) print(z)
    unlist(sapply(1:clen, function(x) unlist(sapply(1:clen, function(y){
      if(x<y){
        #if(length(unique(genemat[[z]][[x]]))>=length(genemat[[z]][[x]]) & 
         #  length(unique(genemat[[z]][[y]]))>=length(genemat[[z]][[y]])){
          from.cell = expand.grid(genemat[[x]][[z]], genemat[[y]][[z]])
          from.cell$Var1 = as.character(from.cell$Var1); 
          from.cell$Var2 = as.character(from.cell$Var2)
          from.cell = from.cell[from.cell$Var1!=from.cell$Var2,]
          if(nrow(from.cell)!=0){
            tmp = sapply(1:nrow(from.cell), function(w) interacting(from.cell, w)) 
            names(tmp) = apply(from.cell, 1, paste, collapse=';'); tmp
          }
        #}
      }
    }))))
  }) 
  
  return(list(cluster=pathway.cluster, cell=pathway.cell))
}

pathway.coherence = function(coherence){
  pathway.cell = coherence$cell; pathway.cluster = coherence$cluster
  from.cell = unlist(pathway.cell); c = sum(from.cell>=1); d = length(from.cell)-c;
  df = matrix(NA, nrow = length(pathway.cluster), ncol = 2); colnames(df) = c('odds', 'p.value')
  for(ci in 1:length(pathway.cluster)){
    if(length(pathway.cluster[[ci]])==0 | length(unlist(pathway.cluster[[ci]]))<=1) next
    from.cluster = unlist(pathway.cluster[[ci]]); a = sum(from.cluster>=1); b = length(from.cluster)-a; 
    ftab = matrix(c(a, b, c, d),nrow=2,dimnames=list(paired = c("yes","no"), from.acc = c("cluster","cell")));
    ftest = fisher.test(ftab, alternative = "two.sided")
    df[ci,'odds'] = round(ftest$estimate, 3); df[ci,'p.value'] = round(ftest$p.value, 100)
  }
  
  return(data.frame(df))
}


