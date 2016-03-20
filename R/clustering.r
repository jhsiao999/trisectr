load_weight_matrix = function(models){
  
  cells = names(models)
  col.cells = unlist(sapply(cells, function(x) rep(x, models[[x]]$n.trees)))
  #oobag.improve = sapply(models, function(x) x$oobag.improve)
  mweights = sapply(models, function(x) sapply(1:x$n.trees, 
                                               function(y) feature_importance(x, n.trees=y, single.tree=T)$rel.inf))
  
  mweights = t(Reduce(cbind, mweights)); colnames(mweights) = models[[cells[1]]]$var.names
  fvals = apply(mweights, 2, sum); mweights = mweights[,colnames(mweights)[fvals>0]]
  return(list(mweights=mweights, col.cells=col.cells))#, oobag.improve=oobag.improve))
}

get_cluster_membership = function(models, clen, ctype='knn', plotit=F){
  
  ld = load_weight_matrix(models); mweights = ld$mweights
  if(ctype=='knn'){
    fit = kmeans(mweights, clen) 
    cluster.membership = as.data.frame.matrix(table(ld$col.cells, fit$cluster))
  }else{
    fit = xyf(mweights, classvec2classmat(ld$col.cells), xweight=xweight,
              somgrid(4, clen/4, "hexagonal"), rlen=500)
    cluster.membership = as.data.frame.matrix(table(ld$col.cells, fit$unit.classif))
  }
  
  if(plotit==T){
    plot_cluster_membership(cluster.membership, desc)
  }
  return(list(fit=fit, cluster.membership=cluster.membership))
}

plot_cluster_membership = function(cluster.membership, desc){
  
  colbreaks = c(seq(0,0.2,length=100), seq(0.2001,1,length=100), 
                seq(1.0001,max(cluster.membership, na.rm=T),length=100))
  pallete = colorRampPalette(brewer.pal(9,"YlGn")[c(1,5,8)])(n = 299)
  #pdf(paste0(getwd(), '/cluster.membership.pdf'))
  heatmap.2(t(as.matrix(cluster.membership)),  main = desc, 
            cellnote=t(as.matrix(cluster.membership)), cexRow=1, notecol="black", 
            density.info="none", trace="none", margins =c(6,7), col=pallete,  dendrogram='col',        
            cexCol=1.2,  Rowv=F, breaks=colbreaks, keysize=0.8)
  #dev.off()
}

make_cluster_ensembles = function(models, fit, clen, ctype='knn'){
  
  ld = load_weight_matrix(models); mweights = ld$mweights
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
