get_targets = function(datalist, clusters, clen){
  
  sequences = lapply(datalist, function(x) rbind(x$train[x$train$label==1,], x$test[x$test$label==1,]))
  scoretab = sapply(1:clen, function(x){ 
    if(class(clusters[[x]])!='gbm'){
      rep(NA, sum(sapply(sequences, nrow)))
    }else{
      unlist(sapply(cells, function(y){
        if(sum(grepl(y, names(clusters[[x]]$trees)))==0){
          return(rep(NA, nrow(sequences[[y]])))
        }else{
          return(get_modelscore(clusters[[x]], sequences[[y]], clusters[[x]]$n.trees))
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

setup_exp_coherence = function(targets, exprsn, exptheK=1, verbose=T){
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

exp_coherence = function(coherence){#do i need target
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

setup_pathway_coherence = function(genemat, pathwaymat, verbose=T){
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

pathway_coherence = function(coherence){
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
