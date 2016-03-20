
get_modelscore = function(submodels, test, best.iter, vtype='response'){
  submodels.pr <- predict(submodels, test[-1], best.iter, type=vtype);
  
  return(submodels.pr)
}

get_modelsize = function(model){
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

get_performance = function(submodels.pred, digits=3, verbose=T){
  auc = performance(submodels.pred,"auc")@y.values[[1]]
  if(verbose==T){
    print(auc) 
  }
  return(round(auc, digits=digits))
}

build_emt = function(dataset, distb='huberized', etype='interaction', interaction.depth=15, 
                         bag.fraction=0.5, shrinkage=0.05, n.minobsinnode=30, n.trees = 100, 
                         ...){
  
  if(etype=='noninteraction'){
    
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
  
  return(list(model=submodels, auc=auc))
}

test_emt = function(submodes, testdata){
  modelscores = get.modelscore(submodels, testdata, submodels$n.trees); 
  model.pred = prediction(modelscores, testdata$label)
  auc = get.performance(model.pred)
  
  return(auc)
}

feature_importance = function (object, n.trees=object$n.trees, single.tree=F, normalize=T) {
  if (n.trees < 1) {
    stop("n.trees must be greater than 0.")
  }
  if (n.trees > object$n.trees) {
    warning("Exceeded total number of GBM terms. Results use n.trees=", object$n.trees, " terms.\n")
    n.trees <- object$n.trees
  }
  rel.inf <- relative_inf(object, n.trees, single.tree=single.tree) #modified for single tree
  rel.inf[rel.inf < 0] <- 0
  if (normalize) 
    rel.inf <- 100 * rel.inf/sum(rel.inf)
  
  return(data.frame(var = object$var.names, rel.inf = rel.inf)) 
}

relative_inf = function (object, n.trees, single.tree=F) {
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



