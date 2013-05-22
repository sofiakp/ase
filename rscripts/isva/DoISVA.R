DoISVA = function(data.mat, pheno, cf.m = NULL, pvthCF = 0.01, th = 0.05, ncomp = NULL, sel.col = 1:ncol(data.matat)){
  
  stopifnot(all(sel.col %in% 1:ncol(data.mat)), length(sel.col) > 0)
  isv = isvaFn(data.mat, pheno, ncomp, sel.col)
  
  if(is.null(cf.m)==FALSE){
    # TODO: THIS PART IGNORES SEL.COL
    ### study pattern of correlation of ISVA components to POI and CFs
    tmp.m <- cbind(pheno,cf.m);
    treatfactor <- c(is.factor(pheno),factor.log);
    pv.m <- matrix(nrow=ncol(isva.o$isv),ncol=1+ncol(cf.m));
    colnames(pv.m) <- c("POI",colnames(cf.m)); ## POI:phenotype of interest
    for(c in 1:ncol(tmp.m)){
      if(treatfactor[c]==FALSE){
        for(sv in 1:ncol(isva.o$isv)){
          lm.o <- lm(isva.o$isv[,sv] ~ as.numeric(tmp.m[,c]));
          pv.m[sv,c] <- summary(lm.o)$coeff[2,4];   
        }
      }
      else {
        for(sv in 1:ncol(isva.o$isv)){
          lm.o <- lm(isva.o$isv[,sv] ~ as.factor(tmp.m[,c]));
          pv.m[sv,c] <- pf(summary(lm.o)$fstat[1],summary(lm.o)$fstat[2],summary(lm.o)$fstat[3],lower.tail=FALSE);   
        }
      }
    }
    
    ### selection of ISVs
    print("Selecting ISVs");
    selisv.idx <- vector();
    for(sv in 1:nrow(pv.m)){
      
      ncf <- length(which(pv.m[sv,2:ncol(pv.m)]< pvthCF)) ## pvth=0.01
      minpv <- min(pv.m[sv,2:ncol(pv.m)]);
      phpv <- pv.m[sv,1];
      if(ncf > 0){
        if(minpv < phpv){
          selisv.idx <- c(selisv.idx,sv);
        }
      }
    }
    if (length(selisv.idx)==0 ){
      print("No ISVs selected because none correlated with the given confounders. Rerun ISVA with cf.m=NULL option"); stop;
    }
    
  }
  else { ### confounder matrix not given, so select all ISVs
    selisv.idx = 1:ncol(isv)
    pv.m = NULL
  }
  
  selisv.m = matrix(isv[, selisv.idx], ncol = length(selisv.idx))
  print("Running final multivariate regressions with selected ISVs")
  mod.sel = model.matrix( ~ pheno[sel.col] + selisv.m[sel.col, ])
  modNULL.sel = model.matrix( ~ selisv.m[sel.col,])
  
  df1 = dim(mod.sel)[2]
  df0 = dim(modNULL.sel)[2]
  
  theta = data.mat[, sel.col] %*% mod.sel %*% solve(t(mod.sel) %*% mod.sel)
  resid = data.mat[, sel.col] - theta %*% t(mod.sel)
  rss1 = rowSums(resid * resid)
  thetaNULL = data.mat[, sel.col] %*% modNULL.sel %*% solve(t(modNULL.sel) %*% modNULL.sel)
  residNULL = data.mat[, sel.col] - thetaNULL %*% t(modNULL.sel)
  rssNULL = rowSums(residNULL * residNULL)
  
  # Get F-statistics for comparing the null and full model
  fstats = ((rssNULL - rss1)/(df1 - df0)) / (rss1/(length(sel.col) - df1))
  pv.v = 1 - pf(fstats, df1 = (df1 - df0), df2 = (length(sel.col) - df1))
  pv.s = sort(pv.v, decreasing = FALSE, index.return = TRUE)
  qv.v = p.adjust(pv.s$x, method = 'BH') #qvalue(pv.s$x)$qvalue # TODO: change this to p.adjust
  ntop = sum(qv.v < th)
  cat('Number of DEGs after ISV adjustment =', ntop, '\n')
  
  if(ntop > 0){
    pred.idx = pv.s$ix[qv.v < th];
    if(ntop == 1){
      lm.o = lm(data.mat[pred.idx, sel.col] ~ pheno[sel.col] + selisv.m[sel.col, ] )
      tstats.v = summary(lm.o)
      tstats.v = tstats.v$coeff[2,3]
    }else{
      # find t-stats of significant ones
      lm.o = lm(t(data.mat[pred.idx, sel.col]) ~ pheno[sel.col] + selisv.m[sel.col, ]) 
      tstats.v = unlist(lapply(summary(lm.o), function(x){ x$coeff[2,3];}))
    }
    lm.m = cbind(tstats.v, pv.s$x[qv.v < th], qv.v[qv.v < th])
    colnames(lm.m) = c("t", "pval", "qval")
  } 
  else {
    pred.idx <- NULL;
    lm.m <- NULL;
  }
  
  if(length(sel.col) < ncol(data.mat)){
    if(ncol(data.mat) - length(sel.col) == 1){
      residNULL.other = data.mat[, -sel.col] - thetaNULL %*% t(model.matrix( ~ t(selisv.m[-sel.col,])))
      resid.other = data.mat[, -sel.col] - theta %*% t(model.matrix( ~ t(pheno[-sel.col]) + t(selisv.m[-sel.col,])))
    }else{
      residNULL.other = data.mat[, -sel.col] - thetaNULL %*% t(model.matrix( ~ selisv.m[-sel.col,]))
      resid.other = data.mat[, -sel.col] - theta %*% t(model.matrix( ~ pheno[-sel.col] + selisv.m[-sel.col,]))
    }
    residNULL.tmp = array(0, dim = c(nrow(residNULL), ncol(data.mat)))
    residNULL.tmp[, sel.col] = residNULL
    residNULL.tmp[, -sel.col] = residNULL.other
    residNULL = residNULL.tmp
    
    resid.tmp = array(0, dim = c(nrow(resid), ncol(data.mat)))
    resid.tmp[, sel.col] = resid
    resid.tmp[, -sel.col] = resid.other
    resid = resid.tmp
  }
  return(list(pval = pv.s$x, qval = qv.v, rank = pv.s$ix, ndeg = ntop, deg = pred.idx, isv = selisv.m, res.null = residNULL, res = resid))
  
}
