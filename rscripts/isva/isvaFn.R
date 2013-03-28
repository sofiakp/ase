require(fastICA)
require(qvalue)

# Runs linear regression of data.mat on pheno and then performs ICA on the residual matrix.
# ncomp is the number of ICA components. Set to NULL for an automatic selection.
# sel.col is a vector of the columns of data.mat that will be used for ICA. The rest
# of the columns will just be projected.
isvaFn = function(data.mat, pheno, ncomp = NULL, sel.col = 1:ncol(data.mat)){

  lm.o = lm(t(data.mat[, sel.col]) ~ pheno[sel.col]) # Fit a linear regression for each gene separately
  res.m = t(lm.o$residuals) 
  
  if(is.null(ncomp)){
    rmt.o <-  EstDimRMT(res.m)
    ncomp <- rmt.o$dim;
    print(paste("Number of candidate ISVs = ",ncomp,sep=""));
  }
  else {
    print("no need to estimate dimensionality");
  }

  # ICA on the residual matrix
  fICA.o = fastICA(res.m, n.comp = ncomp, maxit = 1000)
  isv.tmp = t(fICA.o$A)
  isv = array(0, dim = c(ncol(data.mat), ncomp))
  isv[sel.col, ] = isv.tmp
  if(length(sel.col) < ncol(data.mat)){
    # Residuals for the columns that were not used for fitting
    res.other = data.mat[, -sel.col] - t(model.matrix(~ pheno[-sel.col]) %*% lm.o$coefficients)
    # solve data.mat = SA for the given mixing matrix S 
    isv[-sel.col, ] = as.numeric(t(qr.solve(fICA.o$S, res.other))) 
  }
  return(isv)
  
#   isv.m <- tmp.m;
#   sd <- 1/sqrt(ncol(data.mat)-3);
#   for(k in 1:ncol(tmp.m)){
#    cor.v <- as.vector(cor(t(data.mat),tmp.m[, k]))
#    z.v <- 0.5*log((1+cor.v)/(1-cor.v));
#    pv.v <- 2*pnorm(abs(z.v),0,sd,lower.tail=FALSE)
#    tmp.s <- sort(pv.v,decreasing=FALSE,index.return=TRUE);
#    qv.o <- qvalue(pv.v);
#    nsig <- length(which(qv.o$qvalues<0.01));
#    nsig = max(nsig, 500)
#    red.m <- data.mat[tmp.s$ix[1:nsig],];
#    fICA.o <- fastICA(red.m,n.comp=ncomp);
#    cor.v <- abs(cor(tmp.m[,k],t(fICA.o$A)));
#    kmax <- which.max(cor.v);
#    isv.m[,k] <- t(fICA.o$A)[,kmax];
#    print(paste("Built ISV ",k,sep=""));   
#   }
#   return(list(n.isv=ncomp,isv=isv.m));
}

