binom.val <- function(x, y, p = 0.5, alternative = 'two.sided'){
  nvals <- length(x)
  stopifnot(nvals == length(y))
  pvals <- array(0, dim = c(nvals, 1))
  for(i in 1:nvals){
    if(y[i] > 0){
      pval = binom.test(x[i], y[i], p, alternative = alternative)
      pvals[i] <- pval$p.value
    }else{pvals[i] <- 1.0}
  }
  return(pvals)
}

binom.val.par <- function(x, y, p = 0.5, nchunks = 10, alternative = 'two.sided'){
  nvals <- length(x)
  stopifnot(nvals == length(y))
  chunk.size = ceiling(nvals / nchunks)
  pvals = foreach(i = 1:nchunks, .combine = 'append', .inorder = T) %dopar%{
    chunk.start = (i - 1) * chunk.size + 1
    chunk.end = min(i * chunk.size, nvals)
    binom.val(x[chunk.start:chunk.end], y[chunk.start:chunk.end], p, alternative = alternative)
  }
  return(pvals)
}

cor.par = function(x, y, nchunks = 10, method = 'spearman'){
  stopifnot(all(dim(x) == dim(y)))
  nvals = dim(x)[1]
  chunk.size = ceiling(nvals / nchunks)
  cor.val = foreach(i = 1:nchunks, .combine = 'append', .inorder = T) %dopar%{
    chunk.start = (i - 1) * chunk.size + 1
    chunk.end = min(i * chunk.size, nvals)
    chunk.true.size = chunk.end - chunk.start + 1
    ctmp = array(0, dim = c(chunk.true.size, 1))
    for(j in 0:(chunk.true.size - 1)){
      ctmp[j + 1] = cor(as.numeric(x[chunk.start + j, ]), as.numeric(y[chunk.start + j, ]), method = method, use = 'pairwise.complete.obs')
    }
    ctmp
  }
  return(cor.val)
}

# a = sample(1:100, 10000, replace = T)
# b = sample(101:500, 10000, replace = T)
# start = proc.time()
# foo = binom.val.par(a, b)
# print(proc.time() - start)
# 
# start = proc.time()
# bar = binom.val(a, b)
# print(proc.time() - start)