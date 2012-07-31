binom.val <- function(x, y, p = 0.5){
  nvals <- length(x)
  stopifnot(nvals == length(y))
  pvals <- array(0, dim = c(nvals, 1))
  for(i in 1:nvals){
    if(y[i] > 0){
      pval = binom.test(x[i], y[i], p, alternative = 'two.sided')
      pvals[i] <- pval$p.value
    }else{pvals[i] <- 1.0}
  }
  return(pvals)
}