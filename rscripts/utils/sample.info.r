sample.info <- function(x, suf = '\\.[a-zA-z0-9]*$'){
  sample <- gsub('INPUT$', 'INPUT_1', gsub(paste('SNYDER_HG19_|', suf, sep = ''), '', basename(x), perl = T))
  fields <- unlist(strsplit(sample, '_'))
  if(length(fields) < 3){
    return(NULL)
  }
  indiv <- fields[seq(1, length(fields), 3)]
  mark <- fields[seq(2, length(fields), 3)]
  rep <- fields[seq(3, length(fields), 3)]
  return(data.frame(indiv = indiv, mark = mark, rep = rep))
}