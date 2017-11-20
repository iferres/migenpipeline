### RUN PROGRESSIVEmAUVE ###

runProgressiveMauve <- function(gffs, args='', prefix)
{

  if(Sys.which('progressiveMauve')==''){
    stop('progressiveMauve is not in your $PATH.')
  }

  fastas <- vapply(gffs, extractFastaFromGff, FUN.VALUE = '')

  out <- paste0('--output=',prefix,'_progressiveMauve')
  args <- paste(args, out)
  pMauve <- paste('progressiveMauve', args, fastas)
  system(pMauve)

  file.remove(fastas)

  return(out)
}


extractFastaFromGff <- function(gff)
{

  rl <- readLines(gff)
  st <- grep('##FASTA', rl, fixed = TRUE) + 1L
  en <- length(rl)

  out <- sub('gff$','fasta', gff)

  writeLines(rl[st:en], con = out)

  return(out)
}
