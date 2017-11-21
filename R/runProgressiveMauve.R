### RUN PROGRESSIVEmAUVE ###

runProgressiveMauve <- function(gffs, args='', prefix)
{

  if(Sys.which('progressiveMauve')==''){
    stop('progressiveMauve is not in your $PATH.')
  }

  fastas <- vapply(gffs, extractFastaFromGff, FUN.VALUE = '')

  out <- paste0(prefix, '_progressiveMauve.xmfa')
  oarg <- paste0('--output=',out)
  args <- paste(args, oarg)
  pMauve <- paste('progressiveMauve', args, paste0(fastas, collapse = ' '))
  system(pMauve)

  file.remove(fastas)
  file.remove(paste0(fastas, '.sslist'))

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



stripSubsetLCBs <- function(xmfa, gffs, msi=500L, nco=length(gffs)){

  ncom <- (length(gffs) * 2L) + 2L
  rl <- readLines(xmfa)
  rl <- rl[-(1L:ncom)]

  eq <- grep('^=', rl)
  vp <- vapply(1:length(eq), function(x){

    ini <- ifelse(eq[x]==eq[1], 1L, eq[x-1L]+1L)
    end <- eq[x]-1L
    c(ini, end)

  }, FUN.VALUE = c(1L,1L))

  ap <- t(apply(vp, 2, function(x){

    rr <- rl[x[1]:x[2]]
    gp <- grep('^>', rr)
    ln <- length(gp)
    nch <- nchar(paste0(rr[2L:(ifelse(ln>1L, gp[2L]-1, length(rr)-1L))],
                        collapse = ''))
    c(ln, nch)

  }))

  wh <- which(ap[,1]==nco & ap[,2]>=msi)

  if (length(wh)==0L){
    stop('No LCBs pass filters.')
  }

  ck <- t(vp[, wh])
  ff <- apply(ck, 1, function(x){

    rr <- rl[x[1]:x[2]]
    gp <- grep('^>', rr)
    idx <- vapply(strsplit(sub('> ','',rr[gp]),':'), '[', 1, FUN.VALUE = '')

    ini <- gp+1L
    end <- vapply(1:length(gp), function(y){
      ifelse(gp[y]!=rev(gp)[1], gp[y+1L]-1L, length(rr))
      }, FUN.VALUE = 1L)

    m <- cbind(ini, end)
    fnam <- paste0(idx, '.core_fasta')

    fis <- vapply(1:nrow(m), function(y){
      cat(rr[m[y,1]:m[y,2]], file = fnam[y], sep = '', append = TRUE)
      fnam[y]
    }, FUN.VALUE = '')

    fis

  })

  tmps <- ff[,1]

  out <- sub('xmfa','fasta',xmfa)

  file.create(out)
  he <- sub('[.]gff$', '', gffs[as.integer(sub('[.]core_fasta','',tmps))])

  vapply(1:length(tmps), function(y){
    cat(paste0('>',he[y], '\n'), file = out, append = TRUE)
    file.append(out, tmps[y])
  }, FUN.VALUE = TRUE)

  file.remove(tmps)

  return(out)
}




