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


#' `xmfa` is the progressiveMauve output alignment.
#' `gffs` is a vector with the gff3 files.
#' `pco` and `pal` should be taken together. Only if a LCB has less than
#' `pal` columns (proportion) with a proportion grater than `pco` of gaps, the
#' LCB is then considered in the final alignment. This is suposed to be more
#' 'correct' than passing gblocks or trimAl directly to the progressiveMauve
#' alignment because less artificial junctions are generated. The drawback is
#' that the final alignment may be shorter.
#' `msi` is the minimum LCB size to be considered, and `nco` is the number
#' of conserved sequences a LCB has to have to be considered. For now, nco
#' must be the number of organisms in the alignment, other number is not
#' implemented yet.
#' The output is a FASTA alignment.
stripSubsetLCBs <- function(xmfa,
                            gffs,
                            msi=500L,
                            pco=0.1, #Less is less gaps
                            pal=0.9, #More is less gaps
                            nco=length(gffs)){

  ncom <- (length(gffs) * 2L) + 2L
  rl <- readLines(xmfa)
  rl <- rl[-(1L:ncom)]

  #index of chunks (blocks) of sequences in xmfa
  eq <- grep('^=', rl)
  vp <- vapply(1:length(eq), function(x){

    ini <- ifelse(eq[x]==eq[1], 1L, eq[x-1L]+1L)
    end <- eq[x]-1L
    c(ini, end)

  }, FUN.VALUE = c(1L,1L))

  # Stats about chunks
  ap <- t(apply(vp, 2, function(x){

    rr <- rl[x[1]:x[2]]
    gp <- grep('^>', rr)
    ln <- length(gp)
    nch <- nchar(paste0(rr[2L:(ifelse(ln>1L, gp[2L]-1, length(rr)-1L))],
                        collapse = ''))

    vs <- vapply(1:length(gp), function(y){

      ini <- gp[y] + 1L
      end <- ifelse(gp[y]!=rev(gp)[1], gp[y+1L] - 1L, length(rr))
      c(ini, end)

    }, FUN.VALUE = c(1L,1L))

    cb <- strsplit(apply(vs, 2, function(y){
      paste0(rr[y[1]:y[2]], collapse = '')
    }), '')

    cb <- do.call(rbind, cb)

    #Number of columns with a proportion of gaps lesser than pco.
    prco <- length(which(apply(cb, 2, function(z){
      #proportion of gaps in each column.
      length(which(z=='-'))/ln
      })<=pco))

    #Proportion of columns with less than pco gaps.
    pral <- prco/nch

    c(ln, nch, pral)

  }))

  #Save blocks which have a proportion of columns with less than pco gaps more
  #than pal.
  wh <- which(ap[,1]==nco & ap[,2]>=msi & ap[,3]>=pal)

  if (length(wh)==0L){
    stop('No LCBs pass filters.')
  }

  #Extract selected blocks for each genome.
  ck <- t(vp[, wh])
  ff <- apply(ck, 1, function(x){

    rr <- rl[x[1]:x[2]]
    gp <- grep('^>', rr)

    chu <- sub('[.]xmfa$','.lcbs',xmfa)
    cat(rr[gp], sep = '\n', file = chu, append = TRUE)
    cat('=',sep = '\n', file = chu, append = TRUE)

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
  rm(ff)

  #Create outfile
  out <- sub('xmfa','fasta',xmfa)
  file.create(out)

  #Concatenate blocks
  he <- sub('[.]gff$', '', gffs[as.integer(sub('[.]core_fasta','',tmps))])
  vapply(1:length(tmps), function(y){
    cat(paste0('>',he[y], '\n'), file = out, append = TRUE)
    file.append(out, tmps[y])
    cat('\n', file = out, append = TRUE)
    TRUE
  }, FUN.VALUE = TRUE)

  file.remove(tmps)

  return(out)
}




