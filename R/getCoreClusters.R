getCoreClusters <- function(lcbs,
                            gffs,
                            roary_clusters,
                            prefix,
                            mc.cores = 1L){

  #Read gffs as list with data.frame and gene sequences
  x <- parallel::mclapply(gffs, readGff, mc.cores = mc.cores)
  names(x) <- sub('[.]gff$','', sapply(strsplit(gffs,'/'), function(y){
    rev(y)[1]
    }))

  #Read LCBs
  lcb <- readLCBS(lcbs)

  #Read roary's clusters
  clusters <- readRoaryClusters(roary_clusters)

  #Get which genes appear in the progressiveMauve core-alignment.
  coreGenes <- parallel::mclapply(names(x), function(i){
    x[[i]][[1]][which(apply(x[[i]][[1]], 1, function(y){
      any(vapply(lcb, function(z){
        overlap(a = c(y[[7]],y[[8]]), b = c(z[i,2], z[i,3]))
      }, FUN.VALUE = TRUE))
    })),2]
  }, mc.cores = mc.cores)
  names(coreGenes) <- names(x)

  #Get which roary's clusters appear in the progressiveMauve core-alignment.
  ul <- unlist(coreGenes)
  ulinx <- parallel::mclapply(clusters, function(y){
    any(ul%in%y)
  }, mc.cores = mc.cores)
  uliny <- unlist(uliny)

  clu <- clusters[names(which(uliny))]

  #Just keep clusters with at least one gene per genome.
  clu <- clu[-which(sapply(clu, length)>length(gffs))]

  #Just keep one sequence per genome for each cluster


}

readGff <- function(gff){
  rl <- readLines(gff)
  x <- getGffTable(rl)
  s <- getGffSeqs(x, rl)
  list(x,s)
}

getGffTable <- function(rl){

  # rl <- readLines(gff)

  w <-which(grepl('^\\#\\#',rl))

  re <- rl[w][-c(1, length(w))]
  re <- do.call(rbind,lapply(strsplit(re,' '), '[', c(2:4)))
  re <- as.data.frame(re,stringsAsFactors=FALSE)
  re[,2] <- as.integer(re[,2])
  re[,3] <- as.integer(re[,3])

  re$V5 <- cumsum(re$V3)
  re$V4 <- re$V5 - re$V3 + 1L
  re <- re[, c(1,2,3,5,4)]

  upto <- rev(w)[1] - 1
  from <- rev(w)[2] + 1
  o <- rl[from:upto]

  lst <- strsplit(o,'\t')

  contig <- sapply(lst, function(x){ x[1] })
  type <- sapply(lst, function(x){ x[3] })
  from <- as.integer(sapply(lst, function(x){ x[4] }))
  to <- as.integer(sapply(lst, function(x){ x[5] }))
  strand <- sapply(lst, function(x){ x[7] })
  phase <- sapply(lst, function(x){ x[8] })
  attrib <- sapply(lst, function(x){ x[9] })

  metadata <- strsplit(attrib,';')

  id <- sapply(metadata,function(x){
    gp<-grep('ID=',x,value = T)
    if(length(gp)>0){sub('ID=','',gp)}else{''}
  })
  locustag <- sapply(metadata,function(x){
    gp<-grep('locus_tag=',x,value = T)
    if(length(gp)>0){sub('locus_tag=','',gp)}else{''}
  })
  gene <- sapply(metadata,function(x){
    gp<-grep('gene=',x,value = T)
    if(length(gp)>0){sub('gene=','',gp)}else{''}
  })
  product <- sapply(metadata,function(x){
    gp<-grep('product=',x,value = T)
    if(length(gp)>0){sub('product=','',gp)}else{''}
  })

  out <- data.frame(Contig=contig,
                    ID=id,
                    LocusTag=locustag,
                    Gene=gene,
                    Product=product,
                    Type=type,
                    From=from,
                    To=to,
                    Strand=strand,
                    Phase=phase,
                    stringsAsFactors = F)

  out$From <- apply(out, 1, function(x){
    re[which(re[, 1]==x[1]), 4] + as.integer(x[7]) - 1L
  })

  out$To <- apply(out, 1, function(x){
    re[which(re[ ,1]==x[1]), 4] + as.integer(x[8]) - 1L
  })

  out

}


getGffSeqs <- function(x, rl){
  gp <- grep('##FASTA',rl, fixed = TRUE) + 1L
  fna <- rl[gp:length(rl)]


  gph <- grep('^>', fna)
  fna <- paste0(fna[-gph], collapse = '')
  fna <- strsplit(fna, '')[[1]]

  sqs <- apply(x, 1, function(y){
    gen <- fna[y[[7]]:y[[8]]]
    if(y[[9]]=='-'){
      gen <- rev(seqinr::comp(gen, forceToLower = FALSE))
    }

    se <- paste0(gen, collapse = '')
    seqinr::as.SeqFastadna(se)
  })

  sqs <- as.list(sqs)
  names(sqs) <- x$LocusTag

  return(sqs)
}

readLCBS <- function(lcbs){

  lcb <- readLines(lcbs)
  eq <- grep('^=', lcb)
  vp <- vapply(1:length(eq), function(x){

    ini <- ifelse(eq[x]==eq[1], 1L, eq[x-1L]+1L)
    end <- eq[x]-1L
    c(ini, end)

  }, FUN.VALUE = c(1L,1L))

  apply(vp, 2, function(x){
    sp <- strsplit(lcb[x[1]:x[2]], ':')
    sp2 <- lapply(sp, function(y){strsplit(y,' ')[[2]]})

    df <- vector('list', 5)
    names(df) <- c('fasta', 'start', 'end', 'strand', 'length')
    df$fasta <- sapply(sp2, '[', 3)
    df$start <- as.integer(sapply(sp2, function(y) {
      strsplit(y[1],'-')[[1]][1]
      } ))
    df$end <- as.integer(sapply(sp2, function(y) {
      strsplit(y[1],'-')[[1]][2]
      } ))
    df$strand <- sapply(sp2, '[', 2)
    df$length <- as.integer(sapply(sp, function(y){
      sub(' nch = ','', y[3])
      }))

    df <- as.data.frame(df)
    rownames(df) <- sub('[.]fasta$','', df$fasta)

    return(df)
  })


}


readRoaryClusters <- function(roary_clusters){

  cl <- readLines(roary_clusters)
  sp <- strsplit(cl, ': ')

  out <- lapply(sp, function(x){
    strsplit(x[2],'\t')[[1]]
  })

  names(out) <- sapply(sp, '[', 1)

  return(out)
}



overlap <- function(a=c(1L, 1L), b=c(1L, 1L)){
  if(a[1]>a[2]){
    a <- rev(a)
  }
  if(b[1]>b[2]){
    b <- rev(b)
  }
  return(a[1]<=b[2] & a[2]>=b[1])
}






