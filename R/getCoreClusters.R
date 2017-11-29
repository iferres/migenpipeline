getCoreClusters <- function(lcbs, gffs, roary_clusters, prefix){

  x <- lapply(gffs, readGffTable)
  names(x) <- sapply(strsplit(gffs,'/'), function(y){rev(y)[1]})

  lcb <- readLCBS(lcbs)

  clusters <- readRoaryClusters(roary_clusters)



}


readGffTable <- function(gff){

  rl <- readLines(gff)

  w <-which(grepl('^\\#\\#',rl))
  upto <- rev(w)[1] - 1
  from <- rev(w)[2] + 1
  o <- rl[from:upto]

  lst <- strsplit(o,'\t')

  contig <- sapply(lst, function(x){ x[1] })
  type <- sapply(lst, function(x){ x[3] })
  from <- sapply(lst, function(x){ x[4] })
  to <- sapply(lst, function(x){ x[5] })
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
  out

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

    as.data.frame(df)
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
  a <- sort(a)
  b <- sort(b)

  a[1]<=b[2] & a[2]>=b[1]
}



apply(x[[1]], 1, function(y){
  vapply(lcb, function(z){
    overlap(a = c(y[[7]],y[[8]]), b = c(z[1,2], z[1,3]))
  },FUN.VALUE = TRUE)
})
