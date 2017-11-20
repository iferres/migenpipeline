micropipeline <- function(gffs,
                          map, # Mapping file
                          rargs = '-p 1 -cd 95',
                          pargs = '',
                          gargs = ''
                          # pgpl = 'roary', #pangenome pipeline: roary or pewit
                          # pfam = NULL, # c('pfamA.hmm', 'pfamA.hmm.dat')
                          # genes = NULL # fasta files
)
{

  deps <- checkDependencies()

  if (missing(prefix)){
    prefix <- format(Sys.time(), "%b%d%H%M%S")
  }

  #Run main programs
  # (Capture outputs in lists)

  pang <- runPangenome(prefix, gffs) # o pewit

  progMauve <- runProgressiveMauve(gffs)

  pMauve2

  gubbins <- runGubbins(progMauve)

  pStats <- pangStats(pang, progMauve, gubbins)

  # recSeqs <- extractRecSeqs(gffs, gubbins)

  # domContent <- runDCDA(gffs, pfam, map) # Produce PCoA plot

  # phylo <- makePhylo(coreAl = c("progMauve", "pang"), gubbins) # mask recomb, phangorn

  # genes <- findGenes(pang, genes, map) #Plot ggheatmap-like plot


}
