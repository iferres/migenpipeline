micropipeline <- function(fnas,
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

  gffs <- runProkka(fnas)

  pang <- runPangenome(prefix, gffs, rargs) # o pewit

  xmfa <- runProgressiveMauve(gffs, pargs)

  fastaCore <- stripSubsetLCBs(xmfa = xmfa, gffs, msi=500L, nco=length(gffs))

  gubbins <- runGubbins(fastaCore, gargs)

  pStats <- pangStats(pang[11], progMauve, gubbins[2])

  # recSeqs <- extractRecSeqs(gffs, gubbins)

  # domContent <- runDCDA(gffs, pfam, map) # Produce PCoA plot

  # phylo <- makePhylo(coreAl = c("progMauve", "pang"), gubbins) # mask recomb, phangorn

  # genes <- findGenes(pang, genes, map) #Plot ggheatmap-like plot


}
