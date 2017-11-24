### RUN ROARY ###

runPangenome <- function(
  prefix,
  gffs,
  args=rargs
  # pipeline = 'roary',
  # pfam = NULL
)
{

  # pipeline <- match.arg(pipeline, c('roary','pewit'))

  # if (pipeline == 'roary'){

  p <- runRoary(gffs, args = args, prefix)

  # }else{
  #
  #   p <- pewit::pangenome(gffs, hmmPfam = )
  #
  # }


}

runRoary(gffs, args = "-p 1 -cd 95", prefix){

  if(Sys.which('roary')==''){
    stop('roary is not in your $PATH.')
  }

  out <- paste0(prefix, '_roary')
  outf <- paste0('-f ', out)
  args <- paste(args, outf)
  roary <- paste("roary", args, paste(gffs, collapse = ' '))
  system(roary)

  or <- c("accessory_binary_genes.fa", "accessory_binary_genes.fa.newick",
          "accessory_graph.dot", "accessory.header.embl", "accessory.tab",
          "blast_identity_frequency.Rtab", "clustered_proteins", "core_accessory_graph.dot",
          "core_accessory.header.embl", "core_accessory.tab", "gene_presence_absence.csv",
          "gene_presence_absence.Rtab", "number_of_conserved_genes.Rtab",
          "number_of_genes_in_pan_genome.Rtab", "number_of_new_genes.Rtab",
          "number_of_unique_genes.Rtab", "summary_statistics.txt")

  out <- paste0(out, '/', or)

  return(out)
}

