runGubbins <- function(fastaCore, prefix, gargs='--threads 4 -t fasttree'){

  if (Sys.which('run_gubbins.py')==''){
    stop('"run_gubbins.py" in not in your $PATH.')
  }

  gargs <- paste0(gargs,' -p ', prefix)

  gubbins <- paste('run_gubbins.py', gargs, fastaCore)
  system(gubbins)

  out <- c()
  out[1] <- paste0(prefix, '.branch_base_reconstruction.embl')
  out[2] <- paste0(prefix, '.recombination_predictions.gff')
  out[3] <- paste0(prefix, '.recombination_predictions.embl')
  out[4] <- paste0(prefix, '.filtered_polymorphic_sites.phylip')
  out[5] <- paste0(prefix, '.summary_of_snp_distribution.vcf')
  out[6] <- paste0(prefix, '.per_branch_statistics.csv')
  out[7] <- paste0(prefix, '.filtered_polymorphic_sites.fasta')
  out[8] <- paste0(prefix, '.node_labelled.final_tree.tre')
  out[9] <- paste0(prefix, '.final_tree.tre')

}
