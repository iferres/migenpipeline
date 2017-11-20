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

  out <- paste0('-f ', prefix, '_roary')
  args <- paste(args, out)
  roary <- paste("roary", args, paste(gffs, collapse = ' '))
  system(roary)

  return(out)
}

