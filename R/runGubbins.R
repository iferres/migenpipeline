runGubbins <- function(fastaCore, prefix, gargs=''){

  if (Sys.which('run_gubbins.py')==''){
    stop('"run_gubbins.py" in not in your $PATH.')
  }

  gargs <- paste0(gargs,' -p ', prefix)

  gubbins <- paste('run_gubbins.py', gargs, fastaCore)
  system(gubbins)

}
