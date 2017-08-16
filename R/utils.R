#' @name untargz
#' @title Decompress A Tar.gz File
#' @description Decompress a \code{.tar.gz} file.
#' @param targzfile The full path to the file to decompress.
#' @param exdir The directory to put the decompressed files.
#' @return A \code{character} vector with the decompressed files.
#' @author Ignacio Ferres
untargz <- function(targzfile, exdir){

  if (missing(exdir)){
    exdir <- dirname(targzfile)
    exdir <- normalizePath(exdir)
  }
  if(!dir.exists(exdir)){
    dir.create(exdir)
  }

  unt <- paste0('tar -xzvf ',targzfile,' -C ',exdir)
  fils <- system(unt, intern = TRUE)
  fils <- paste0(exdir, '/', fils)

  return(fils)
}

#' @name setHMM
#' @title Prepare Compressed Hmm Files To Hmmsearch
#' @description Decompress, concatenate and press (\code{hmmpress}) hmm files
#' to allow \code{hmmsearch} to run.
#' @param hmmTarGz The path to the \code{.hmm.tar.gz} file downloaded from
#' EggNOG website.
#' @return A \code{character} vector with the names of the hmm file and their
#' indices.
#' @author Ignacio Ferres
setHmm <- function(hmmTarGz){

  cat('Decompressing.. ')
  hmms <- untargz(targzfile = hmmTarGz)
  cat('DONE!\n')

  cat('Concatenating.. ')
  nam <- sub('.tar.gz', '', rev(strsplit(hmmTarGz, '/')[[1]])[1], fixed = T)
  cate <- paste0('cat ', paste(hmms[1:5], collapse = ' '), ' > ',nam)
  system(cate)
  file.remove(hmms)
  cat('DONE!\n')

  cat('Pressing.. ')
  hmm <- hmmPress(model = nam)
  cat('DONE!\n')

  return(hmm)
}
