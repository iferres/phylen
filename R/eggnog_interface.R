#' @name list_eggnogdb
#' @title List available EggNOG db datasets
#' @description This function list all the available EggNOG datasets to
#' download, associated with their nog prefix.
#' @return A \code{data.frame} with the available datasets and their respective
#' nog prefixes.
#' @author Ignacio Ferres
#' @importFrom httr GET
#' @importFrom XML getHTMLLinks
#' @importFrom utils read.csv
#' @export
list_eggnogdb <- function(){

  b <- 'http://eggnogdb.embl.de/download/latest/'
  ln <- getHTMLLinks(b)

  tx <- grep('_species_by_taxonomic_level.tsv',
             ln,
             fixed = TRUE,
             value = TRUE)
  ta <- paste0(b, tx, collapse = '')
  tsv <- GET(ta)

  if(httr::status_code(tsv)==200){

    prs <- unlist(strsplit(rawToChar(tsv$content), '\n'))
    res <- read.csv(text=prs, sep='\t')[, c(3,2)]
    res

  }else{
    stop('Something went wrong with the conection.')
  }

}


#' @name download_nog_hmm
#' @title Download specific hmm datasets from EggNOG.
#' @description Takes a nog prefix and downloads a compressed (.tar.gz) file
#' from EggNOG db containing the hmm models of the specified group.
#' @param nog.prefix \code{character}, a nog prefix as returned by
#' \link{list_eggnogdb} function.
#' @param onDir The directory name where to download the file.
#' @return A tar.gz compressed file on the specified directory. On console it
#' returns the path to the downloaded file.
#' @author Ignacio Ferres
#' @importFrom utils download.file
#' @export
download_nog_hmm <- function(nog.prefix='proNOG', onDir='.'){

  onDir <- normalizePath(dirname(onDir))
  b <- 'http://eggnogdb.embl.de/download/latest/data/'
  hmm <- paste0(nog.prefix, '.hmm.tar.gz', collapse = '')
  dwn <- paste0(b,
                nog.prefix,
                '/',
                hmm,
                collapse = '')

  destfile <- paste0(onDir, '/', hmm, collapse = '')
  download.file(url = dwn, destfile = hmm)
  destfile

}
