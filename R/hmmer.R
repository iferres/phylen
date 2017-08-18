#' @name hmmStat
#' @title Compute Hmm Stats
#' @description Compute hmm stats.
#' @param hmmfile \code{character} The path to the hmm file.
#' @return The file name of the temporary file where stats were written.
hmmStat <- function(hmmfile){
  tmp <- tempfile()
  hmmstat <- paste0('hmmstat ', hmmfile, ' > ', tmp)
  system(hmmstat)
  return(tmp)
}

#' @name hmmPress
#' @title hmmPress
#' @description Wrapper function of \code{hmmpress} (HMMER 3).
#' @param model \code{character} The name of the hmm file.
#' @return The names of indexed files.
hmmPress <- function(model){
  hmmpress <- paste('hmmpress -f', model)
  system(hmmpress, ignore.stdout = TRUE)
  o <- paste0(model, c('','.h3f', '.h3i', '.h3m', '.h3p'))
  o
}

#' @name hmmSearch
#' @title Run hmmsearch (HMMER 3)
#' @description Takes a fasta file and a Hidden Markov Model profile and
#' performs a search of the former over the latter.
#' @param fasta A protein fasta file.
#' @param hmm A hmm file. Must be pressed (see hmmpress from HMMER manual).
#' @param eval evalue threshold.
#' @param oty The \code{hmmsearch} output type.
#' @param n_threads An \code{integer}. The number of cores to use.
#' @return The path to a temporary file where the hmmsearch output is placed.
hmmSearch <- function(fasta,
                      hmm,
                      eval = '1e-10',
                      oty = 'domtblout',
                      n_threads = 1L){

  oty <- match.arg(oty, c('tblout', 'domtblout'))

  #run hmmsearch
  blout <- tempfile(pattern = 'tmpo', fileext = '.tab')
  hmmse <- paste0('hmmsearch -o /dev/null --noali -E',
                  eval,
                  paste0(' --', oty, ' '),
                  blout,
                  paste0(' --cpu ', n_threads),
                  ' ',
                  hmm,
                  ' ',
                  fasta)

  system(hmmse)

  return(blout)

}



#' @name readDomtblout
#' @title Process hmmsearch output to make it 'R'eadable.
#' @description Process hmmsearch output to make it readable by R.
#' @param domtblout \code{character}. hmmsearch output file, formated as a
#' \code{--domtblout} table. See HMMER 3.1b2.
#' @return A \code{data.frame} with the hmmsearch output.
#' @note Taken and adapted from \code{micropan} package (Lars Snipen and
#' Kristian Hovde Liland).
readDomtblout <- function(domtblout){

  rl <- readLines(domtblout)
  rl <- rl[which(!grepl("^\\#",rl))]
  rl <- gsub("[ ]+"," ",rl)
  lst <- strsplit(rl," ")

  query <- sapply(lst, function(x){x[1]})
  hit <- sapply(lst, function(x){x[4]})
  pfmID <- sapply(lst, function(x){x[5]})
  eval <- as.numeric(sapply(lst, function(x){x[13]}))
  score <- as.numeric(sapply(lst, function(x){x[14]}))
  st <- as.numeric(sapply(lst, function(x){x[18]}))
  en <- as.numeric(sapply(lst, function(x){x[19]}))
  desc <- sapply(lst, function(x){paste(x[23:length(x)],collapse = " ")})

  hmmer.table <- data.frame(Query=query,
                            Hit=hit,
                            PfamID=pfmID,
                            Evalue=eval,
                            Score=score,
                            Start=st,
                            End=en,
                            Description=desc,
                            stringsAsFactors = F)


  return(hmmer.table)
}

