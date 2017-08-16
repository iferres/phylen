#' @name coreAlign
#' @title Compute Core Genome Alignment
#' @description Identify and align core genes, and concatenate the alignments
#' in a single file suitable for phylogenetic analyses.
#' @param gffs A \code{character} vector with the gff file paths.
#' @param hmm A \code{character} string with the path to the \code{.hmm.tar.gz}
#' file downloaded from EggNOG database website.
#' @param outfile A \code{character} string with the coregenome alignment file
#' name. (Default: coregenome.aln).
#' @param level \code{numeric}
#' @param n_threads \code{integer} The number of cpus to use.
#' @return A core genome alignment file.
#' @export
#' @author Ignacio Ferres
coreAlign <- function(gffs = character(),
                      hmm = character(),
                      outfile = 'coregenome.aln',
                      level = 1,
                      n_threads = 1L){

  #Err

  #Extract aa seqs from gffs

  #Search and filter

  #Merge

  #Process

  #Extract cds from gffs

  #Align

  #Creates supergene

  #Trim?

  #Out

}




