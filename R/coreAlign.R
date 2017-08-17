#' @name coreAlign
#' @title Compute Core Genome Alignment
#' @description Identify and align core genes, and concatenate the alignments
#' in a single file suitable for phylogenetic analyses.
#' @param gffs A \code{character} vector with the gff file paths.
#' @param hmmFile The path to the \code{.hmm.tar.gz} file downloaded from
#' EggNOG website, or the already prepared \code{hmm} text file. See
#' \code{isCompressed} below.
#' @param isCompressed \code{logical()} If the \code{hmm} param points to the
#' "hmm.tar.gz" file, it should be set to \code{TRUE}. If the pipeline has been
#' already ran, the \code{hmm} parameter should point to the ".hmm" file that
#' was generated before, and the \code{isCompressed} param should be set to
#' \code{FALSE}. If the pipeline was ran before, the function will also check
#' for index files and will produce them if any of the required is missing. See
#' "hmmpress" from HMMER 3.1b2 manual.
#' @param outfile A \code{character} string with the coregenome alignment file
#' name. (Default: coregenome.aln).
#' @param level \code{numeric}
#' @param n_threads \code{integer} The number of cpus to use.
#' @return A core genome alignment file.
#' @importFrom parallel mclapply splitIndices
#' @importFrom seqinr write.fasta
#' @export
#' @author Ignacio Ferres
coreAlign <- function(gffs = character(),
                      hmmFile = character(),
                      isCompressed = TRUE,
                      outfile = 'coregenome.aln',
                      level = 1,
                      n_threads = 1L){

  #Err

  #Decompress hmm.tar.gz, concatenate models, hmmpress
  hmm <- setHmm(hmm = hmmFile, isCompressed)

  #Extract aa seqs from gffs
  hits <- mclapply(gffs, function(x){

    tmp <- tempdir()
    aas <- extractSeqsFromGff3(x, keep = 'none', in.path = tmp, write.in.path = 'aa')
    aas <- paste0(tmp,'/', sub('gff$','faa',rev(strsplit(x,'/')[[1]])[1]))

    blout <- hmmSearch(aas, hmm = hmm[1], n_threads = 0)
    m <- readDomtblout(domtblout = blout)



  }, mc.cores = n_threads, mc.preschedule = FALSE)



  #Merge

  #Process

  #Extract cds from gffs

  #Align

  #Creates supergene

  #Trim?

  #Out

}




