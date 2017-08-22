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
#' @param eval Consider hits with and evalue less than \code{eval}.
#' @param outfile A \code{character} string with the coregenome alignment file
#' name. (Default: coregenome.aln).
#' @param mafftMode Alignment accuracy. One of "mafft", "ginsi", "linsi" or
#' "einsi". The first one is the default MAFFT mode, very fast. The second uses
#' mafft options "--maxiterate 1000 --globalpair". The third uses "--maxiterate
#' 1000 --localpair" (phylen DEFAULT). The fourth uses "--ep 0 --maxiterate
#' 1000 --genafpair". See MAFFT manual for more details.
#' @param ogsDirNam \code{character()} The directory name where to put
#' intermediate files. If not provided or already exists, it will be
#' automatically generated. This directory is removed at the end of the
#' pipeline, unless \code{keepOgs = TRUE} (see below).
#' @param keepOgs \code{logical()} If want to keep an intermediate directory
#' fasta files containing the orthologous groups (DEFAULT: \code{FALSE}).
#' @param level \code{numeric}
#' @param n_threads \code{integer} The number of cpus to use.
#' @return A core genome alignment file.
#' @importFrom parallel mclapply
#' @importFrom seqinr write.fasta
#' @importFrom graphics plot
#' @export
#' @author Ignacio Ferres
coreAlign <- function(gffs = character(),
                      hmmFile = character(),
                      isCompressed = TRUE,
                      eval = 1e-30,
                      outfile = 'coregenome.aln',
                      mafftMode = 'linsi',
                      ogsDirNam,
                      keepOgs = FALSE,
                      level,
                      n_threads = 1L){

  #Err

  #wd
  wd <- paste0(normalizePath(dirname(outfile)), '/')

  #Decompress hmm.tar.gz, concatenate models, hmmpress
  hmm <- setHmm(hmm = hmmFile, isCompressed)

  cat('Getting information from hmms.. ')
  stats <- hmmStat(hmm[1])
  lev <- getIdsFromStats(stats)
  cat('DONE!\n')


  #Extract aa seqs from gffs
  cat('Searching.. ')
  hits <- mclapply(gffs, function(x){

    tmp <- tempdir()
    aas <- extractSeqsFromGff3(x, keep = 'none', in.path = tmp, write.in.path = 'aa')
    aas <- paste0(tmp,'/', sub('gff$','faa',rev(strsplit(x,'/')[[1]])[1]))

    blout <- hmmSearch(aas, hmm = hmm[1], n_threads = 0)
    file.remove(aas)
    m <- readDomtblout(domtblout = blout)
    file.remove(blout)
    m <- m[-which(m$Evalue>=eval),]
    sp <- split(m, m$Query)
    assig <-lapply(sp, function(y){
      sp2 <- split(y, y$Hit)
      hi <- do.call(rbind,lapply(sp2, function(z){sum(z$Score) / sum(z$End-z$Start)}))
      ma <- which.max(hi[,1])
      c(rownames(hi)[ma], hi[ma, 1])
    })
    vs <- as.data.frame(do.call(rbind, assig))
    colnames(vs) <- c('Model', 'MeanScore')
    sp <- split(vs, vs$Model)
    lp <- lapply(sp, function(y){rownames(y)[which.max(y$MeanScore)]})
    do.call(c,lp)

  }, mc.cores = n_threads, mc.preschedule = FALSE)
  cat('DONE!\n')

  #Merge
  cat('Computing panmatrix.. ')
  pm <- lapply(hits, function(x){ table(factor(names(x), levels = lev)) })
  ntb <- unlist(lapply(hits, function(x){ strsplit(x[1], ';')[[1]][1] }))
  names(pm) <- ntb
  pm <- do.call(rbind, pm)
  pm <- pm[, -which(colSums(pm)==0)]
  cat('DONE!\n')

  #Identify core-models
  if (missing(level)){
    sq <- seq(1, 0.85, -0.01)
    ev <- sapply(sq, function(x){
      length(which(colSums(pm) >= (nrow(pm)*x)))
      })
    plot(cbind(sq*100, ev),
         ylab = 'Number of core-genes',
         xlab = 'Percentage of genomes a gene must be in to be core',
         main = 'Choose:',
         xlim = rev(range(sq*100)),
         las=1)
    level <- as.integer(readline(prompt = 'Choose a percentage:'))
  }
  level <- level/100
  ge <- names(which(colSums(pm) >= nrow(pm) * level))

  #Identify core-genes
  cat('Getting core-genes.. ')
  ges <- lapply(ge, function(x){
    unlist(lapply(hits, function(y){
      as.character(y[which(names(y)%in%x)])
    }))
  })
  names(ges) <- ge

  #Extract cds from gffs
  ffns <- mclapply(gffs, function(x){
    extractSeqsFromGff3(x, keep = 'dna', write.in.path = 'none')
  }, mc.cores = n_threads, mc.preschedule = FALSE)
  ffns <- unlist(ffns, recursive = FALSE)
  cat('DONE!\n')

  #Write groups of orthologous
  cat('Writting fastas.. ')
  if (missing(ogsDirNam) || dir.exists(ogsDirNam)){
    now <- format(Sys.time(), "%b%d%H%M%S")
    ogsDirNam <- paste0(wd, 'phylen_', now, '/')
    dir.create(ogsDirNam)
  }
  gfi <- lapply(names(ges), function(x){
    ofi <- paste0(ogsDirNam, x,'.fasta')
    write.fasta(ffns[ges[[x]]],
                names = ges[[x]],
                file.out = ofi)
    ofi
  })
  gfi <- unlist(gfi)
  rm(ges)
  cat('DONE!\n')

  #Align
  cat('Aligning.. ')
  afi <- mclapply(gfi, function(x){
    mafft(infile = x, mode = mafftMode)
  }, mc.cores = n_threads, mc.preschedule = FALSE)
  afi <- unlist(afi)
  cat('DONE!\n')

  #Creates supergene
  cat('Concatenating.. ')
  #Concatenates Horizontal
  ch <- catHoriz(rn = rownames(pm),
                 ogsDirNam = ogsDirNam,
                 afi = afi,
                 extension = '_supergene.fasta',
                 n_threads = n_threads)


  #Concatenates vertical
  cv <- catVert(outfile = outfile,
                sos = ch)
  file.remove(ch)
  cat('DONE!\n')

  if (!keepOgs){
    cat('Removing intermediate files.. ')
    unlink(ogsDirNam, recursive = TRUE)
    cat('DONE!\n')
  }

  #Trim?

  #Out
  fin <- paste0('Finished: ',
                length(ge),
                ' groups of orthologous from ',
                nrow(pm),
                ' isolates have been used in the alignment.\n')
  cat(fin)
  return(cv)

}




