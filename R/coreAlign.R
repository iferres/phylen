#' @name phylen
#' @title Compute Core Genome Alignment And Phylogeny
#' @description Identify and align core genes, concatenate the alignments
#' in a single file, and use it to compute a phylogeny.
#' @param gffs A \code{character} vector with the gff file paths.
#' @param hmmFile The path to the \code{.hmm.tar.gz} file downloaded from
#' EggNOG website or using \link{list_eggnogdb} and \link{download_nog_hmm}
#' functions on this package. The already prepared \code{hmm} text file can
#' also be provided (see \code{isCompressed} below). Alternatively a custom set
#' of hmm files can be passed as a concatenated single file.
#' @param isCompressed \code{logical()} If the \code{hmm} param points to the
#' "hmm.tar.gz" file downloaded from EggNOG, it should be set to \code{TRUE}.
#' If the pipeline has been already ran or a custom set of hmm is used, the
#' \code{hmm} parameter should point to the ".hmm" file, and the
#' \code{isCompressed} param should be set to \code{FALSE}. If the pipeline was
#' ran before, the function will also check for index files and will produce
#' them if any of the required is missing. See "hmmpress" from HMMER 3.1b2
#' manual.
#' @param eval Consider hits with an evalue less than \code{eval}.
#' @param level \code{numeric} The percentage of isolates a gene must be in to
#' be considered part of the coregenome. If nothing is specified, a plot is
#' generated at the middle of the process showing number of core genes vs a
#' percentage (from 100 to 85%), and the user is asked to choose an integer
#' in that range. The process wont continue until the user choose a level.
#' @param phyloMode One of "nj" (Neighbour-joining) of "ml" (Maximum
#' likelihood).
#' @param nbs Number of bootstrap. If \code{phyloMode} is set to "nj", this
#' parameter is ignored. If \code{phyloMode} is set to "ml", and nbs is set to
#' 0, no bootstrap is performed. If \code{phyloMode} = "ml", and \code{nbs}>0,
#' then bootstrap is performed and 2 newick files are generated, one with the
#' ML optimized tree (subffix "_ml.nwk"), and another with the bootstrap trees
#' (subffix "_ml_\code{nbs}.nwk").
#' @param outDir Where to put the output files. If \code{outDir} do not exists,
#' then a directory with the specified name is created.
#' @param aliPfx A \code{character} string with the coregenome alignment file
#' prefix. (Default: coregenome).
#' @param treePfx A \code{character} string with the newick trees files
#' prefixes. (Default: phylo).
#' @param mafftMode Alignment accuracy. One of "mafft", "ginsi", "linsi" or
#' "einsi". The first one is the default MAFFT mode, very fast. The second uses
#' mafft options "--maxiterate 1000 --globalpair". The third uses "--maxiterate
#' 1000 --localpair" (phylen DEFAULT). The fourth uses "--ep 0 --maxiterate
#' 1000 --genafpair". See MAFFT manual for more details.
#' @param keepOgs \code{logical()} If want to keep an intermediate directory
#' fasta files containing the orthologous groups (DEFAULT: \code{FALSE}). If
#' \code{TRUE}, then a directory called "orthogroups" is kept inside the
#' \code{outDir} directory.
#' @param n_threads \code{integer} The number of cpus to use.
#' @return A core genome alignment file, a phylogenetic tree in newick format
#' (or two, see \code{nbs} parameter), and an object of class "phylo" on
#' console. Optionally, a directory with the orthologous groups used for the
#' alignment (see \code{keepOgs} parameter).
#' @details This function takes gff files as returned by prokka (Seemann T,
#' 2014) and a set of hmm models, search the models in the genomes, identifies
#' the "core" set of genes, align and concatenates them into a "super gene"
#' alignemnt. Once this alignment is built, a phylogeny is inferred.
#'
#' HMMER 3.1b2 is used as search engine, and MAFFT aligner is used to align the
#' orthologous groups. Both software must be installed before running this
#' pipeline.
#'
#' \link{phangorn} package is used to perform the phylogenetic inference.
#' @importFrom parallel mclapply
#' @importFrom seqinr write.fasta
#' @importFrom graphics plot
#' @author Ignacio Ferres
#' @references Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and
#' evolution in R language. Bioinformatics 20: 289-290.
#'
#' Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4)
#' 592-593.
#'
#' Katoh K, Standley DM. 2013. MAFFT Multiple Sequence Alignment Software Version 7:
#' Improvements in Performance and Usability. Molecular Biology and Evolution
#' 30(4):772-780.
#'
#' Jensen LJ, Julien P, Kuhn M, et al. eggNOG: automated construction and annotation of
#' orthologous groups of genes. Nucleic Acids Research 36(Database issue):D250-D254.
#'
#' S. R. Eddy. 2011. Accelerated profile HMM searches. PLoS Comp. Biol. 7:e1002195.
#' @export

phylen <- function(gffs = character(),
                   hmmFile = character(),
                   isCompressed = TRUE,
                   eval = 1e-30,
                   level,
                   phyloMode = 'ml',
                   nbs = 100L,
                   outDir = 'phylen',
                   aliPfx = 'coregenome',
                   treePfx = 'phylo',
                   mafftMode = 'linsi',
                   keepOgs = FALSE,
                   n_threads = 1L){

  #Err

  if (.Platform$OS.type!='unix'){
    stop('Sorry, this package works only on unix-like platforms.')
  }

  if(Sys.which('hmmsearch')==''){
    stop('hmmsearch is not in $PATH. Please install it before running this function.')
  }

  if(Sys.which('mafft')==''){
    stop('mafft is not in $PATH. Please install it before running this function.')
  }

  if (any(!file.exists(gffs))){
    stop("One or more gff files doesn't exists.")
  }

  if(length(gffs)<2){
    stop('At least 2 gff files must be provided.')
  }

  if(!file.exists(hmmFile)){
    stop("The hmm file doesn't exists in the specified path.")
  }

  mafftMode <- match.arg(mafftMode, choices = c('mafft',
                                                'ginsi',
                                                'linsi',
                                                'einsi'))

  #wd
  if (dir.exists(outDir)){
    wd <- paste0(normalizePath(outDir), '/')
  } else{
    dir.create(outDir)
    wd <- paste0(normalizePath(outDir), '/')
  }

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
         main = 'Choose: ',
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
  ogsDirNam <- paste0(wd, 'orthogroups/')
  dir.create(ogsDirNam)
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
  outAli <- paste0(wd, aliPfx, '.aln')
  cv <- catVert(outfile = outAli,
                sos = ch)
  file.remove(ch)
  cat('DONE!\n')

  if (!keepOgs){
    cat('Removing intermediate files.. ')
    unlink(ogsDirNam, recursive = TRUE)
    cat('DONE!\n')
  }

  #phylo
  cat('Generating phylogeny..\n')
  phylo <- computePhylo(ali = cv,
                        mode = phyloMode,
                        nbs = nbs,
                        n_threads = n_threads,
                        outPrefix = treePfx,
                        outDir = wd)
  cat('Generating phylogeny.. DONE!\n\n')

  #Out
  fin <- paste0('Finished: ',
                length(ge),
                ' groups of orthologous from ',
                nrow(pm),
                ' isolates have been used in the alignment.\n',
                'Returning an object of class "phylo" with ',
                length(phylo$tree$tip.label), 'tips and ',
                phylo$tree$Nnode, ' nodes.\n')
  cat(fin)
  return(phylo)

}




