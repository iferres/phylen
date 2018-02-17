
#' @name computePhylo
#' @title Compute Phylogeny
#' @description Takes the core genes alignment and produce a phylogeny.
#' @param ali The path to the alignment file.
#' @param mode One of "nj" (Neighbour-joining) of "ml" (Maximum likelihood).
#' @param nbs The number of bootstraps.
#' @param n_threads The number of cpus to use. Just used if nbs > 0.
#' @param outPrefix A prefix to use on the output files.
#' @param outDir Where to put the newick files.
#' @param ... Further arguments to pass to \link[phangorn]{optim.pml}.
#' @details Takes the core gene alignment and produces newick files. If \code{mode}
#' is set to "nj", the outfile will have a "_nj.nwk" subffix. If set to "ml" and
#' nbs = 0L, then it will be "_ml.nwk", and if "ml" and nbs > 0L, then it will
#' output a second file with the bootstrap trees with the suffix
#' "_ml_\code{nbs}.nwk".
#'
#' Multiple cpus can be used when performing bootstrap.
#'
#' It uses the "GTR" model to optimize the tree.
#' @return An object of class "phylo" (phangorn and ape packages), and the
#' specified newick files.
#' @author Ignacio Ferres
#' @importFrom phangorn read.phyDat dist.hamming NJ pml optim.pml bootstrap.pml plotBS
#' @importFrom ape write.tree
#' @importFrom grDevices pdf dev.off
#' @references Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and
#' evolution in R language. Bioinformatics 20: 289-290.
#' Schliep K.P. 2011. phangorn: phylogenetic analysis in R. Bioinformatics, 27(4)
#' 592-593.
computePhylo <- function(ali,
                         mode = 'ml',
                         nbs = 100L,
                         n_threads = 1,
                         outPrefix='phylo',
                         outDir='.',
                         ...){

  outDir <- normalizePath(outDir)

  mode <- match.arg(mode, c('ml','nj'))

  dat <- read.phyDat(ali, type = 'DNA', format = 'fasta')
  dm <- dist.hamming(dat)
  nj.tree <- NJ(dm)

  if (mode=='ml'){
    nwk <- paste0(outDir, '/', outPrefix, '_ml.nwk')
    fitNJ <- pml(nj.tree, dat, model='GTR')
    fit <- optim.pml(object = fitNJ,
                     model = 'GTR',
                     ...)

    if (nbs > 0L){
      nwks <- paste0(outDir, '/', outPrefix, '_mlbs_', nbs, '.nwk')
      bs <- bootstrap.pml(fit,
                          bs = nbs,
                          optNni = TRUE,
                          multicore = TRUE,
                          mc.cores = n_threads)

      write.tree(fit$tree, file = nwk)
      write.tree(bs, file = nwks)

      pdf(file = NULL)
      res <- plotBS(fit$tree,
                    BStrees = bs)
      dof <- dev.off()

      return(res) #class 'phylo'


    }else{

      write.tree(fit$tree, file = nwk)
      return(fit$tree) #class 'phylo'

    }

  }else{

    nwk <- paste0(outDir, '/', outPrefix, '_nj.nwk')
    write.tree(nj.tree, file = nwk)
    return(nj.tree) #Class 'phylo'

  }


}
