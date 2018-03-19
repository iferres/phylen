#' @name extractSeqsFromGff3
#' @title Extract Sequences from Gff3 file
#' @description Extract Sequences from Gff3 file and outputs a list of element,
#' one per each CDS. Each element contain the DNA gene sequence and the AA
#' protein sequence.
#' @param infile \code{character}. The gff3 filename.
#' @param in.path Where the output (fasta files) will be written.
#' @param keep \code{character}. What to keep in memory. One of "aa", "dna",
#' "both" or "none".
#' @param write.in.path \code{character}. What to write \code{in.path}. One of
#' "aa", "dna", "both" or "none".
#' @return A \code{list} with CDS gene and/or protein sequences, if specifyied.
#' @author Gregorio Iraola and Ignacio Ferres.
#' @importFrom seqinr as.SeqFastadna s2c
extractSeqsFromGff3 <- function(infile,
                                in.path=NULL,
                                keep = 'aa',
                                write.in.path='dna'){

  if (is.null(in.path)){
    in.path <- './'
  }

  if(!dir.exists(in.path)){
    stop('Directory does not exist.')
  }

  write.in.path <- match.arg(write.in.path, c('aa','dna','both','none'))

  keep <- match.arg(keep, c('aa', 'dna', 'both', 'none'))

  readLines(infile) -> rl
  if(rl[1]!="##gff-version 3"){
    stop('One or more gff files seems not to be gff3 (version 3).')
  }

  if(grepl('#',infile)){
    warning('Conflictive character "#", will be substituted by "_".')
    gsub('#','_',infile) -> infile
  }

  #Save sequence in SeqFastadna object
  if(!any(grepl("##FASTA",rl))){
    stop("One or more gff files do(es) not contain the fasta genome sequences.")
  }
  grep("##FASTA",rl) + 1 -> fasta.ini
  length(rl) -> fasta.end
  rl[fasta.ini:fasta.end] -> fasta
  grep('>',fasta,fixed = T) -> fheaders
  matrix(nrow = length(fheaders),ncol = 2) -> t
  fheaders + 1 -> t[,1]
  c((fheaders - 1)[-1],length(fasta)) -> t[,2]
  apply(t,1,function(x){
    paste(fasta[x[1]:x[2]],collapse = '')
  }) -> ap
  names(ap) <- gsub('>','',fasta[fheaders])
  lapply(tolower(ap),function(x){
    as.SeqFastadna(s2c(x),name = names(x))
  }) -> fnas

  extractGffTable(rl = rl) -> gfftable

  gfftable[which(gfftable$Type=='CDS'),] -> cds

  #Patch to avoid wrong formatted cds in gff3 file (predictions with artemis
  # return a different format than Prodigal or Aragorn, so they are discarded.
  # Usually there are very few (if any) proteins predicted with this method,
  # nothing to worry about).
  rl[which(grepl('^\\#\\#sequence-region',rl))] -> sqreg
  do.call(rbind,strsplit(sqreg,' '))[,2] -> contigs
  which(!cds$Contig%in%contigs) -> bad
  if (length(bad)>0){
    cds[-bad,] -> cds
  }


  apply(cds,1,function(x){
    getFfnFaa(fnas=fnas,
              contig = x[1],
              strand = x[9],
              from = as.numeric(x[7]),
              to = as.numeric(x[8]),id = x[2],
              product = x[5])
  }) -> fin
  nam <- sapply(fin,function(x){attr(x[[1]],'name')})

  if(any(grepl('#',nam))){
    warning('Conflictive character "#", will be substituted by "_".')
    gsub('#','_',nam) -> nam
  }
  names(fin) <- nam

  if (write.in.path == 'dna'){

    ffn <- sapply(fin,function(x){x[[1]]})
    names(ffn) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(ffn))
    write.fasta(ffn,
                names = names(ffn),
                file.out = paste0(in.path,
                                  '/',
                                  sub('.gff','.ffn',rev(strsplit(infile,'/')[[1]])[1]),
                                  collapse = '/'))

  }else if(write.in.path == 'aa'){

    faa <- sapply(fin,function(x){x[[2]]})
    names(faa) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(faa))
    write.fasta(faa,
                names = names(faa),
                file.out = paste0(in.path,
                                  '/',
                                  sub('.gff','.faa',rev(strsplit(infile,'/')[[1]])[1]),
                                  collapse = '/'))

  }else if(write.in.path == 'both'){

    ffn <- sapply(fin,function(x){x[[1]]})
    names(ffn) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(ffn))
    write.fasta(ffn,
                names = names(ffn),
                file.out = paste0(in.path,
                                  '/',
                                  sub('.gff','.ffn',rev(strsplit(infile,'/')[[1]])[1]),
                                  collapse = '/'))

    faa <- sapply(fin,function(x){x[[2]]})
    names(faa) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(faa))
    write.fasta(faa,
                names = names(faa),
                file.out = paste0(in.path,
                                  '/',
                                  sub('.gff','.faa',rev(strsplit(infile,'/')[[1]])[1]),
                                  collapse = '/'))

  }else if(write.in.path == 'none'){
    ffn <- sapply(fin,function(x){x[[1]]})
    names(ffn) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(ffn))
    faa <- sapply(fin,function(x){x[[2]]})
    names(faa) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(faa))
  }

  if (keep == 'dna'){

    if(!exists('ffn')){
      ffn <- sapply(fin,function(x){x[[1]]})
      names(ffn) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(ffn))
    }
    ffn <- lapply(ffn, function(x){paste0(x, collapse = '')})
    ffn

  }else if(keep == 'aa'){

    if (!exists('faa')){
      faa <-  sapply(fin,function(x){x[[2]]})
      names(faa) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(faa))
    }
    faa <- lapply(faa, function(x){paste0(x, collapse = '')})
    faa

  }else if (keep == 'both'){

    fin <- lapply(fin,function(x){
      lapply(x, function(y){
        paste0(y, collapse = '')
      })
    })
    names(fin) <- paste0(sub('.gff$',';',rev(strsplit(infile,'/')[[1]])[1]),names(fin))
    fin

  }else if (keep == 'none'){

    # cat('Nothing to return.\n')

  }

}

#' @name getFfnFaa
#' @title Extract CDS gene (DNA) and protein (AA) sequences
#' @description Extract CDS gene (DNA) and protein (AA) sequences from both a
#' gff3 table and the fasta sequence.
#' @param fnas A \code{list} of \code{SeqFastadna} genome sequences.
#' @param contig \code{character}. The name of the contig.
#' @param strand \code{character}. Either "+" or "-".
#' @param from \code{numeric}. The coordinate "from".
#' @param to \code{numeric}. The coordinate "to".
#' @param id \code{character}. The ID of the CDS.
#' @param product \code{character}. The gene product.
#' @return A \code{list} with two entries. The first is the CDS's DNA sequence,
#' and the last is the protein sequence. Object classes are \code{SeqFastadna}
#' and \code{SeqFastaAA}, respectively.
#' @author Ignacio Ferres
#' @importFrom seqinr as.SeqFastadna as.SeqFastaAA
getFfnFaa <- function(fnas,contig,strand,from,to,id,product){
  # seqinr::getFrag.SeqFrag(fnas[contig][[1]],begin = from,end = to) -> ffn
  fnas[contig[1]][[1]][from:to] -> ffn
  # attr(ffn,'begin') <- NULL
  # attr(ffn,'end') <- NULL
  if(strand=='-'){
    rev(comp(ffn)) -> ffn
  }
  seqinr::as.SeqFastadna(ffn,
                         name = id,
                         Annot = product) -> ffn
  seqinr::as.SeqFastaAA(translate(ffn,numcode = 11),
                        name = id,
                        Annot = product) -> faa

  list(ffn,faa) -> out
  out
}

#' @name extractGffTable
#' @title Extract the gff3 table and make it 'R'eadable.
#' @description Read a gff3 file and transforms the table in a \code{data.frame}.
#' @param rl \code{character}. A vector of character strings as passed by
#' \code{readLines()}, reading the gff3 file.
#' @return A \code{data.frame}.
#' @author Ignacio Ferres
extractGffTable <- function(rl){
  which(grepl('^\\#\\#',rl)) -> w
  rev(w)[1] - 1 -> upto
  rev(w)[2] + 1 -> from
  rl[from:upto] -> o

  strsplit(o,'\t') -> lst

  contig <- sapply(lst,function(x){x[1]})
  type <- sapply(lst,function(x){x[3]})
  from <- sapply(lst,function(x){x[4]})
  to <- sapply(lst,function(x){x[5]})
  strand <- sapply(lst,function(x){x[7]})
  phase <- sapply(lst,function(x){x[8]})
  attrib <- sapply(lst,function(x){x[9]})

  metadata <- strsplit(attrib,';')

  id <- sapply(metadata,function(x){
    gp<-grep('ID=',x,value = T)
    if(length(gp)>0){sub('ID=','',gp)}else{''}
  })
  locustag <- sapply(metadata,function(x){
    gp<-grep('locus_tag=',x,value = T)
    if(length(gp)>0){sub('locus_tag=','',gp)}else{''}
  })
  gene <- sapply(metadata,function(x){
    gp<-grep('gene=',x,value = T)
    if(length(gp)>0){sub('gene=','',gp)}else{''}
  })
  product <- sapply(metadata,function(x){
    gp<-grep('product=',x,value = T)
    if(length(gp)>0){sub('product=','',gp)}else{''}
  })

  out <- data.frame(Contig=contig,
                    ID=id,
                    LocusTag=locustag,
                    Gene=gene,
                    Product=product,
                    Type=type,
                    From=from,
                    To=to,
                    Strand=strand,
                    Phase=phase,
                    stringsAsFactors = F)
  out
}

#' @name translate
#' @title Translate
#' @description Slightly modified version of \link[seqinr]{translate} function
#' of the \code{seqinr} package.
#' @param seq See \link[seqinr]{translate}.
#' @param frame See \link[seqinr]{translate}.
#' @param sens See \link[seqinr]{translate}.
#' @param numcode See \link[seqinr]{translate}.
#' @param NAstring See \link[seqinr]{translate}.
#' @param ambiguous See \link[seqinr]{translate}.
#' @return See \link[seqinr]{translate}.
#' @importFrom seqinr s2n s2c amb
translate <- function(seq,
                      frame = 0,
                      sens = "F",
                      numcode = 1,
                      NAstring = "X",
                      ambiguous = FALSE){


  seqn <- seqinr::s2n(seq, levels = c('t', 'c', 'a', 'g'), forceToLower = FALSE) #######
  l <- 3 * ((length(seq) - frame)%/%3)
  c1 <- seq(from = frame + 1, to = frame + l, by = 3)
  tra <- 16 * seqn[c1] + 4 * seqn[c1 + 1] + seqn[c1 + 2] + 1
  code <- seqinr::s2c(seqinr::SEQINR.UTIL$CODES.NCBI$CODES[numcode])
  result <- code[tra]
  result[is.na(result)] <- NAstring
  if (ambiguous) {
    toCheck <- which(result == NAstring)
    for (i in toCheck) {
      codon <- seq[c1[i]:(c1[i] + 2)]
      allcodons <- as.vector(outer(as.vector(outer(seqinr::amb(codon[1],forceToLower = F),
                                                   seqinr::amb(codon[2], forceToLower = F), paste, sep = "")), amb(codon[3]),
                                   paste, sep = ""))
      allaminoacids <- sapply(allcodons, function(x){
        translate(seqinr::s2c(x), numcode = numcode, ambiguous = FALSE)
      })
      if (all(allaminoacids == allaminoacids[1]))
        result[i] <- allaminoacids[1]
    }
  }
  return(result)
}


#' @name comp
#' @title Complements a nucleic acid sequence
#' @description Sligthly modified version of \link[seqinr]{comp}, from
#' \code{seqinr} package.
#' @param seq See \link[seqinr]{comp}.
#' @param forceToLower See \link[seqinr]{comp}.
#' @param ambiguous See \link[seqinr]{comp}.
#' @return See \link[seqinr]{comp}.
#' @importFrom seqinr n2s s2n
comp <- function (seq,
                  forceToLower = TRUE,
                  ambiguous = FALSE){

  if (all(seq %in% LETTERS)) {
    isUpper <- TRUE
  }
  else {
    isUpper <- FALSE
  }
  # seq <- tolower(seq) #######################################################
  result <- as.vector(seqinr::n2s((3 - seqinr::s2n(seq, forceToLower=FALSE))))########
  if (ambiguous) {
    result[which(seq == "b")] <- "v"
    result[which(seq == "d")] <- "h"
    result[which(seq == "h")] <- "d"
    result[which(seq == "k")] <- "m"
    result[which(seq == "m")] <- "k"
    result[which(seq == "s")] <- "s"
    result[which(seq == "v")] <- "b"
    result[which(seq == "w")] <- "w"
    result[which(seq == "n")] <- "n"
    result[which(seq == "y")] <- "r"
    result[which(seq == "r")] <- "y"
  }
  result[which(seq == "n")] <- "n"
  if (isUpper && !forceToLower) {
    result <- toupper(result)
  }
  return(result)
}
