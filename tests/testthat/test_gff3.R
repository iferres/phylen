context('gff3')

source(system.file('R', 'gff3_Handling.R', package = 'phylen'))

tgz <- system.file('extdata', 'toydata.tar.gz', package = 'phylen')

file <- 'H_bilis_str_AAQJH.gff'
untar(tgz, files = , exdir = tempdir())
gff <- list.files(path = tempdir(),
                  pattern = file,
                  full.names = TRUE)

rl <- readLines(gff)
fasta.ini <- grep("##FASTA",rl) + 1
fasta.end <- length(rl)
fasta <- rl[fasta.ini:fasta.end]
fheaders <- grep('>',fasta,fixed = T)
t <- matrix(nrow = length(fheaders),ncol = 2)
t[,1] <- fheaders + 1
t[,2] <- c((fheaders - 1)[-1],length(fasta))
ap <- apply(t,1,function(x){
  paste(fasta[x[1]:x[2]],collapse = '')
})
names(ap) <- gsub('>','',fasta[fheaders])
fnas <- lapply(tolower(ap),function(x){
  seqinr::as.SeqFastadna(seqinr::s2c(x),name = names(x))
})

test_that("extracting table from gff works",{
  tt <- phylen:::extractGffTable(rl)
  d <- dim(tt)
  d1 <- d[1]
  d2 <- d[2]
  expect_is(tt, class = "data.frame")
  expect_equal(d1, 2365)
  expect_equal(d2, 10)
})

test_that("extracting sequences from gff table works",{
  cds <- gfftable[1, ]
  sq <- phylen:::getFfnFaa(fnas,
                           contig = cds$Contig,
                           strand = cds$Strand,
                           from = cds$From,
                           to = cds$To,
                           id = cds$ID,
                           product = cds$Product)

  expect_is(sq, 'list')
  expect_length(sq, 2)
  expect_is(sq[[1]], "SeqFastadna")
  expect_is(sq[[2]], "SeqFastaAA")
})


file.remove(gff)
