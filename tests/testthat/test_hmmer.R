context('hmmer')


example_file <- system.file('testdata', 'example_domtblout.tar.gz', package = 'phylen')
untar(example_file, files = 'example_domtblout.tab', exdir = tempdir())
example_domtblout <- list.files(path = tempdir(),
                                pattern = '^example_domtblout.tab$',
                                full.names = TRUE)

xcln <- c("Query", "Hit", "PfamID", "Evalue", "Score", "Start", "End",
          "Description")

xclc <- c("character", "character", "character", "numeric", "numeric",
          "numeric", "numeric", "character")

test_that('readDomtblout works',{
  x <- phylen:::readDomtblout(example_domtblout)
  expect_is(x, 'data.frame')
  d <- dim(x)
  expect_identical(d[1], 26L)
  expect_identical(d[2], 8L)
  expect_identical(colnames(x), xcln)
  clc <- sapply(1:ncol(x), function(y){ class(x[, y]) })
  expect_equivalent(clc, xclc)
})

file.remove(example_domtblout)
