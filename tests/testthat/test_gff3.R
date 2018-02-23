context('gff3')

source(system.file('R', 'gff3_Handling.R', package = 'phylen'))

tgz <- system.file('extdata', 'toydata.tar.gz', package = 'phylen')

file <- 'H_bilis_str_AAQJH.gff'
untar(tgz, files = , exdir = tempdir())
gff <- list.files(path = tempdir(),
                  pattern = file,
                  full.names = TRUE)

rl <- readLines(gff)


test_that("extracting table from gff",{
  tt <- phylen:::extractGffTable(rl)
  d <- dim(tt)
  d1 <- d[1]
  d2 <- d[2]
  expect_is(tt, class = "data.frame")
  expect_equal(d1, 2365)
  expect_equal(d2, 10)
})

test_that("")


file.remove(gff)
