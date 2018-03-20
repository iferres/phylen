context('utils')

stat_file <- system.file('testdata', 'hmmstat.tar.gz', package = 'phylen')
untar(stat_file, files = 'hmmstat.txt', exdir = tempdir())
stat <- list.files(path = tempdir(),
                   pattern = '^hmmstat.txt$',
                   full.names = TRUE)

xpec <- c("thaNOG.ENOG411CCBB.meta_raw", "thaNOG.ENOG411CCBA.meta_raw",
          "thaNOG.ENOG411CCB9.meta_raw", "thaNOG.ENOG411CCB8.meta_raw")

test_that('getIdsFromStats works',{
  x <- phylen:::getIdsFromStats(stats = stat)
  expect_is(x, 'character')
  expect_length(x, 4L)
  expect_identical(x, xpec)
})

file.remove(stat)
