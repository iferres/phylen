context('Compute phylogenetic tree')

xp <- structure(list(edge = structure(c(8L, 8L, 7L, 7L, 6L, 6L, 6L,
                                        4L, 5L, 8L, 1L, 3L, 2L, 7L),
                                      .Dim = c(7L, 2L)),
                     edge.length = c(0.0833333333333333,
                                     0.0416666666666667,
                                     0.1875,
                                     0.1875,
                                     0.125,
                                     0,
                                     0.0625),
                     tip.label = c("A","B", "C", "D", "E"),
                     Nnode = 3L),
                .Names = c("edge", "edge.length", "tip.label", "Nnode"),
                class = "phylo", order = "postorder")


ali_file <- system.file('testdata', 'ali.tar.gz', package = 'phylen')
untar(ali_file, files = 'ali.fasta', exdir = tempdir())
ali <- list.files(path = tempdir(),
                  pattern = '^ali.fasta$',
                  full.names = TRUE)

test_that('phylogeny inference works', {
  p <- phylen:::computePhylo(ali, mode = 'nj', nbs = 0, outDir = tempdir())
  expect_is(p, 'phylo')
  expect_identical(p$tip.label, c('A', 'B', 'C', 'D', 'E'))
  expect_identical(p$Nnode, 3L)
  expect_equivalent(p, xp)
})

file.remove(ali)
nwk <- grep('/phylo_nj.nwk$', list.files(path = tempdir(), full.names = TRUE), value = T)
file.remove(nwk)
