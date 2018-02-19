# phylen
Automatic phylogenetic reconstruction using the EggNOG database (or others) hmm files.

## Dependencies
To work, it is required to have installed on your $PATH variable the following software:
 * [Mafft](http://mafft.cbrc.jp/alignment/software/)
 * [HMMER 3.1b2](http://hmmer.org/download.html)
### Other external dependencies

 * A [EggNOG](http://eggnogdb.embl.de/#/app/downloads) hmm database file. Its not required to decompress (`untar`) those files before running the pipeline.
or..
 * A custom set of hmms concatenated into one single file.
 
## Installation
```
devtools::install_github('iferres/phylen')
```

## `phylen`: Compute Core Genome Alignment And Phylogeny
### Usage
```r
phylen(gffs = character(), hmmFile = character(), isCompressed = TRUE,
  eval = 1e-30, level, phyloMode = "ml", nbs = 100L, outDir = "phylen",
  aliPfx = "coregenome", treePfx = "phylo", mafftMode = "linsi",
  keepOgs = FALSE, n_threads = 1L)
```
### Description
Identify and align core genes, concatenate the alignments in a single file, and use it to compute a phylogeny.

### Arguments
 * `gffs`: A `character` vector with the gff file paths.
 * `hmmFile`: The path to the ".hmm.tar.gz" file downloaded from EggNOG website or using `list_eggnogdb` and `download_nog_hmm` functions on this package. The already prepared hmm text file can also be provided (see `isCompressed` below). Alternatively a custom set of hmm files can be passed as a concatenated single file.
 * `isCompressed`: `logical` If the hmm param points to the "hmm.tar.gz" file downloaded from EggNOG, it should be set to `TRUE`. If the pipeline has been already ran or a custom set of hmm is used, the hmm parameter should point to the ".hmm" file, and the `isCompressed` param should be set to `FALSE`. If the pipeline was ran before, the function will also check for index files and will produce them if any of the required is missing. See "hmmpress" from HMMER 3.1b2 manual.
 * `eval`: Consider hits with an evalue less than eval.
 * `level`: `numeric` The percentage of isolates a gene must be in to be considered part of the coregenome. If nothing is specified, a plot is generated at the middle of the process showing number of core genes vs a percentage (from 100% to 85%) in that range. The process wont continue until the user choose a level.
 * `phyloMode`: One of `"nj"` (Neighbour-joining) of `"ml"` (Maximum likelihood).
 * `nbs`: Number of bootstrap. If `phyloMode` is set to `"nj"`, this parameter is ignored. If `phyloMode` is set to `"ml"`, and `nbs` is set to `0`, no bootstrap is performed. If `phyloMode` = `"ml"`, and `nbs` > 0, then bootstrap is performed and 2 newick files are generated, one with the ML optimized tree (subffix "\_ml.nwk"), and another with the bootstrap trees (subffix "\_ml\_nbs.nwk").
 * `outDir`: Where to put the output files (DEFAULT: "."). If `outDir` do not exists, then a directory with the specified name is created.
 * `aliPfx`: A character string with the coregenome alignment file prefix. (Default: `"coregenome"`).
 * `treePfx`: A character string with the newick trees files prefixes. (Default: `"phylo"`).
 * `mafftMode`: Alignment accuracy. One of `"mafft"`, `"ginsi"`, `"linsi"` or `"einsi"`. The first one is the default MAFFT mode, very fast. The second uses mafft options `"–maxiterate 1000 –globalpair"`. The third uses `"–maxiterate 1000 –localpair"` (accurate, phylen DEFAULT). The fourth uses `"–ep 0 –maxiterate 1000 –genafpair"`. See MAFFT manual for more details.
 * `keepOgs`: `logical` If want to keep an intermediate directory fasta files containing the orthologous groups (DEFAULT: `FALSE`). If `TRUE`, then a directory called "orthogroups" is kept inside the `outDir` directory.
 * `n_threads`: `integer` The number of cpus to use.
 
### Value
A core genome alignment file, a phylogenetic tree in newick format (or two, see `nbs` parameter), and an object of class `phylo` on console (see [phangorn]{https://github.com/KlausVigo/phangorn} and [ape]{http://ape-package.ird.fr/}). Optionally, a directory with the orthologous groups used for the alignment (see `keepOgs` parameter).

## Other useful functions
The following package functions conects with EggNOG webpage and return information from it:

## `list_eggnogdb()`: List available EggNOG db datasets
### Usage
```r
list_eggnogdb()
```
### Description
This function list all the available EggNOG datasets to download, associated with their nog prefix.

### Value
A `data.frame` with the available datasets and their respective "nog" prefixes.

## `download_nog_hmm()`: Download specific hmm datasets from EggNOG.
### Usage
```r
download\_nog\_hmm(nog.prefix = "proNOG", onDir = ".")
```
### Description
Takes a nog prefix and downloads a compressed (.tar.gz) file from EggNOG db containing the hmm models of the specified group.

### Value
A tar.gz compressed file on the specified directory. On console it returns the path to the downloaded file.

## Note:

The `phylen` package has been designed for UNIX-like platforms only.



