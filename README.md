# phylen
Automatic phylogenetic reconstruction using the EggNOG database (or others) hmm files.

## Dependencies
To work, it is required to have installed on your $PATH variable the following software:
 * [Mafft](http://mafft.cbrc.jp/alignment/software/)
 * [HMMER 3.1b2](http://hmmer.org/download.html)
### Other external dependencies

 * A [EggNOG](http://eggnogdb.embl.de/#/app/downloads) hmm database file. Its not required to decompress (`untar`) those files before running the pipeline.
 
## Installation
```
devtools::install_github('iferres/phylen')
```

## Usage
```r
coreAlign(gffs = character(), hmmFile = character(), isCompressed = TRUE,
  eval = 1e-30, outfile = "coregenome.aln", mafftMode = "linsi",
  ogsDirNam, keepOgs = FALSE, level, n_threads = 1L)
```
### Arguments
 * `gffs`:     A `character()` vector with the gff file paths.
 * `hmmFile`:  The path to the `.hmm.tar.gz` file downloaded from EggNOG website, or the already prepared `.hmm` text file. See `isCompressed` below.
 * `isCompressed`: `logical()` If the `hmm` param points to the "hmm.tar.gz" file, it should be set to `TRUE`. If the pipeline has been already ran, the `hmm` parameter should point to the ".hmm" file that was generated before, and the `isCompressed` parameter should be set to `FALSE`. If the pipeline was ran before, the function will also check for index files and will produce them if any of the required is missing. See "hmmpress" from HMMER 3.1b2 manual.
 * `eval`:     Consider hits with and evalue less than this.
 * `outfile`:  A `character()` string with the coregenome alignment file name. (Default: "coregenome.aln").
 * `mafftMode`: Alignment accuracy. One of "mafft", "ginsi", "linsi" or "einsi". The first one is the default MAFFT mode, very fast but not so accurate. The second uses mafft options "--maxiterate 1000 --globalpair"; accurate but slow. The third uses "--maxiterate 1000 --localpair" (phylen DEFAULT); very accurate but slow. The fourth uses "--ep 0 --maxiterate 1000 --genafpair". See MAFFT manual for more details.
 * `ogsDirNam`: `character()` The directory name where to put intermediate files. If not provided or already exists, it will be automatically generated. This directory is removed at the end of the pipeline, unless `keepOgs = TRUE` (see below).
 * `keepOgs`:  `logical()` If want to keep an intermediate directory fasta files containing the orthologous groups (DEFAULT: `FALSE`).
 * `level`:    `numeric()` The percentage of the genomes at which a gene must be present to be considered as part of the core genome. If nothing is provided, a plot will be generated showing the number of core-genes at different levels, and the user will be asked to choose an integer between 100 - 85.
 * `n_threads`: `integer()` The number of cpus to use.
 
### Value
A core genome alignment file.

### Description
Identify and align core genes, and concatenate the alignments in a single file suitable for phylogenetic analyses.

## Note:

The `phylen` package has been designed for UNIX-like platforms only.
