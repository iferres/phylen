[![Travis-CI Build Status](https://travis-ci.org/iferres/phylen.svg?branch=master)](https://travis-ci.org/iferres/phylen)

phylen
================

Phylen is an R package that performs automatic phylogenetic reconstruction given a set of Hidden Markov Models (HMMs). Genomes are screened against these HMMs, genes found in all genomes ("core genes") are aligned individually, those alignments are concatenated into a single supergene alignment, and a phylogenetic reconstruction is performed and returned as an object of class "phylo" so it can be further analysed using ape/phangorn framework in R. Functions to download well curated HMMs from clade-specific orthologous sets from the EggNOG database are provided although any custom set of HMMs can be used as well.

## Installation

The easiest way to install this package is using `devtools` package:

``` r
if(!require(devtools)){
   install.packages("devtools")
}
devtools::install_github("iferres/phylen")
```

### Requirements

`phylen` depends on [HMMER 3.1b2](http://hmmer.org/) and [MAFFT](https://mafft.cbrc.jp/alignment/software/). You should have them installed on you `$PATH` variable prior to using this software. Besure that the following MAFFT aliases are also installed: `linsi` (alias for `mafft --maxiterate 1000 --localpair`), `ginsi` (alias for `mafft --maxiterate 1000 --globalpair`) and `einsi` (alias for `mafft --ep 0 --maxiterate 1000 --genafpair`). This shortcuts are installed together with `mafft` if you download and install the software from the MAFFT webpage (link above), but probably not if use a package manager as `apt` or `brew`.

It also depends on [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html) package, which in turn depends on `igraph`. Some system requirements are needed to install the latter, please check them [here](http://igraph.org/r/).

## Standard workflow

This tutorial begins with the extraction of toy data attached on this package. It consists in 10 genomes of the *Campylobacterales* order, 5 of them are *Campylobacter* species and the other 5 are *Helicobacter* species. 

**Note:** `phylen` input files must be in **gff3** format, as returned by the [prokka](https://github.com/tseemann/prokka) annotation software.

``` r
# List the attached tar.gz file
tgz <- system.file('extdata', 'toydata.tar.gz', package = 'phylen')
# List the files inside it
gff <- untar(tarfile = tgz, exdir = getwd(), list = T)
# Decompress on current working directory
untar(tarfile = tgz,exdir = getwd())
# List the decompressed files
gff
```

    ##  [1] "C_fetus_subsp_testudinum_Sp3.gff"                   
    ##  [2] "C_fetus_subsp_venerealis_str_84-112.gff"            
    ##  [3] "C_hyointestinalis_subsp_lawsonii_CCUG_27631.gff"    
    ##  [4] "C_iguaniorum_str_RM11343.gff"                       
    ##  [5] "C_pinnipediorum_subsp_pinnipediorum_str_RM17262.gff"
    ##  [6] "H_bilis_str_AAQJH.gff"                              
    ##  [7] "H_himalayensis_str_YS1.gff"                         
    ##  [8] "H_pylori_J166.gff"                                  
    ##  [9] "H_pylori_str_2018.gff"                              
    ## [10] "H_typhlonius_str_1.gff"

Let's load the `phylen` package and list the available Hidden Markov Models (HMMs) on the [EggNOG](http://eggnogdb.embl.de/#/app/home) database.

``` r
# Load the phylen package
library(phylen)
```

    ## Loading required package: phangorn

    ## Loading required package: ape

``` r
# List first 50 available sets of EggNOG
list_eggnogdb()[1:50, ]
```

    ##                level.name nog.prefix
    ## 1                    LUCA        NOG
    ## 2           Acidobacteria     aciNOG
    ## 3          Acidobacteriia    acidNOG
    ## 4            Aconoidasida     acoNOG
    ## 5          Actinobacteria     actNOG
    ## 6              Agaricales     agaNOG
    ## 7          Agaricomycetes    agarNOG
    ## 8             Apicomplexa     apiNOG
    ## 9    Proteobacteria_alpha    aproNOG
    ## 10              Aquificae     aquNOG
    ## 11                Archaea      arNOG
    ## 12           Archaeoglobi     arcNOG
    ## 13             Arthropoda     artNOG
    ## 14      Arthrodermataceae    arthNOG
    ## 15             Ascomycota     ascNOG
    ## 16                   Aves     aveNOG
    ## 17                Bacilli     bacNOG
    ## 18               Bacteria    bactNOG
    ## 19            Bacteroidia   bacteNOG
    ## 20          Basidiomycota     basNOG
    ## 21          Bacteroidetes    bctoNOG
    ## 22              Bilateria      biNOG
    ## 23    Proteobacteria_beta    bproNOG
    ## 24            Brassicales     braNOG
    ## 25              Carnivora     carNOG
    ## 26          Chaetomiaceae     chaNOG
    ## 27               Chlorobi     chlNOG
    ## 28             Chlamydiae    chlaNOG
    ## 29            Chloroflexi    chloNOG
    ## 30            Chloroflexi   chlorNOG
    ## 31            Chlorophyta  chloroNOG
    ## 32               Chordata    chorNOG
    ## 33            Chromadorea     chrNOG
    ## 34             Clostridia     cloNOG
    ## 35               Coccidia     cocNOG
    ## 36          Crenarchaeota     creNOG
    ## 37      Cryptosporidiidae     cryNOG
    ## 38          Cyanobacteria     cyaNOG
    ## 39             Cytophagia     cytNOG
    ## 40      Debaryomycetaceae     debNOG
    ## 41        Deferribacteres     defNOG
    ## 42      Dehalococcoidetes     dehNOG
    ## 43     Deinococcusthermus     deiNOG
    ## 44          delta/epsilon     delNOG
    ## 45                Diptera     dipNOG
    ## 46        Dothideomycetes     dotNOG
    ## 47   Proteobacteria_delta    dproNOG
    ## 48          Drosophilidae     droNOG
    ## 49 Proteobacteria_epsilon    eproNOG
    ## 50        Erysipelotrichi     eryNOG

The order *Campylobacterales* belongs to the class *Epsilonproteobacteria*, which is listed at row number 49. The corresponding set of HMMs is, then, `eproNOG`.

We will download the HMMs directly from the EggNOG servers using the `download_nog_hmm()` function.

``` r
hmm <- download_nog_hmm(nog.prefix = "eproNOG", onDir = getwd())
hmm
```

    ## [1] "/home/iferres/Documents/eproNOG.hmm.tar.gz"

Now we have the HMMs, so we can proceed to run the main function in order to obtain a core genome alignment, and a phylogeny. This would take ~20-30 min using 4 CPUs.

``` r
p <- phylen(gffs = gff, # The gff files extracted on the first step
            hmmFile = hmm, # The downloaded HMMs
            isCompressed = TRUE, # The downloaded HMMs are compressed (tar.gz)
            phyloMode = "ml", # nj or ml, in this case we will use maximum likelihood
            nbs = 100, # The number of bootstrap.
            level = 100, # The percentage of genomes a gene must be in to be considered as part of the coregenome.
            outDir = "phylen", # Lets create a directory called "phylen" to put the output files
            n_threads = 4) # Use 4 threads
```

    ## Decompressing.. DONE!
    ## Concatenating.. 
    ## ===========================================================================
    ## DONE!
    ## Pressing.. DONE!
    ## Getting information from hmms.. DONE!
    ## Searching.. DONE!
    ## Computing panmatrix.. DONE!
    ## Getting core-genes.. DONE!
    ## Writting fastas.. DONE!
    ## Aligning.. DONE!
    ## Concatenating.. DONE!
    ## Removing intermediate files.. DONE!
    ## Generating phylogeny..
    ## optimize edge weights:  -5958637 --> -5731557 
    ## optimize base frequencies:  -5731557 --> -5686007 
    ## optimize rate matrix:  -5686007 --> -5621016 
    ## optimize edge weights:  -5621016 --> -5619537 
    ## optimize base frequencies:  -5619537 --> -5616203 
    ## optimize rate matrix:  -5616203 --> -5615171 
    ## optimize edge weights:  -5615171 --> -5615158 
    ## optimize base frequencies:  -5615158 --> -5614921 
    ## optimize rate matrix:  -5614921 --> -5614863 
    ## optimize edge weights:  -5614863 --> -5614861 
    ## optimize base frequencies:  -5614861 --> -5614846 
    ## optimize rate matrix:  -5614846 --> -5614842 
    ## optimize edge weights:  -5614842 --> -5614842 
    ## optimize base frequencies:  -5614842 --> -5614841 
    ## optimize rate matrix:  -5614841 --> -5614841 
    ## optimize edge weights:  -5614841 --> -5614841 
    ## optimize base frequencies:  -5614841 --> -5614841 
    ## optimize rate matrix:  -5614841 --> -5614841 
    ## optimize edge weights:  -5614841 --> -5614841

    ## Generating phylogeny.. DONE!
    ## 
    ## Finished: 658 groups of orthologous from 10 isolates have been used in the alignment.
    ## Returning an object of class "phylo" with 10 tips and 8 nodes.


Finally we have a "phylo" object, and we can continue analizing it with `phangorn` and `ape` packages.

``` r
# Print it
p
```

    ## 
    ## Phylogenetic tree with 10 tips and 8 internal nodes.
    ## 
    ## Tip labels:
    ##  C_fetus_subsp_testudinum_Sp3, C_fetus_subsp_venerealis_str_84-112, C_hyointestinalis_subsp_lawsonii_CCUG_27631, C_iguaniorum_str_RM11343, C_pinnipediorum_subsp_pinnipediorum_str_RM17262, H_bilis_str_AAQJH, ...
    ## Node labels:
    ##  100, 100, 100, 100, 100, 100, ...
    ## 
    ## Unrooted; includes branch lengths.

``` r
# Class?
class(p)
```

    ## [1] "phylo"

``` r
# Plot it
plot(p, type = 'unrooted', cex = 0.7, lab4ut = 'axial')
```

![](vignettes/readme_img1.png)

## Any questions?

If you have any question about usage, implementation, or result interpretation, please feel free to ask. Send them to the [issue tracker](https://github.com/iferres/phylen/issues) or directly to my e-mail (iferres@pasteur.edu.uy). If you think your concern could help others, we encourage you to use the first channel of communication.

## Contributing

If you want to contribute to the development of this software, please refer to the [contributing guidelines](CONTRIBUTING.md).

## Citation

Phylen manuscript is under review in [JOSS](https://joss.theoj.org/).

This package has been successfully used in already published papers, see `paper.md` on this repository for details.
