gaiah-wiwa
================
23 January, 2017

-   [Preliminaries](#preliminaries)
-   [Running the analyses](#running-the-analyses)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is a repository that holds code to perform the analyses in Ruegg, Anderson, Harrigan, Paxton, Kelly, Moore, and Smith (2017) "New Title" in ...

Preliminaries
-------------

First, you gotta get the `gaiah` R package. Hopefully that will be up on CRAN soon, but until then, get it from GitHub:

``` r
devtools::install_github("eriqande/gaiah")
```

Then make sure that you have all the necesary packages. Here they are, listed in their `library` calls:

    #> library(grid)
    #> library(gridExtra)
    #> library(raster)  # call this before dplyr so it doesn't mask select
    #> library(dplyr)
    #> library(stringr)
    #> library(ggplot2)
    #> library(gaiah)
    #> library(forcats)
    #> library(tikzDevice)

The `tikzDevice` package is there for putting nice LaTeX typesetting on the figures.

Running the analyses
--------------------

All the major analyses are in the file R-script at `R-main/01-wiwa-analysis.R`.

You can just run it. It can take a while (a few hours, depending on how fast your machine is). Along the way, the code will create and cache some intermediate outputs that you can use to more quickly reproduce parts of the code. This can be controlled by changing TRUE to FALSE at the lines:

``` r
#### CHOOSE WHETHER TO USE PRE-STORED VALUES OR RECOMPUTE THINGS####
COMPUTE_ISO_POSTERIORS <- TRUE  # if false then it will just load these up from a cache
RECOMPUTE_MHIA_GRID <- TRUE
RECOMPUTE_PMGCD_GRID <- TRUE
REMAKE_ALL_SUPP_MAPS <- TRUE
```

in the file `R-main/01-wiwa-analysis.R`
