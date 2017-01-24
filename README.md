gaiah-wiwa
================
24 January, 2017

-   [Preliminaries](#preliminaries)
-   [Running the analyses](#running-the-analyses)
-   [Transferring figures for the publication](#transferring-figures-for-the-publication)
-   [Making supplement 1](#making-supplement-1)

<!-- README.md is generated from README.Rmd. Please edit that file -->
This is a repository that holds code to perform the analyses in Ruegg, Anderson, Harrigan, Paxton, Kelly, Moore, and Smith (2017) "New Title" in ...

Preliminaries
-------------

First, you gotta get the `gaiah` R package. Hopefully that will be up on CRAN soon, but until then, get it from GitHub:

``` r
devtools::install_github("eriqande/gaiah")
```

By the way, this is known to run on `gaiah` commit c63cce5f63d49ce8010fd560cb2394df6fb96c8e.

Just in case `devtools::install_github()` doesn't recursively get all the dependencies, you may as well install them all afresh like this:

``` r
install.packages(c("dplyr", 
                   "geosphere", 
                   "ggplot2", 
                   "lazyeval", 
                   "magrittr", 
                   "maptools", 
                   "raster", 
                   "rgeos", 
                   "sp", 
                   "stringr", 
                   "tidyr"))
```

Then make sure that you have all the necesary additional packages to run the analysis. Here they are, listed in their `library` calls:

    #> library(grid)
    #> library(gridExtra)
    #> library(raster)  # call this before dplyr so it doesn't mask select
    #> library(dplyr)
    #> library(stringr)
    #> library(ggplot2)
    #> library(gaiah)
    #> library(forcats)
    #> library(tikzDevice)

So, you can get all those like this:

``` r
install.packages(c("grid", "gridExtra", "forcats", "tikzDevice"))
```

The `tikzDevice` package is there for putting nice LaTeX typesetting on the figures. In fact, some of the code in `R-main/01-wiwa-analysis.R` makes system calls to LaTeX to typeset those figures, so if you want reproduce everything entirely, you are going to need `pdflatex` on your system path and the capacity for `R` to make calls using the `system()` function.

Running the analyses
--------------------

All the major analyses are in the R-script at `R-main/01-wiwa-analysis.R`.

You can just run it. This script creates some intermediate output files so that the computationally intensive portions of the analysis can be avoided should you need to re-run everything.
After you have run through the code once, you can change the TRUEs to FALSEs in the following lines in `R-main/01-wiwa-analysis.R` to use the cached results rather than the recomputing everything from the beginning:

``` r
#### CHOOSE WHETHER TO USE PRE-STORED VALUES OR RECOMPUTE THINGS####
COMPUTE_ISO_POSTERIORS <- TRUE  # if false then it will just load these up from a cache
REMAKE_ALL_SUPP_MAPS <- TRUE
RECOMPUTE_PMGCD_GRID <- TRUE
```

The intermediate files that are produced in each of the relevant code sections are as follows, along with the approximate compute time:

<table style="width:100%;">
<colgroup>
<col width="22%" />
<col width="56%" />
<col width="20%" />
</colgroup>
<thead>
<tr class="header">
<th>Flag Option</th>
<th>output file</th>
<th>approx compute time</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>COMPUTE_ISO_POSTERIORS</td>
<td>outputs/isotope_ref_bird_results.rda</td>
<td>1 hour</td>
</tr>
<tr class="even">
<td>REMAKE_ALL_SUPP_MAPS</td>
<td>outputs/birdmaps/bird_XX.pdf (where XX is the bird's ID)</td>
<td>20 minutes</td>
</tr>
<tr class="odd">
<td>RECOMPUTE_PMGCD_GRID</td>
<td>outputs/128_pmgcd_vals.rds</td>
<td>10 minutes</td>
</tr>
</tbody>
</table>

Setting each of those flags to FALSE (once the intermediate files in `R-main/01-wiwa-analysis.R` makes it possible to get through the whole analysis in about 2 to 2.5 minutes on a reasonably fast laptop.

Transferring figures for the publication
----------------------------------------

I include `R-main/02-move-figures.R` so that the interested can see which figures produced in the previous step get placed in the document itself. In order for this to actually work you will need to have the appropriate TeX source files, etc in the right directory inside `git-overleaf-repos`. You also need to have `pdfcrop` installed on your system.

Making supplement 1
-------------------

This is done with `R-main/03-compile-supp-1.R`. It just makes a LaTeX file that can then be typeset into the supplement.
