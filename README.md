# sybilMosek

sybilMosek is a R package, which connects [RMosek](https://docs.mosek.com/9.2/rmosek/index.html) `Mosek v9.2 x64 bit`solver to the [sybil](https://cran.r-project.org/web/packages/sybil/index.html)  package (Efficient Constrained Based Modelling in R). The package provides optimization support for the following problem types:

+ Linear Programming
+ Mixed Integer Programming
+ Quadratic Programming

SybilMosek draws inspiration from the package [sybilGurobi](https://www.cs.hhu.de/lehrstuehle-und-arbeitsgruppen/computational-cell-biology-prof-dr-martin-lercher/software-contributions/sybil.html), and interfaces RMosek optimization package to Sybil package for solving optimization problems related to Constraint Based Reconstruction and Analysis (COBRA) models.

## Prerequisite Software

These packages must be installed before installing sybilMosek

+ [sybil](http://cran.r-project.org/web/packages/sybil/index.html)
+ [Rmosek](https://docs.mosek.com/9.2/rmosek/install-interface.html). The package is tested on the **64-bit** Mosek Solver. So, please install 64 bit Mosek solver.
+ [devtools](http://cran.r-project.org/web/packages/devtools/index.html)

## Installation Instructions
SybilMosek package can be installed through _devtools_ package. The following code will check Sybil and devtools installation. If the packages are not installed, the code will install the packages.

```
list_packages <- c("sybil", "devtools")
if (length(setdiff(list_packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(list_packages, rownames(installed.packages())),
  repos = "https://cloud.r-project.org")  
}
```
The package assumes that the _licence_ for the Mosek solver is obtained and all the license setup is completed. Please check Mosek website for the [license](https://www.mosek.com/) details. To install RMosek follow the instructions provided on the Mosek website [Link](https://docs.mosek.com/9.2/rmosek/install-interface.html). Again, the package is tested on **64 bit** solver, so install x64 solver.

Now, everything is setup for installing __sybilMosek__. 

```
library(devtools)

install_github('r2rahul/sybilMosek')

```  

Alternatively,
```
devtools::install_github('r2rahul/sybilMosek')
```
## Quick Start Guide
A quick linear optimization (LP) example from Sybil vignette [1]. The problem statement

![](inst/req.png)

The equivalent sybilMosek implementation

```
lp <- optObj(solver = "sybilMosek", method = "mosek")
lp <- initProb(lp)
cm <- Matrix(c(0.5, 2, 1, 1), nrow = 2)
loadLPprob(lp, nCols = 2, nRows = 2, mat = cm,
 lb = c(0, 0), ub = rep(1000, 2), obj = c(1, 1),
 rlb = c(0, 0), rub = c(4.5, 9), rtype = c("U", "U"),
 lpdir = "max")
status <- solveLp(lp)
getObjVal(lp)
```

### Example Flux Balance Analysis (FBA)
Here is the example code for using sybilMosek for solving Flux Balance Analysis (FBA) using sybil. 

+ Let us initialize the libraries

```
library(sybil)
library(sybilMosek)
```

+ Next, we will load the example model from the package sybil. 


```
mp  <- system.file(package = "sybil", "extdata")
mod <- readTSVmod(prefix = "Ec_core", fpath = mp, quoteChar = "\"")
```

+ Finally, execute the model 

```
optL <- optimizeProb(Ec_core, solver = "sybilMosek", method = "mosek")
summaryOptsol
```

To further explore the package please go through the package vignette.

```

library(sybilMosek)

browseVignettes(package = "sybilMosek")

vignette('quickstartsybilMosekR', package = "sybilMosek")

```
## TODO
+ Update Vignettes
+ Add more COBRA tests.

## Acknowledgement
The package was originally conceived during my post-doctoral research. Initially, the package was for my personal research. Over the years I further refined the code to be used by the users of Sybil ecosystem of packages. 

## References
[1] Gelius-Dietrich, G., Amer Desouki, A., Fritzemeier, C.J. and Lercher, M.J. "sybil – Efficient constraint-based modelling in R". BMC Systems Biology, 2013. 7:125. doi:10.1186/1752-0509-7-125

Article: http://www.biomedcentral.com/1752-0509/7/125
