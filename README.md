## COHCAP: Bioconductor package for DNA Methylation Analysis (and Gene Expression Integration) for Illumina Infinium Array and BS-Seq Data

For most users, the easiest way to install COHCAP is to use the Bioconductor Release Branch:

https://www.bioconductor.org/packages/release/bioc/html/COHCAP.html

However, if you are trying to use some functions in the development version and are having issues installing the development version of COHCAP, you can try directly installing from this location using [devtools](https://github.com/hadley/devtools):

If you haven't already installed devtools, you can install the package from CRAN using `install.packages("devtools")`.  Once installed, you can run the following commands:

 ```R
 library("devtools")
 devtools::install_github("cwarden45/COHCAP")
 ```

 If you ever need to verbally describe COHCAP, it is pronounced "co-cap."