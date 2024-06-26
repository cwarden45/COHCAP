## COHCAP: R package for DNA Methylation Analysis (and Gene Expression Integration) for Illumina Infinium Array and BS-Seq Data

```diff
+ The COHCAP package was first developed within the Bioinformatics Core, under the direction of Yate-Ching Yuan.
+ Charles Warden also provided assistance with updates and maintance that include a Bioconductor version of the package
+ while working in the Integrative Genomics Core (with Xiwei Wu as the director) as well as when he was/is not working at City of Hope.
+
+ Accordingly, Charles Warden, Yate-Ching Yuan, and Xiwei Wu were listed as authors (when a Bioconductor version was available). 
```

You can directly install from this location using [devtools](https://github.com/hadley/devtools):

If you haven't already installed devtools, you can install the package from CRAN using `install.packages("devtools")`.  Once installed, you can run the following commands:

 ```R
 library("devtools")
 devtools::install_github("cwarden45/COHCAP")
 ```

### Publications

[Warden et al. 2013](https://academic.oup.com/nar/article/41/11/e117/2411160) - original publication; not to be confused with [citation for corrigendum in 2019](https://academic.oup.com/nar/article/47/15/8335/5538015).

[Warden 2014](https://www.researchsquare.com/article/nprot-2965/v1) - Bioconductor version (most closely related to this code).

### Additional Tutorials (including Post-2016 Updates)

[COHCAP Demo for Illumina EPIC Array](https://github.com/cwarden45/DNAmethylation_templates/tree/master/EPIC_COHCAP_Demo)

### Additional Acknowledgements

 ```
For Post-2016 function additions / updates, we would like to acknowledge:
 - Yuan Yuan MD/PhD for having a project used for initial EPIC testing in the standalone version (which also caused me to debug the custom annotation function in COHCAP).
 - Susan Neuhausen / Ding Yuan Chun for projects that led to the addition of continuous variable linear regression analysis.
 - Yapeng Su (Caltech, at time of support) for having a project to test application of COHCAP to RRBS data
 - D. Joe Jerry (UMass-Amherst) for a project that was helpful for general debugging.
 - Feng Miao for having a project where I wrote partially similar code to the COHCAP de novo region boundaries for a 3C-Seq dataset (from the Natarajan Lab)
 - Shiuan Chen (application to EPIC dataset, collaboration with Susan Neuhausen)
```
