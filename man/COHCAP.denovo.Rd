\name{COHCAP.denovo}
\alias{COHCAP.denovo}

\title{
De-Novo Identification of Differentially Methylated CpG Site Clusters
}
\description{
Identifies differentially methylated CpG sites that form clusters.

Only works with 2-group comparison result.  Does not map regions to genes or integrate with gene expression data.
}
\usage{
COHCAP.denovo(site.table, project.name, project.folder,
				min.sites = 4, max.dist = 500, output.format = "xls")
}
\arguments{
 \item{site.table}{
	Data frame with CpG site statistics (one row per CpG site) and CpG site annotations (in columns 2-5).
	
	The COHCAP.site function automatically creates this file.
}
  \item{project.name}{
	Name for COHCAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for COHCAP output files.
}
 \item{min.sites}{
	Minimum number of differentially methylated sites to define a differentially methylated CpG island.
}
 \item{max.dist}{
	Maximum distance between ordered differentially methylated CpG sites to be considered in a cluster.
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'csv' for comma-separated file, or 'txt' for tab-delimited text file.
}}

\value{
	Returns a table of clusters of differentially methylated sites, defined using criteria specified to COHCAP.site()
}

\examples{
library("COHCAP")

dir <- system.file("extdata", package="COHCAP")
beta.file <- file.path(dir,"GSE42308_truncated.txt")
sample.file <- file.path(dir,"sample_GSE42308.txt")
project.folder <- getwd()
project.name <- "450k_test"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name,
				project.folder, ref="parental")
COHCAP.denovo(filtered.sites, project.name, project.folder)
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}