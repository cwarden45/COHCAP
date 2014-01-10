\name{COHCAP.integrate.avg.by.island}
\alias{COHCAP.integrate.avg.by.island}

\title{
Integration with Gene Expression Data(Average by Island Workflow).
}
\description{
Provides lists of genes with a significant negative correlation between DNA methylation and gene expression data.

A table of normalizated intensity / expression values is provided in the gene expression table and a table of filtered beta values is provided for the DNA methylation data.

Lists of genes with negative expression trends will be created in the "Integrate" folder, along with scatter plots (if descired).
All correlation stats are provided in the "Raw_Data" folder.
}
\usage{
COHCAP.integrate.avg.by.island(beta.table, project.name, project.folder,
				expr.file, sample.file, cor.pvalue.cutoff=0.05,
				cor.fdr.cutoff = 0.05, cor.cutoff = -0.2, plot.scatter=TRUE,
				output.format = "xls")
}
\arguments{
 \item{beta.table}{
	Data frame with beta values averaged across differentially methylated sites (islands in rows, samples in columns).  This table is already filtered for differentially methylated CpG islands.
	
	The COHCAP.avg.by.island function automatically creates this file.
}
  \item{project.name}{
	Name for COHCAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for COHCAP output files
}
  \item{expr.file}{
	Table of normalized expression or intensity values (can be for either microarray or RNA-Seq data).
	
	Sample IDs (listed in column header) must match the sample IDs used for the DNA methylation data (e.g. those listed in beta.table)
}
  \item{sample.file}{
	Tab-delimited text file providing group attributions for all samples considered for analysis.
	
	Only used if plot.scatter=TRUE
} 
  \item{cor.cutoff}{
	The minimum negative correlation coefficient to define a differentially expressed gene.
}
  \item{cor.pvalue.cutoff}{
	Maximum p-value allowed to define a gene as differentially expressed.
}
  \item{cor.fdr.cutoff}{
	Maximum False Discovery Rate (FDR) allowed to define a gene as differentially expressed.
}
 \item{plot.scatter}{
	A logical value: Create scatter plot for genes with a significant negative correlation?
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'txt' for tab-delimited text file
}
}

\examples{
library("COHCAP")

dir <- system.file("extdata", package="COHCAP")
beta.file <- file.path(dir,"GSE42308_truncated.txt")
sample.file <- file.path(dir,"sample_GSE42308.txt")
project.folder <- getwd()
expression.file <- file.path(dir,"expression-Average_by_Island_truncated.txt")
project.name <- "450k_avg_by_island_test"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name,
				project.folder, ref="parental")
filtered.islands <- COHCAP.avg.by.island(sample.file, filtered.sites,
				beta.table, project.name, project.folder, ref="parental")
COHCAP.integrate.avg.by.island(filtered.islands, project.name,
				project.folder, expression.file, sample.file)
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}