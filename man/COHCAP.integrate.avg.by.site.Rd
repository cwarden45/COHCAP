\name{COHCAP.integrate.avg.by.site}
\alias{COHCAP.integrate.avg.by.site}

\title{
Integration with Gene Expression Data(Average by Site Workflow).
}
\description{
Provides lists of genes with an inverse CpG island methylation trend (Methylation Down, Expression Up and Methylation Up, Expression Down).

Lists of genes with negative expression trends will be created in the "Integrate" folder.

The "Average by Site" workflow requires that genes already have fold-change, p-value, and FDR values calculated.  There many tools available for this type of analysis (limma, DEseq2, edgeR, etc.)

This function will only work for 2-group comparisons.
}
\usage{
COHCAP.integrate.avg.by.site(island.table, project.name, project.folder,
				expr.file, expr.pvalue=0.05, expr.fdr = 0.05, expr.fc = 1.5,
				output.format = "txt")
}
\arguments{
 \item{island.table}{
	Data frame with CpG island statistics (one row per CpG island) for differentially methylated CpG islands.
	
	The COHCAP.avg.by.site function automatically creates this file.
}
  \item{project.name}{
	Name for COHCAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for COHCAP output files
}
  \item{expr.file}{
	Table of differential expression statistics.  Gene symbols must be in the first column, fold-change values must be in the second column,
	p-values must be in the third column, and false discovery rate (FDR) values must be in the fourth column.

	These statistics must be calculated outside of COHCAP.  Duplicate gene symbols are OK - statistics will be averaged among duplicate gene symbols.
}
  \item{expr.fc}{
	The minimum absolute value for fold-change values (treatment versus reference) to define gene as differentially expressed (from gene expression table).
	Only used for 2-group comparison.  Fold-change is expected to be on a linear scale.
}
  \item{expr.pvalue}{
	Maximum p-value allowed to define a gene as differentially expressed (from gene expression table).
}
  \item{expr.fdr}{
	Maximum False Discovery Rate (FDR) allowed to define a gene as differentially expressed (from gene expression table)
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'csv' for comma-separated file, or 'txt' for tab-delimited text file
}
}

\examples{
library("COHCAP")

dir = system.file("extdata", package="COHCAP")
beta.file = file.path(dir,"GSE42308_truncated.txt")
sample.file = file.path(dir,"sample_GSE42308.txt")
project.folder = tempdir()#you may want to use getwd() or specify another folder
expression.file = file.path(dir,"expression-Average_by_Site_truncated.txt")
project.name = "450k_avg_by_site_test"

beta.table = COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
filtered.sites = COHCAP.site(sample.file, beta.table, project.name,
				project.folder, ref="parental")
filtered.islands = COHCAP.avg.by.site(filtered.sites, project.name,
				project.folder)
COHCAP.integrate.avg.by.site(filtered.islands, project.name, project.folder,
				expression.file)
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}
