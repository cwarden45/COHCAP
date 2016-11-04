\name{COHCAP.site}
\alias{COHCAP.site}

\title{
CpG Site Differential Methylation Analysis
}
\description{
Provides statistics for CpG sites as well as a list of differentially methylated sites.  Can also provide .wig files for visualization in IGV, UCSC Genome Browser, etc.

List of differentially methylated sites and .wig files will be created in the "CpG_Site" folder.  Table of statistics for all CpG sites will be created in the "Raw_Data" folder.
}
\usage{
COHCAP.site(sample.file, beta.table, project.name, project.folder,
				methyl.cutoff=0.7, unmethyl.cutoff = 0.3, paired=FALSE,
				delta.beta.cutoff = 0.2, pvalue.cutoff=0.05,
				fdr.cutoff=0.05, ref="none", num.groups=2,
				create.wig = "avg",plot.heatmap=TRUE, output.format = "xls")
}
\arguments{
  \item{sample.file}{
	Tab-delimited text file providing group attributions for all samples considered for analysis.
} 
 \item{beta.table}{
	Data frame with CpG sites in columns (with DNA methylation represented as beta values or percentage methylation),
	samples in columns, and CpG site annotations are included (in columns 2-5).
	
	The COHCAP.annotate function automatically creates this file.
}
  \item{project.name}{
	Name for COHCAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for COHCAP output files
}
  \item{methyl.cutoff}{
	Minimum beta or percentage methylation value to be used to define a methylated CpG site.
	Default is 0.7 (used for beta values), which would correspond to 70 percent methylation.
	Used for either 1-group or 2-group comparison.
}
  \item{unmethyl.cutoff}{
	Minimum beta or percentage methylation value to be used to define an unmethylated CpG site.
	Default is 0.3 (used for beta values), which would correspond to 30 percent methylation.
	Used for either 1-group or 2-group comparison.
}
  \item{delta.beta.cutoff}{
	The minimum absolute value for delta-beta values (mean treatment beta - mean reference beta) to define a differentially methylated CpG site.
	Only used for 2-group comparison.
}
  \item{pvalue.cutoff}{
	Maximum p-value allowed to define a site as differentially methylated.
	Used only for comparisons with at least 2 groups (with 3 replicates per group)
}
  \item{fdr.cutoff}{
	Maximum False Discovery Rate (FDR) allowed to define a site as differentially methylated.
	Used only for comparisons with at least 2 groups (with 3 replicates per group)
}
  \item{ref}{
	Reference group used to define baseline methylation levels.
	Set to "continuous" for a continuous primary variable.
	Otherwise, only used for 2-group comparison.
}
 \item{num.groups}{
	Number of groups described in sample description file.  COHCAP algorithm differs when analysing 1-group, 2-group, or >2-group comparisons.
	COHCAP cannot currently analyze continuous phenotypes.
}
 \item{create.wig}{
	Set to "avg" to create average beta (per group) and delta-beta values.
	Set to "sample" to create .wig files for each sample.
	Set to "avg.and.sample" to create average, delta-beta, and per-sample .wig files.
	Seto to "none" to avoid creating .wig files
	In the standalone version of COHCAP, this was only an option when using the "Average by Site" workflow (because
	that was the only situation where the analysis method matched the visualization).
	.wig files are defined with respect to hg19 (for pre-defined annotation files) and can be visualized using IGV, UCSC Genome Browser, etc.
}
 \item{plot.heatmap}{
	Logical value: Should heatmap be created to visualize CpG site differential methylation?
	For best comparison to island heatmap, please use same parameters at site and island level.
}
 \item{paired}{
	A logical value: Is there any special pairing between samples in different groups?
	If so, the pairing variable must be specified in the 3rd column of the sample description file.
	Used for p-value calculation, so this only applies to comparisons with at least 2 groups.

	If you have a secondary continuous variable (like age), you can set paired to "continuous".
	COHCAP will then perform linear-regression analysis (converting primary categorical variable into continuous variable, if necessary)
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'csv' for comma-separated file, or 'txt' for tab-delimited text file.
}
}

\value{
	Data frame of average beta (or percentage methylation) statistics and/or p-value / false discovery rate statistics.
	
	The content of the data frame depends upon the number of groups specified for analysis (avg.beta only for 1-group; avg.beta, delta.beta, p-value,
	and FDR for 2-group; p-value and FDR only for >2 groups).
	
	This data frame is used for CpG island analysis.
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
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}