\name{COHCAP.avg.by.island}
\alias{COHCAP.avg.by.island}

\title{
CpG Island Differential Methylation Analysis (Average by Island Workflow).
}
\description{
Provides statistics for CpG islands as well as a list of differentially methylated sites.  CpG Island statistics are calculated by averaging beta values among samples per site
and comparing the average beta values across groups (considering the pairing between sites).

List of differentially methylated islands will be created in the "CpG_Island" folder.  Table of statistics for all CpG islands will be created in the "Raw_Data" folder.
}
\usage{
COHCAP.avg.by.island(sample.file, site.table, beta.table, project.name, project.folder,
				methyl.cutoff=0.7, unmethyl.cutoff = 0.3, delta.beta.cutoff = 0.2,
				pvalue.cutoff=0.05, fdr.cutoff=0.05,
				num.groups=2, num.sites=4, plot.box=TRUE, plot.heatmap=TRUE,
				paired=FALSE, ref="none", lower.cont.quantile=0, upper.cont.quantile=1,
				max.cluster.dist = NULL, ttest.sub="none",
				output.format = "xls", gene.centric=TRUE)
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
 \item{site.table}{
	Data frame with differentially methylated CpG site statistics (one row per CpG site) and CpG site annotations (in columns 2-5).
	
	The COHCAP.site function automatically creates this file.
}
  \item{project.name}{
	Name for COHCAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for COHCAP output files
}
  \item{methyl.cutoff}{
	Minimum beta or percentage methylation value to be used to define a methylated CpG island.
	Default is 0.7 (used for beta values), which would correspond to 70% methylation.
	Used for either 1-group or 2-group comparison.
}
  \item{unmethyl.cutoff}{
	Minimum beta or percentage methylation value to be used to define an unmethylated CpG island.
	Default is 0.3 (used for beta values), which would correspond to 30% methylation.
	Used for either 1-group or 2-group comparison.
}
  \item{max.cluster.dist}{
	Update annotations by running de-novo clustering within each set of provided annotations.
	This is the maximum distance (in bp) between filtered sites with a consistent delta-beta sign.
	This can produce more than one cluster per pre-existing annotation.
	Set to NULL by default, to run standard COHCAP algorithm.
	If you would like to test this function, I would recommend a value between 50 and 500 bp.
	Only valid for 2-group or continuous comparison.
}
  \item{delta.beta.cutoff}{
	The minimum absolute value for delta-beta values (mean treatment beta - mean reference beta) to define a differentially methylated CpG island.
	Only used for 2-group comparison (and continuous comparison, where delta beta is max - min beta, with sign based upon correlation coefficient).
}
  \item{pvalue.cutoff}{
	Maximum p-value allowed to define an island as differentially methylated.
	Used only for comparisons with at least 2 groups (with 3 replicates per group)
}
  \item{fdr.cutoff}{
	Maximum False Discovery Rate (FDR) allowed to define an island as differentially methylated.
	Used only for comparisons with at least 2 groups (with 3 replicates per group)
}
 \item{num.groups}{
	Number of groups described in sample description file.  COHCAP algorithm differs when analysing 1-group, 2-group, or >2-group comparisons.
	COHCAP cannot currently analyze continuous phenotypes.
}
 \item{num.sites}{
	Minimum number of differentially methylated sites to define a differentially methylated CpG island.
}
  \item{ref}{
	Reference group used to define baseline methylation levels.
	Set to "continuous" for a continuous primary variable.
	Otherwise, only used for 2-group comparison.
}
 \item{plot.box}{
	Logical value: Should box-plots be created to visualize CpG island differential methylation?
	If using a continuous primary variable, line plots are provided instead of box plots.
}
 \item{plot.heatmap}{
	Logical value: Should heatmap be created to visualize CpG island differential methylation?
}
  \item{lower.cont.quantile}{
	For continuous analysis, what beta quantile should be the lower threshold for calculating delta-beta values?
	Default = 0 (minimum)
}
  \item{upper.cont.quantile}{
	For continuous analysis, what beta quantile should be the upper threshold for calculating delta-beta values?
	Default = 1 (maximum)
}
 \item{paired}{
	A logical value: Is there any special pairing between samples in different groups?
	
	If so, the pairing variable must be specified in the 3rd column of the sample description file.
	Used for p-value calculation, so this only applies to comparisons with at least 2 groups.

	If you have a secondary continuous variable (like age), you can set paired to "continuous".
	COHCAP will then perform linear-regression analysis (converting primary categorical variable into continuous variable, if necessary)
}
  \item{ttest.sub}{
	In original workflow, use an alternative test instead of t-test (for
		1-variable, 2-group comparison).
	
	Can be 'none' or 'ANOVA'.  Might be helpful when SD is 0.
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'csv' for comma-separated file, or 'txt' for tab-delimited text file
}
  \item{gene.centric}{
	Should CpG islands not mapped to genes be ignored?
	Default: TRUE (Recommended setting for integration with gene expession data)
}
}

\value{
	List (island.list) with 2 data frames to be used for integration analysis:
	
	beta.table = data frame of average beta (or percentage methylation) values across differentially methylated sites within a differentially methylated CpG island
	filtered.island.stats = differential methylation results for CpG Islands
}

\examples{
library("COHCAP")

dir <- system.file("extdata", package="COHCAP")
beta.file <- file.path(dir,"GSE42308_truncated.txt")
sample.file <- file.path(dir,"sample_GSE42308.txt")
project.folder <- getwd()
project.name <- "450k_avg_by_island_test"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name,
				project.folder, ref="parental")
island.list <- COHCAP.avg.by.island(sample.file, filtered.sites,
				beta.table, project.name, project.folder, ref="parental")
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}