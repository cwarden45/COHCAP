\name{COHCAP.avg.by.site}
\alias{COHCAP.avg.by.site}

\title{
CpG Island Differential Methylation Analysis (Average by Site Workflow).
}
\description{
Provides statistics for CpG islands as well as a list of differentially methylated sites.  CpG Island statistics are calculated by averaging beta values among samples per site
and comparing the average beta values across groups (considering the pairing between sites).

List of differentially methylated islands will be created in the "CpG_Island" folder.  Table of statistics for all CpG islands will be created in the "Raw_Data" folder.
}
\usage{
COHCAP.avg.by.site(site.table, project.name, project.folder,
				methyl.cutoff=0.7, unmethyl.cutoff = 0.3,
				delta.beta.cutoff = 0.2, pvalue.cutoff=0.05,
				fdr.cutoff=0.05, num.groups=2, num.sites=4,
				max.cluster.dist = NULL,
				output.format = "txt")
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
	Folder for COHCAP output files
}
  \item{methyl.cutoff}{
	Minimum beta or percentage methylation value to be used to define a methylated CpG site.
	Default is 0.7 (used for beta values), which would correspond to 70% methylation.
	Used for either 1-group or 2-group comparison.
}
  \item{unmethyl.cutoff}{
	Minimum beta or percentage methylation value to be used to define an unmethylated CpG site.
	Default is 0.3 (used for beta values), which would correspond to 30% methylation.
	Used for either 1-group or 2-group comparison.
}
  \item{delta.beta.cutoff}{
	The minimum absolute value for delta-beta values (mean treatment beta - mean reference beta) to define a differentially methylated CpG site.
	Only used for 2-group comparison.
}
  \item{pvalue.cutoff}{
	Maximum p-value allowed to define a CpG island as differentially methylated.
}
  \item{fdr.cutoff}{
	Maximum False Discovery Rate (FDR) allowed to define CpG island as differentially methylated.
}
 \item{num.groups}{
	Number of groups described in sample description file.  COHCAP algorithm differs when analysing 1-group, 2-group, or >2-group comparisons.
	COHCAP cannot analyze continuous phenotypes using the Average-by-Site workflow.
}
 \item{num.sites}{
	Minimum number of differentially methylated sites to define a differentially methylated CpG island.
}
  \item{max.cluster.dist}{
	Update annotations by running de-novo clustering within each set of provided annotations.
	This is the maximum distance (in bp) between filtered sites with a consistent delta-beta sign.
	This can produce more than one cluster per pre-existing annotation.
	Set to NULL by default, to run standard COHCAP algorithm.
	If you would like to test this function, I would recommend a value between 50 and 500 bp.
	Only valid for 2-group or continuous comparison.
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'csv' for comma-separated file, or 'txt' for tab-delimited text file.
}
}

\value{
	Data frame of average beta (or percentage methylation) statistics and/or p-value / false discovery rate statistics (per CpG island).
	
	The content of the data frame depends upon the number of groups specified for analysis.  All workflows provide p-values and FDR values.
	1 and 2 group comparisons provide counts for methylated and unmethylated sites as well as an overall methylation status per island.  >2 group comparisons only provide counts for the total number of differentially methylated sites.
	
	This data frame can be used for integration analysis.
}

\examples{
library("COHCAP")

dir = system.file("extdata", package="COHCAP")
beta.file = file.path(dir,"GSE42308_truncated.txt")
sample.file = file.path(dir,"sample_GSE42308.txt")
project.folder = tempdir()#you may want to use getwd() or specify another folder
project.name = "450k_avg_by_site_test"

beta.table = COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
filtered.sites = COHCAP.site(sample.file, beta.table, project.name,
				project.folder, ref="parental")
filtered.islands = COHCAP.avg.by.site(filtered.sites, project.name,
				project.folder)
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}
