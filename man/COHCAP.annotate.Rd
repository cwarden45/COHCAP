\name{COHCAP.annotate}
\alias{COHCAP.annotate}

\title{
Annotation for CpG Sites
}
\description{
Provides annotations (chromosome, location, gene, and CpG island) for CpG sites from a specified annotation file.

Annotations for common platforms are provided with COHCAP (with respect to hg19).  Custom annotation files can also be provided.

Output files will be created in the "Raw_Data" subfolder.
}
\usage{
COHCAP.annotate(beta.file, project.name, project.folder,
				platform, annotation.file = NULL,
				output.format = "txt")
}
\arguments{
  \item{beta.file}{
	Table of beta / percentage methylation values.  CpG sites are represented in rows.
	
	Samples are represented in columns.
}
  \item{project.name}{
	Name for COHCAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for COHCAP output files
}
  \item{platform}{
	Annotation file to be used.  Enter "450k-UCSC" for UCSC CpG Islands for 450k array probes, "450k-HMM" for HMM CpG Islands for 450k array probes,
	"27k" for UCSC CpG Islands for 27k array probes.
	
	If none of these pre-defined annotations are acceptable, please enter "custom" for the platform and provide an annotation file.
}
  \item{annotation.file}{
	Annotation file to be used for a custom platform.  This variable is not used for common, pre-defined annotation files.
	
	The annotation file should be a tab-delimited text file with the header "SiteID	Chr	Loc	Gene	Island", with columns respectively specifying the CpG identifier (must match beta / percent
	methylation file), chromosome for CpG site, position for CpG site (preferably in hg19 coordinates), nearest gene mapping for CpG site, nearest CpG island mapping for CpG site.
}
  \item{output.format}{
	Format for output tables: 'xls' for Excel file, 'csv' for comma-separated file, or 'txt' for tab-delimited text file.
}
}
\value{
	Data frame of beta values (must be between 0 and 1) or percentage methylation values (must be between 0 and 100).
	
	Just like the input table, the first column specifies the SiteID, CpG sites are represented on rows,
	samples are represented in samples (starting with the 6th column).  Additionally, the 2nd column now specifies the CpG site chromosome,
	the 3rd column now specifies the CpG site position (in hg19 coordinates, for pre-defined annotation files), the 4th column lists the nearest gene mapping,
	and the 5th column lists the nearest CpG island mapping.
	
	This data frame is used for quality control and differential methylation analysis.
}

\examples{
library("COHCAP")

dir = system.file("extdata", package="COHCAP")
beta.file = file.path(dir,"GSE42308_truncated.txt")
project.folder = tempdir()
project.name = "450k_test"

beta.table = COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}
