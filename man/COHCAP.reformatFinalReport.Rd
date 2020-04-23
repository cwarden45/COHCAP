\name{COHCAP.reformatFinalReport}
\alias{COHCAP.reformatFinalReport}

\title{
Prepare COHCAP beta file from Illumina GenomeStudio Final Report
}
\description{
Reformats FinalReport file into format for COHCAP.annotate() function.

Optionally allows user to re-name samples, using a table where the first column is the chipID and the 2nd column is the name that should appear in heatmaps, etc.
}
\usage{
COHCAP.reformatFinalReport(FinalReport, beta.file, renaming.file=NULL,
										detection.pvalue.cutoff=0.01)
}
\arguments{
  \item{FinalReport}{
	FinalReport from Illumina GenomeStudio.  "AVG_Beta" and "Detection Pval" columns should be exported for "Sample Methylation Profile".
	
	The goal of this function is to create a single table with beta values, with the Detection P-value being used to remove measurements with low intensities.
}
  \item{beta.file}{
	Table of beta / percentage methylation values.  CpG sites are represented in rows.
	
	Samples are represented in columns.
	
	Samples with high detection p-values will be censored as NA values.
}

  \item{renaming.file}{
	If you would like to provide more descriptive names for your samples, uprovide a table where the first column is the chipID and the 2nd column is the name that should appear in heatmaps, etc.

	Please make sure your new names don't begin with numbers, or contain certain special characters (such as dashes, spaces, etc) that will be reformated as periods in column names when imported into R.
}

  \item{detection.pvalue.cutoff}{
	Maximum p-value to keep beta values in reformated table.
}

}
\value{
	Function does not return value.  Instead, it writes the censored beta table for COHCAP input as a text file (specified as 'beta.file').
}

\examples{
library("COHCAP")

dir = system.file("extdata", package="COHCAP")
FinalReport = file.path(dir,"EPIC_DEMO.txt")
beta.file = "COHCAP_beta.txt"
sample.relabel.file = file.path(dir,"EPIC_DEMO_samples.txt")

COHCAP.reformatFinalReport(FinalReport, beta.file, renaming.file=sample.relabel.file)

project.folder = tempdir()
project.name = "EPIC_FinalReport_test"
annotation.file = file.path(dir,"EPIC_DEMO_mapping.txt")
beta.table = COHCAP.annotate(beta.file, project.name, project.folder,
				platform="custom",annotation.file=annotation.file,
				output.format = "txt")
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}