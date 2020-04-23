\name{COHCAP.qc}
\alias{COHCAP.qc}

\title{
DNA Methylation Quality Control Statistics
}
\description{
Provides descriptive statistics (median, top/bottom quartiles, mininum,maximum),
sample histograms, sample dendrogram, principal component analysis plot.

Output files will be created in the "QC" subfolder.
}
\usage{
COHCAP.qc(sample.file, beta.table, project.name, project.folder,
				plot.legend=TRUE, color.palette = c("red","blue",
				"green","orange","purple","cyan","pink","maroon",
				"yellow","grey","black",colors()))
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
  \item{plot.legend}{
	A logical value: Should legend be plotted within QC figures?
}
  \item{color.palette}{
	Colors for primary variable (specified in the second column of the sample file).
	
	Remember, COHCAP can only analyze discrete variables categoried with groups (preferably with replicates).
}
}

\examples{
library("COHCAP")

dir = system.file("extdata", package="COHCAP")
beta.file = file.path(dir,"GSE42308_truncated.txt")
sample.file = file.path(dir,"sample_GSE42308.txt")
project.folder = tempdir()
project.name = "450k_test"

beta.table = COHCAP.annotate(beta.file, project.name, project.folder,
				platform="450k-UCSC")
COHCAP.qc(sample.file, beta.table, project.name, project.folder)
}

\seealso{
COHCAP Discussion Group: http://sourceforge.net/p/cohcap/discussion/general/
}