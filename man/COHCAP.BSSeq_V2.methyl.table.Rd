\name{COHCAP.BSSeq_V2.methyl.table}
\alias{COHCAP.BSSeq_V2.methyl.table}

\title{
Create Percent Methylation Table from Bismark Coverage files for Targeted BS-Seq data
}
\description{
Creates percent methylation for COHCAP input (for COHCAP.annotate).

This function is not necessary for Illumina methylation array analysis.

NOTE: This does not create an annotation file.
There is a script that can create an annotation file using the output of this function
and an Ensembl .gtf (for hg19): *inst/extdata/Perl/downloaded_Ensembl_annotation.pl*

However, it is possible that the downloaded_Ensembl_annotation.pl script may need to be edited for other builds / organisms.
}
\usage{
COHCAP.BSSeq_V2.methyl.table(cov.files, sampleIDs, percent.table.output,
				min.cov = 10, min.percent.observed = 0.75,
				chr.index = 1, pos.index = 2, percent.index = 4,
				methylated.count.index = 5, unmethylated.count.index = 6,
				read.gz=FALSE)
}
\arguments{
  \item{cov.files}{
	Array of Bismark coverage files created using 'bismark_methylation_extractor' (following Bismark alignment)
	
	Files can be compressed.  However, to double-check files, please set read.gz to TRUE if using the original, compressed files.
}
  \item{sampleIDs}{
	Array of sample IDs to be used as column heads in percent methylation table.
	
	Paired determined by order (which needs to be the same) in 'cov.files' and 'sampleIDs'
}
  \item{percent.table.output}{
	Tab-delimited text output file for this function.
	
	This can be used as the input for COHCAP.annotate(), but an additional annotation file is also required for that function.
}
  \item{min.cov}{
	Minimum coverage to list methylation at site for a given samples.
	
	Default is 10 reads.
	
	If site is present at this coverage in other samples, this value will be represented as an "NA"
}
  \item{min.percent.observed}{
	Minimum percent of samples with at least 'min.cov' at a particular site.
	
	Default is 0.75 (75\%)
}
  \item{chr.index}{
	1-based index for chromosome in input file (element within 'cov.files' array).
	
	So, if this value is changed, function can potentially be applied to any table with all necessary values.
}
  \item{pos.index}{
	1-based index for chromosome position in input file (element within 'cov.files' array).
	
	So, if this value is changed, function can potentially be applied to any table with all necessary values.
}
  \item{percent.index}{
	1-based index for percent methylation in input file (element within 'cov.files' array).
	
	So, if this value is changed, function can potentially be applied to any table with all necessary values.
	
	Not strictly necessary if methylated and unmethylated counts also present, but this is currently required to be defined.
}
  \item{methylated.count.index}{
	1-based index for methylated counts in input file (element within 'cov.files' array).
	
	So, if this value is changed, function can potentially be applied to any table with all necessary values.
}
  \item{unmethylated.count.index}{
	1-based index for unmethylated counts in input file (element within 'cov.files' array).
	
	So, if this value is changed, function can potentially be applied to any table with all necessary values.
}
  \item{read.gz}{
	Logical: are all files in 'cov.files' compressed?
	
	One value, so compression status must be the same for all files.
	
	May not be strictly necessary, but probably good to check input file type.
}
}
\value{
	This function creates a tab-delimited text containing methylation values
	for samples present in at least 'min.percent.observed' samples.
	
	This table would be comparable to the beta input file for COHCAP.annotate().
	
	The 1st column from this table could also be used to create an annotation file.
}

\examples{
library("COHCAP")

dir = system.file("extdata/BSSeq", package="COHCAP")
rep1 = file.path(dir,"SRR096437_truncated.bismark.cov")
rep2 = file.path(dir,"SRR096438_truncated.bismark.cov")

cov.files = c(rep1,rep2)
sampleIDs = c("MCF7.Rep1","MCF7.Rep2")
percent.table="BS_Seq_combined_V2.txt"

COHCAP.BSSeq_V2.methyl.table(cov.files, sampleIDs, percent.table, read.gz=FALSE, min.percent.observed = 0.45)
}

\seealso{
Bismark: http://www.bioinformatics.babraham.ac.uk/projects/bismark/

Raw Data for Demo Dataset (different samples within study described in COHCAP paper):
https://www.ncbi.nlm.nih.gov/sra/SRR096437
https://www.ncbi.nlm.nih.gov/sra/SRR096438

BioPerl: https://bioperl.org/

Ensembl: https://ensembl.org

For hg19, there are steps to find and download the Ensembl .gtf are within:
			inst/extdata/Perl/downloaded_Ensembl_annotation.pl
}