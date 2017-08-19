\name{COHCAP.BSSeq.preprocess}
\alias{COHCAP.BSSeq.preprocess}

\title{
Preprocessing for Targeted BS-Seq data
}
\description{
Creates custom annotation file as well as COHCAP input file (for COHCAP.annotate).

This function is not necessary for Illumina methylation array analysis.

Output files will be created in specified locations
}
\usage{
COHCAP.BSSeq.preprocess(methyl.folder=getwd(),
				cohcap.inputfile = file.path(getwd(),"BS_Seq_combined.txt"),
				gene.table = file.path(getwd(),"GENCODE_Genes.bed"),
				targeted.regions = file.path(getwd(),"UCSC_CpG_Islands.bed"),
				annotation.file = file.path(getwd(),"COHCAP.targeted.BSSeq.anno.txt"),
				shore.length=2000)
}
\arguments{
  \item{methyl.folder}{
	Folder containing .bed files created using genome_methylation_bismark2bedGraph_v3.pl (following Bismark alignment)
}
  \item{cohcap.inputfile}{
	Output file containing a tab-delimited table of percentage methylation values.
	
	This table will compatible with the custom annotation file created by this function (annotation.file)
}
  \item{gene.table}{
	.bed file containing gene names and coordinates
}
  \item{targeted.regions}{
	.bed file containing regions selected for targeted BS-Seq
}
  \item{annotation.file}{
	Custom annotation file providing gene and targeted region mappings for CpG sites specifically covered in your Bismark alignment
}
  \item{shore.length}{
	Length of shores considered to be part of the CpG island (in bp upstream and downstream of targeted region coordinates)
}
}
\value{
	This function creates two tab-delimited text files.
	
	One is to be used to as a custom annotation file (annotation.file).
	
	The other is used to create an appropriate input file for COHCAP (cohcap.inputfile).
	
	This function will likely take several hours to run.  However, it only needs to be run once.
}

\examples{
library("COHCAP")

dir = system.file("extdata", package="COHCAP")
bed.folder = file.path(dir,"BSSeq")
gene.table = file.path(dir,"GENCODE_Genes_truncated.bed")
targeted.regions = file.path(dir,"UCSC_CpG_Islands_truncated.bed")

output.folder = getwd()
annotation.file = file.path(output.folder,"COHCAP.targeted.BSSeq.anno.txt")
cohcap.inputfile = file.path(output.folder,"BS_Seq_combined.txt")

COHCAP.BSSeq.preprocess(bed.folder, cohcap.inputfile, gene.table,
				targeted.regions,annotation.file)
}

\seealso{
Useful Example Files: http://sourceforge.net/projects/cohcap/files/COHCAP_BSSEQ_anno.zip/download
*Default settings utilize these files (in current working directory)
*These files were created using the UCSC Genome Browser (build hg19)

Raw Data for Demo Dataset: http://www.ncbi.nlm.nih.gov/sra/SRX084504
Full, Formatted Demo Dataset (in standalone package): http://sourceforge.net/projects/cohcap/

Bismark: http://www.bioinformatics.babraham.ac.uk/projects/bismark/

UCSC Genome Browser: http://genome.ucsc.edu/cgi-bin/hgTables?command=start
}