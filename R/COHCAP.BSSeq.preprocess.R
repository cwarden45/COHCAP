`COHCAP.BSSeq.preprocess` <-function (methyl.folder=getwd(), cohcap.inputfile = file.path(getwd(),"BS_Seq_combined.txt"), gene.table = file.path(getwd(),"GENCODE_Genes.bed"), targeted.regions = file.path(getwd(),"UCSC_CpG_Islands.bed"), annotation.file = file.path(getwd(),"COHCAP.targeted.BSSeq.anno.txt"), shore.length=2000)
{
	Perl.Path <- file.path(path.package("COHCAP"), "Perl")
	perl.script <- file.path(Perl.Path , "find_bismark_CpG_sites.pl")
	cmd <- paste("perl \"",perl.script,"\" \"", methyl.folder,"\" \"", cohcap.inputfile,"\" \"", gene.table,"\" \"", targeted.regions,"\" \"", annotation.file,"\" ", shore.length, sep="")
	res <- system(cmd, intern=TRUE, wait=TRUE)
	message(res)
}#end def BSSeq.anno