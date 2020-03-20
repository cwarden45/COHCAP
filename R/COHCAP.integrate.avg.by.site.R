`COHCAP.integrate.avg.by.site` <-function (island.table, project.name, project.folder, expr.file, expr.pvalue=0.05, expr.fdr = 0.05, expr.fc = 1.5, output.format = "xls")
{
	integrate.folder<-file.path(project.folder,"Integrate")
	dir.create(integrate.folder, showWarnings=FALSE)
	
	data.folder<-file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	temp.mU.eD <- file.path(data.folder,"mUeD.txt")
	temp.mD.eU <- file.path(data.folder,"mDeU.txt")
	temp.methyl <- file.path(data.folder,"temp_methyl.txt")
	write.table(island.table, temp.methyl, quote=F, row.names=F, sep="\t")

	
	Perl.Path <- file.path(path.package("COHCAP"), "Perl")
	perl.script <- file.path(Perl.Path , "integrate_site.pl")
	cmd <- paste("perl \"",perl.script,"\" \"", expr.file,"\" \"",temp.methyl,"\" \"",temp.mD.eU,"\" \"",temp.mU.eD,"\" ",expr.pvalue," ", expr.fdr," ",expr.fc, sep="")
	res <- system(cmd, intern=TRUE, wait=TRUE)
	message(res)
	
	mU.eD <- read.table(temp.mU.eD, header=T, sep="\t", stringsAsFactors=TRUE)
	if(output.format == 'xls'){
		xlsfile <- file.path(integrate.folder, paste(project.name,"_Methyl_Up_Expression_Down-Avg_by_Site.xlsx",sep=""))
		WriteXLS("mU.eD", ExcelFileName = xlsfile)
	} else if(output.format == 'csv') {
		txtfile <- file.path(integrate.folder, paste(project.name,"_Methyl_Up_Expression_Down-Avg_by_Site.csv",sep=""))
		write.table(mU.eD, file=txtfile, quote=F, row.names=F, sep=",")
	} else if(output.format == 'txt') {
		txtfile <- file.path(integrate.folder, paste(project.name,"_Methyl_Up_Expression_Down-Avg_by_Site.txt",sep=""))
		write.table(mU.eD, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	mD.eU <- read.table(temp.mD.eU, header=T, sep="\t", stringsAsFactors=TRUE)
	if(output.format == 'xls'){
		xlsfile <- file.path(integrate.folder, paste(project.name,"_Methyl_Down_Expression_Up-Avg_by_Site.xlsx",sep=""))
		WriteXLS("mD.eU", ExcelFileName = xlsfile)
	} else if(output.format == 'csv'){
		txtfile <- file.path(integrate.folder, paste(project.name,"_Methyl_Down_Expression_Up-Avg_by_Site.csv",sep=""))
		write.table(mD.eU, file=txtfile, quote=F, row.names=F, sep=",")
	} else if(output.format == 'txt'){
		txtfile <- file.path(integrate.folder, paste(project.name,"_Methyl_Down_Expression_Up-Avg_by_Site.txt",sep=""))
		write.table(mD.eU, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	unlink(temp.mU.eD)
	unlink(temp.mD.eU)
	unlink(temp.methyl)
}#end def RNA.deg
