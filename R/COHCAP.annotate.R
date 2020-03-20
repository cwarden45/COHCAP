`COHCAP.annotate` <-function (beta.file, project.name, project.folder, platform, annotation.file = NULL, output.format = "xls")
{
	data.folder<-file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	beta.table <- read.table(beta.file, header=T, sep="\t", stringsAsFactors=TRUE)
	print(dim(beta.table))
	data.mat <- beta.table[,2:ncol(beta.table)]
	
	if(platform == "450k-UCSC"){
		data(COHCAP.450k.UCSC)
		annotation.table <- COHCAP.450k.UCSC
	} else if(platform == "450k-HMM"){
		data(COHCAP.450k.HMM)
		annotation.table <- COHCAP.450k.HMM
	} else if(platform == "27k"){
		data(COHCAP.27k)
		annotation.table <- COHCAP.27k
	} else if((platform == "custom") & (length(annotation.file) == 1)) {
		print(paste("Using custom island/gene annotations from : ",annotation.file, sep=""))
		annotation.table <- read.table(annotation.file, sep="\t", header=T, stringsAsFactors=TRUE)
		print(dim(annotation.table))
	} else{
		stop("You need to provide a valid annotation platform and/or custom annotation file!")
	}
	
	matched.annotations <- annotation.table[match(beta.table$SiteID, annotation.table$SiteID, nomatch=0),]
	print(dim(matched.annotations))
	data.mat <- data.mat[match(matched.annotations$SiteID, beta.table$SiteID, nomatch=0),]
	if(nrow(data.mat) != nrow(matched.annotations))
		{
			stop("Annotations not present for all CpG sites.  Please check your annotation file!")
		}
	ann.mat <- cbind(matched.annotations,data.mat)
	
	print(dim(ann.mat))
	#print(head(ann.mat))
	if(output.format == 'xls'){
		xlsfile <- file.path(data.folder, paste(project.name,"_annotated_CpG_sites.xlsx",sep=""))
		WriteXLS("ann.mat", ExcelFileName = xlsfile)
	} else if(output.format == 'csv') {
		txtfile <- file.path(data.folder, paste(project.name,"_annotated_CpG_sites.csv",sep=""))
		write.table(ann.mat, file=txtfile, quote=F, row.names=F, sep=",")
	} else if(output.format == 'txt') {
		txtfile <- file.path(data.folder, paste(project.name,"_annotated_CpG_sites.txt",sep=""))
		write.table(ann.mat, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	return(ann.mat)
}#end def COHCAP.annotate