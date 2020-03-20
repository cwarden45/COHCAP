`COHCAP.reformatFinalReport` = function(FinalReport, beta.file, renaming.file=NULL, detection.pvalue.cutoff=0.01)
{
	input.table = read.table(FinalReport, header=T, sep="\t", skip=8, stringsAsFactors=TRUE)
	probeID = input.table[,1]
	#if(delete.last.column){
	#	input.table = input.table[,-ncol(input.table)]
	#}
	#beta.mat = round(input.table[,seq(2,ncol(input.table),2)], digits=2)
	#pvalue.mat = input.table[,seq(3,ncol(input.table),2)]
	
	beta.mat = round(input.table[,grep("AVG_Beta",colnames(input.table))], digits=2)
	pvalue.mat = input.table[,grep("Detection.Pval",colnames(input.table))]
	
	rm(input.table)
	
	beta.names = colnames(beta.mat)
	#beta.names = beta.names[grep("AVG_Beta",beta.names)]
	beta.names = gsub(".AVG_Beta","",beta.names)
	pval.names = colnames(pvalue.mat)
	#pval.names = pval.names[grep("Detection.Pval",pval.names)]
	pval.names = gsub(".Detection.Pval","",pval.names)
	
	if((length(beta.names)!= length(pval.names))&any(beta.names != pval.names)){
		print("Discordance between Beta matrix and Detection P-value matrix!")
		print(paste("Beta names: ",beta.names,sep=""))
		print(paste("P-value names: ",pval.names,sep=""))
		print("If only last p-value is missing, consider setting 'delete.last.column=FALSE'")
		stop()
	}else{
		filtered.beta.mat = beta.mat
		filtered.beta.mat[pvalue.mat > detection.pvalue.cutoff] = NA
		colnames(filtered.beta.mat)=beta.names
		rm(beta.mat)
		rm(pvalue.mat)
		
		if(!is.null(renaming.file)){
			label.table = read.table(renaming.file,header=T, sep="\t", stringsAsFactors=TRUE)
			chipID = as.character(label.table[,1])
			if(length(grep("^\\d",chipID)) > 0){
				print("Adding X to beginning of samples:")
				print(chipID)
				chipID[grep("^\\d",chipID)]=paste("X",chipID[grep("^\\d",chipID)],sep="")
				print(chipID)
			}
			
			new.labels = label.table[,2]
			new.labels = new.labels[match(beta.names,chipID)]
			names(new.labels)=beta.names
			if(length(new.labels[is.na(new.labels)]) != 0){
				print("Sample mapping issue:")
				print(new.labels)
				stop()
			}else{
				colnames(filtered.beta.mat)=new.labels
			}
		}
		
		filtered.beta.mat = data.frame(SiteID = probeID, filtered.beta.mat)
		write.table(filtered.beta.mat, beta.file, quote=F, row.names=F, sep="\t")
	}#end else
}#end def COHCAP.reformatFinalReport
