COHCAP.BSSeq_V2.methyl.table = function(cov.files, sampleIDs, percent.table.output,
										min.cov = 10, min.percent.observed = 0.75,
										chr.index = 1, pos.index = 2, percent.index = 4,
										methylated.count.index = 5, unmethylated.count.index = 6,
										read.gz=FALSE){
										
	##just create table with this function : leave annotation separate
	
	common.sites=c()
	
	for (i in 1:length(cov.files)){
		print(paste(sampleIDs[i],cov.files[i],sep=" --> "))

		input.file = cov.files[i]
		
		gz.result = input.file[grep(".gz$",input.file)]
		if((length(gz.result)==1) & (read.gz==FALSE)){
			print(paste(".gz file extension for ",input.file,sep=""))
			print("...but 'read.gz=FALSE'")
			print("Please double-check files or set 'read.gz=TRUE'")
			stop()
		}#end if((length(gz.result)==1) & (read.gz==FALSE))
			
		#This can work without the extra gzfile() function, but seems good to have
		if(read.gz){
			print("Reading compressed coverage file...")
			cov.table = read.table(gzfile(cov.files[i]), header=F, sep="\t", stringsAsFactors=TRUE)
		}else{
			print("Reading uncompressed coverage file...")
			cov.table = read.table(cov.files[i], header=F, sep="\t", stringsAsFactors=TRUE)
		}#end else
		
		temp.total.cov = cov.table[,methylated.count.index] + cov.table[,unmethylated.count.index]
		
		cov.table = cov.table[temp.total.cov >= min.cov,]
		temp.siteID = paste(as.character(cov.table[,chr.index]),as.character(cov.table[,pos.index]),sep=":")
		
		if(i == 1){
			common.sites = temp.siteID
			percent.table = data.frame(V4=cov.table[,percent.index])
			colnames(percent.table)=as.character(sampleIDs[i])
			print(head(percent.table))
		}else{
			prev.common.sites = common.sites
			common.sites = union(common.sites, temp.siteID)
			prev.ids = colnames(percent.table)
			print(prev.ids)
			percent.table = percent.table[match(common.sites,prev.common.sites),]
			new.percent =cov.table[,percent.index]
			new.percent=new.percent[match(common.sites, temp.siteID)]
			
			percent.table = data.frame(percent.table, new.percent)
			col.ids =c(prev.ids,as.character(sampleIDs[i]))
			colnames(percent.table) = col.ids
			print(head(percent.table))
		}#end else
		
	}#end for (i in 1:length(cov.files))
	
	count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
	}#end count.covered
	site.covered.by.sample = apply(percent.table, 1, count.covered)
	
	print(dim(percent.table))
	common.sites = common.sites[site.covered.by.sample > min.percent.observed * ncol(percent.table)]
	percent.table = percent.table[site.covered.by.sample > min.percent.observed * ncol(percent.table),]
	percent.table = round(percent.table / 100, digits=2)
	print(dim(percent.table))
	percent.table = data.frame(SiteID=common.sites, percent.table)
	write.table(percent.table, percent.table.output, quote=F, sep="\t", row.names=F)
}#end def COHCAP.BSSeq_V2.methyl.table