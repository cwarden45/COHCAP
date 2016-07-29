`COHCAP.denovo` =function(site.table, project.name, project.folder, min.sites = 4, max.dist = 500, output.format = "xls")
{
	site.folder<-file.path(project.folder,"CpG_Site")
	
	print(dim(site.table))
	site.table <- site.table[!is.na(site.table[[8]]), ]
	print(dim(site.table))
	
	cluster.chr = c()
	cluster.start = c()
	cluster.stop = c()
	cluster.num.sites = c()

	site.chrs = as.character(levels(as.factor(as.character(site.table$Chr))))
	
	for (site.chr in site.chrs){
		chr.table = site.table[as.character(site.table$Chr) == site.chr,]
		if (nrow(chr.table) >= min.sites){
			print(paste("Checking for clusters on chr",site.chr,sep=""))
			print(dim(chr.table))
			chr.table = chr.table[order(chr.table$Loc),]
			probe.pos = chr.table$Loc
			
			temp.region.chr = as.character(chr.table$Chr[1])
			temp.region.start = chr.table$Loc[1]
			temp.region.stop = chr.table$Loc[1]
			temp.delta.beta = chr.table[1,8]
			temp.site.count = 1

			for (i in 2:length(probe.pos)){
				test.pos = probe.pos[i]
				test.delta.beta = chr.table[i,8]
				sign.check = FALSE
				if((temp.delta.beta > 0) & (test.delta.beta > 0)){
					sign.check = TRUE
				}else if((temp.delta.beta < 0) & (test.delta.beta < 0)){
					sign.check = TRUE
				}
				
				if(sign.check){
					if(((test.pos - temp.region.stop) <= max.dist) & (temp.site.count == 1)){
						#start cluster
						temp.region.stop = test.pos				
						temp.site.count = 2
					}else if (temp.site.count > 1){
						if((test.pos - temp.region.stop) <= max.dist){
							#expand region
							temp.region.stop = test.pos	
							temp.site.count = temp.site.count + 1
						} else {
							#define stop
							if(temp.site.count >= min.sites){
								cluster.chr = c(cluster.chr, paste("chr",temp.region.chr,sep=""))
								cluster.start = c(cluster.start, temp.region.start)
								cluster.stop = c(cluster.stop, temp.region.stop)
								
								cluster.num.sites = c(cluster.num.sites, temp.site.count)
							}#end if(temp.site.count >= min.sites)

							temp.region.chr = as.character(chr.table$Chr[i])
							temp.region.start = chr.table$Loc[i]
							temp.region.stop = chr.table$Loc[i]
							temp.site.count = 1
						}#end else
					}else{
						#reset
						temp.region.chr = as.character(chr.table$Chr[i])
						temp.region.start = chr.table$Loc[i]
						temp.region.stop = chr.table$Loc[i]
						temp.site.count = 1					
					}#end else
				}else{
					#reset
					if (temp.site.count >= min.sites){
						cluster.chr = c(cluster.chr, paste("chr",temp.region.chr,sep=""))
						cluster.start = c(cluster.start, temp.region.start)
						cluster.stop = c(cluster.stop, temp.region.stop)
										
						cluster.num.sites = c(cluster.num.sites, temp.site.count)
					}#end if (temp.hits > 1)

					temp.region.chr = as.character(chr.table$Chr[i])
					temp.region.start = chr.table$Loc[i]
					temp.region.stop = chr.table$Loc[i]
					temp.site.count = 1	
				}#end else
				
				temp.delta.beta = chr.table[i,8]
			}#end for (i in 2:length(test.pos))

			if (temp.site.count >= min.sites){
				cluster.chr = c(cluster.chr, paste("chr",temp.region.chr,sep=""))
				cluster.start = c(cluster.start, temp.region.start)
				cluster.stop = c(cluster.stop, temp.region.stop)
								
				cluster.num.sites = c(cluster.num.sites, temp.site.count)
			}#end if (temp.hits > 1)
		}#end if (nrow(chr.table) > min.sites)
	}#end for (site.chr in site.chrs)
	cluster.table = data.frame(hg19.chr=cluster.chr, hg19.start=cluster.start, hg19.stop=cluster.stop, cluster.num.sites=cluster.num.sites)
	
	if(nrow(cluster.table) > 0){
		print(dim(cluster.table))
		if(output.format == 'xls'){
			xlsfile <- file.path(site.folder, paste(project.name,"_CpG_site_clusters.xlsx",sep=""))
			WriteXLS("cluster.table", ExcelFileName = xlsfile)
		} else if(output.format == 'txt'){
			txtfile <- file.path(site.folder, paste(project.name,"_CpG_site_clusters.txt",sep=""))
			write.table(cluster.table, file=txtfile, quote=F, row.names=F, sep="\t")
		} else {
			warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
		}
	}else{
		print("Sorry, no clusters were defined.")
	}
}#end def COHCAP.denovo