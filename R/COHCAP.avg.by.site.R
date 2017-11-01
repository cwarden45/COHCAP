custom.mean = function(arr)
{
	return(mean(as.numeric(arr), na.rm = T))
}#end def custom.mean

annova.pvalue = function(arr, grp.levels)
{
	#print(arr)
	#print(grp.levels)
	
	grp.no.na = as.factor(as.character(grp.levels[!is.na(arr)]))
	if(length(levels(grp.no.na)) >= 2)
		{
			test = aov(arr ~ grp.levels) 
			result = summary(test)[[1]][['Pr(>F)']][1]
			if(is.null(result))
				{
					return(1)
				}
			else
				{
					return(result)
				}
		}
	else
		{
			return(1)
		}
}#end def annova.pvalue


annova.2way.pvalue = function(arr, grp.levels, pairing.levels)
{
	#print(arr)
	#print(grp.levels)
	#print(pairing.levels)
	
	grp.no.na = as.factor(as.character(grp.levels[!is.na(arr)]))
	
	rep.flag = 1
	if((length(levels(grp.no.na)) >= 2) && (rep.flag == 1))
		{
			test = aov(arr ~ grp.levels + pairing.levels) 
			result = summary(test)[[1]][['Pr(>F)']][1]
			if(is.null(result))
				{
					return(1)
				}
			else
				{
					return(result)
				}
		}
	else
		{
			return(1)
		}
}#end def annova.pvalue

`COHCAP.avg.by.site` = function (site.table, project.name, project.folder, methyl.cutoff=0.7, unmethyl.cutoff = 0.3, delta.beta.cutoff = 0.2, pvalue.cutoff=0.05, fdr.cutoff=0.05, num.groups=2, num.sites=4, max.cluster.dist = NULL, output.format = "xls")
{
	island.folder=file.path(project.folder,"CpG_Island")
	dir.create(island.folder, showWarnings=FALSE)
	
	data.folder=file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	print("Checking CpG Site Stats Table")
	print(dim(site.table))
	site.table = site.table[!is.na(site.table[[5]]), ]
	print(dim(site.table))
	
	if(!is.null(max.cluster.dist)){
		if ((num.groups==2)|(ref=="continuous")){
			print("Looking for clusters within provided CpG Island Annotations")
			island.clusters = rep(NA,nrow(site.table))
			
			cpg.islands = levels(as.factor(as.character(site.table[[5]])))

			for (i in 1:length(cpg.islands))
				{
					island = cpg.islands[i]
					cpg.sites = as.character(site.table[site.table[[5]] == island,1])
					ann.mat = site.table[match(cpg.sites, as.character(site.table[,1]),nomatch=0),]
					if(nrow(ann.mat) >= num.sites)
						{
							#print(island)
							ann.mat = ann.mat[order(ann.mat$Loc),]
							probe.pos = ann.mat$Loc
						
							temp.cluster.probes = c(as.character(ann.mat$SiteID[1]))
							temp.region.chr = as.character(ann.mat$Chr[1])
							temp.region.start = ann.mat$Loc[1]
							temp.region.stop = ann.mat$Loc[1]
							temp.delta.beta = ann.mat[1,8]
							temp.site.count = 1

							for (i in 2:length(probe.pos)){
								test.pos = probe.pos[i]
								test.delta.beta = ann.mat[i,8]
								sign.check = FALSE
								if((temp.delta.beta > 0) & (test.delta.beta > 0)){
									sign.check = TRUE
								}else if((temp.delta.beta < 0) & (test.delta.beta < 0)){
									sign.check = TRUE
								}
								
								if(sign.check){
									if(((test.pos - temp.region.stop) <= max.cluster.dist) & (temp.site.count == 1)){
										#start cluster
										temp.region.stop = test.pos				
										temp.site.count = 2
										temp.cluster.probes = c(temp.cluster.probes, as.character(ann.mat$SiteID[i]))
									}else if (temp.site.count > 1){
										if((test.pos - temp.region.stop) <= max.cluster.dist){
											#expand region
											temp.region.stop = test.pos	
											temp.site.count = temp.site.count + 1
											temp.cluster.probes = c(temp.cluster.probes, as.character(ann.mat$SiteID[i]))
										} else {
											#define stop
											if(temp.site.count >= num.sites){
												cluster.chr = paste("chr",temp.region.chr,sep="")
												cluster.start = temp.region.start
												cluster.stop = temp.region.stop
												
												new.island = paste(cluster.chr,":",cluster.start,"-",cluster.stop,sep="")
												island.clusters[match(temp.cluster.probes,as.character(site.table$SiteID), nomatch=0)]=new.island
											}#end if(temp.site.count >= num.sites)

											temp.cluster.probes = c(as.character(ann.mat$SiteID[i]))
											temp.region.chr = as.character(ann.mat$Chr[i])
											temp.region.start = ann.mat$Loc[i]
											temp.region.stop = ann.mat$Loc[i]
											temp.site.count = 1
										}#end else
									}else{
										#reset
										temp.cluster.probes = c(as.character(ann.mat$SiteID[i]))
										temp.region.chr = as.character(ann.mat$Chr[i])
										temp.region.start = ann.mat$Loc[i]
										temp.region.stop = ann.mat$Loc[i]
										temp.site.count = 1					
									}#end else
								}else{
									#reset
									if (temp.site.count >= num.sites){
										cluster.chr = paste("chr",temp.region.chr,sep="")
										cluster.start = temp.region.start
										cluster.stop = temp.region.stop
														
										new.island = paste(cluster.chr,":",cluster.start,"-",cluster.stop,sep="")
										island.clusters[match(temp.cluster.probes,as.character(site.table$SiteID), nomatch=0)]=new.island
									}#end if (temp.site.count >= num.sites)
									
									temp.cluster.probes = c(as.character(ann.mat$SiteID[i]))
									temp.region.chr = as.character(ann.mat$Chr[i])
									temp.region.start = ann.mat$Loc[i]
									temp.region.stop = ann.mat$Loc[i]
									temp.site.count = 1	
								}#end else
								
								temp.delta.beta = ann.mat[i,8]
							}#end for (i in 2:length(test.pos))

							if (temp.site.count >= num.sites){
								cluster.chr = paste("chr",temp.region.chr,sep="")
								cluster.start = temp.region.start
								cluster.stop =  temp.region.stop
												
								new.island = paste(cluster.chr,":",cluster.start,"-",cluster.stop,sep="")
								island.clusters[match(temp.cluster.probes,as.character(site.table$SiteID), nomatch=0)]=new.island
							}#end if (temp.hits > 1)
						}#end if(nrow(data.mat) >= num.sites)
				}#end for (i in 1:length(cpg.islands))
			#print(island.clusters)
			
			mapping.table = data.frame(SiteID = as.character(site.table$SiteID), Chr =as.character(site.table$Chr),
											Loc=as.character(site.table$Loc), updated.island=island.clusters)
			if(output.format == 'xls'){
				xlsfile = file.path(data.folder, paste(project.name,"_DeNovo_Site_to_Island_Mapping-Avg_by_Site.xlsx",sep=""))
				WriteXLS("mapping.table", ExcelFileName = xlsfile)
			} else if(output.format == 'csv'){
				txtfile = file.path(data.folder, paste(project.name,"_DeNovo_Site_to_Island_Mapping-Avg_by_Site.csv",sep=""))
				write.table(mapping.table, file=txtfile, quote=F, row.names=F, sep=",")
			} else if(output.format == 'txt'){
				txtfile = file.path(data.folder, paste(project.name,"_DeNovo_Site_to_Island_Mapping-Avg_by_Site.txt",sep=""))
				write.table(mapping.table, file=txtfile, quote=F, row.names=F, sep="\t")
			} else {
				warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
			}
			rm(mapping.table)
			
			site.table[,5] = island.clusters
			site.table = site.table[!is.na(site.table[[5]]), ]
		} else{
			print("Can only define clusters with 2-group (or continuous) comparison.  Sorry.")
		}
	}#revise island annotations
	
	#print(site.table[[5]])
	
	cpg.islands = levels(as.factor(as.character(site.table[[5]])))
	print(paste("Analyzing",length(cpg.islands),"CpG Islands...",sep=" "))
	genes = array(dim=length(cpg.islands))
	num.methyl = array(dim=length(cpg.islands))
	num.unmethyl = array(dim=length(cpg.islands))
	delta.beta = array(dim=length(cpg.islands))
	cpg.island.pvalue = array(dim=length(cpg.islands))
	cpg.status = array(dim=length(cpg.islands))

	if(num.groups == 1) {
	avg.beta = site.table[[ncol(site.table)]]
	bg.up = length(avg.beta[avg.beta > methyl.cutoff])
	bg.down = length(avg.beta[avg.beta < unmethyl.cutoff])
	#print(bg.up)
	#print(bg.down)
	for (i in 1:length(cpg.islands))
		{
			island = cpg.islands[i]
			data.mat = site.table[site.table[5] == island,]
			genes[i] = as.character(data.mat[1,4])
			#print(data.mat)
			#print(genes[i])
			
			island.beta = data.mat[[ncol(data.mat)]]
			num.methyl[i] = length(island.beta[island.beta > methyl.cutoff])
			num.unmethyl[i] = length(island.beta[island.beta < unmethyl.cutoff])
			#print(num.methyl[i])
			#print(num.unmethyl[i])
			if(num.methyl[i] > num.unmethyl[i])
				{
					cpg.status[i] = "Methylated"
				}#end if(mean(trt.beta) > mean(ref.beta))
			if(num.unmethyl[i] > num.methyl[i])
				{
					cpg.status[i] = "Unmethylated"
				}#end if(mean(trt.beta) > mean(ref.beta))
			#print(cpg.status[i])
			#stop()
			#print(summary(test))
			if((num.methyl[i] + num.unmethyl[i]) > 2)
				{
					methyl.mat = matrix(c(num.methyl[i] , num.unmethyl[i], bg.up, bg.down), nrow=2, ncol=2,
							 dimnames = list(c("Methyl","Unmethyl"),c("Sample","Background")))
					cpg.island.pvalue[i] = fisher.test(methyl.mat)$p.value
				}
			else
				{
					cpg.island.pvalue[i] =1
				}
			#print(length(summary(test)[[1]][['Pr(>F)']][1]))

			#print(cpg.island.pvalue)
			#stop()
		}#end for (i in 1:length(cpg.islands))
	cpg.island.fdr = p.adjust(cpg.island.pvalue, method="fdr")
	island.table = data.frame( cpg.island = cpg.islands,
								genes = genes,
								num.methyl = num.methyl,
								num.unmethyl = num.unmethyl,
								island.pvalue = cpg.island.pvalue,
								island.fdr = cpg.island.fdr,
								methylation.status = cpg.status)
	} else if (num.groups == 2) {
	for (i in 1:length(cpg.islands))
		{
			island = cpg.islands[i]
			#print(island)
			data.mat = site.table[site.table[5] == island,]
			genes[i] = as.character(data.mat[1,4])
			#print(data.mat)
			#print(genes[i])
			
			trt.beta = data.mat[[6]]
			ref.beta = data.mat[[7]]
			mean.trt.beta = mean(trt.beta)
			mean.ref.beta = mean(ref.beta)
			num.methyl[i] = length(trt.beta[mean.trt.beta > mean.ref.beta])
			num.unmethyl[i] = length(trt.beta[mean.trt.beta < mean.ref.beta])
			delta.beta[i] = mean.trt.beta - mean.ref.beta
			if(mean.trt.beta > mean.ref.beta)
				{
					cpg.status[i] = "Increased Methylation"
				}#end if(mean(trt.beta) > mean(ref.beta))
			if(mean.trt.beta < mean.ref.beta)
				{
					cpg.status[i] = "Decreased Methylation"
				}#end if(mean(trt.beta) > mean(ref.beta))
			#print(num.methyl[i])
			#print(num.unmethyl[i])
			comb.beta = c(trt.beta, ref.beta)
			comb.group = c(rep(1, length=nrow(data.mat)), rep(2, length=nrow(data.mat)))
			comb.site = c(1:nrow(data.mat),1:nrow(data.mat))
			#print(comb.group)
			#print(comb.site)
			test = aov(comb.beta ~ comb.group + comb.site ) 
			#print(summary(test))
			#stop()
			if(length(comb.site) > 2)
				{
					result = summary(test)[[1]][['Pr(>F)']][1]
					if(is.null(result))
						{
							cpg.island.pvalue[i] =1
						}
					else
						{
							cpg.island.pvalue[i] = result 
						}
				}
			else
				{
					cpg.island.pvalue[i] =1
				}
			#print(length(summary(test)[[1]][['Pr(>F)']][1]))

			#print(cpg.island.pvalue)
			#stop()
		}#end for (i in 1:length(cpg.islands))
	cpg.island.fdr = p.adjust(cpg.island.pvalue, method="fdr")
	island.table = data.frame( cpg.island = cpg.islands,
								genes = genes,
								num.methyl = num.methyl,
								num.unmethyl = num.unmethyl,
								island.delta.beta = delta.beta,
								island.pvalue = cpg.island.pvalue,
								island.fdr = cpg.island.fdr,
								methylation.status = cpg.status)
	} else {
	for (i in 1:length(cpg.islands))
		{
			island = cpg.islands[i]
			data.mat = site.table[site.table[5] == island,]
			genes[i] = as.character(data.mat[1,4])
			#print(data.mat)
			#print(genes[i])

			#print(ncol(data.mat))
			data.beta = data.mat[,6:(ncol(data.mat)-2)]
			#print(dim(data.beta))
			comb.beta = as.vector(as.matrix(data.beta))
			#print(length(comb.beta))
			comb.group = array()
			for (j in 1:ncol(data.beta))
				{
					comb.group = c(comb.group,rep(j, length=nrow(data.mat)))
				}
			comb.group = comb.group[-1]
			#print(length(comb.group))
			#print(comb.group)
			comb.site = rep(1:nrow(data.mat),ncol(data.beta))
			#print(length(comb.site))
			num.methyl[i] = nrow(data.mat)
			test = aov(comb.beta ~ comb.group + comb.site ) 
			#print(summary(test))
			if(length(comb.site) > 2)
				{
					result = summary(test)[[1]][['Pr(>F)']][1]
					if(is.null(result))
						{
							cpg.island.pvalue[i] =1
						}
					else
						{
							cpg.island.pvalue[i] = result 
						}
				}
			else
				{
					cpg.island.pvalue[i] =1
				}
			#print(length(summary(test)[[1]][['Pr(>F)']][1]))

			#print(cpg.island.pvalue)
			#stop()
		}#end for (i in 1:length(cpg.islands))
	cpg.island.fdr = p.adjust(cpg.island.pvalue, method="fdr")
	island.table = data.frame( cpg.island = cpg.islands,
								genes = genes,
								num.sites = num.methyl,
								island.pvalue = cpg.island.pvalue,
								island.fdr = cpg.island.fdr)
	}

	if(output.format == 'xls'){
		xlsfile = file.path(data.folder, paste(project.name,"_CpG_island_stats-Avg_by_Site.xlsx",sep=""))
		WriteXLS("island.table", ExcelFileName = xlsfile)
	} else if(output.format == 'csv'){
		txtfile = file.path(data.folder, paste(project.name,"_CpG_island_stats-Avg_by_Site.csv",sep=""))
		write.table(island.table, file=txtfile, quote=F, row.names=F, sep=",")
	} else if(output.format == 'txt'){
		txtfile = file.path(data.folder, paste(project.name,"_CpG_island_stats-Avg_by_Site.txt",sep=""))
		write.table(island.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	island.table = island.table[island.table$island.pvalue <= pvalue.cutoff,]
	island.table = island.table[island.table$island.fdr <= fdr.cutoff,]
	if(num.groups ==1){
		island.table = island.table[(island.table$num.methyl + island.table$num.unmethyl)>= num.sites,]
	}else if(num.groups == 2) {
		island.table = island.table[(island.table$num.methyl + island.table$num.unmethyl)>= num.sites,]
		island.table = island.table[abs(island.table$island.delta.beta) >= delta.beta.cutoff,]
	} else{
		island.table = island.table[island.table$num.sites >= num.sites,]
	}
	
	if(output.format == 'xls'){
		xlsfile = file.path(island.folder, paste(project.name,"_CpG_island_filtered-Avg_by_Site.xlsx",sep=""))
		WriteXLS("island.table", ExcelFileName = xlsfile)
	} else if(output.format == 'csv'){
		txtfile = file.path(island.folder, paste(project.name,"_CpG_island_filtered-Avg_by_Site.csv",sep=""))
		write.table(island.table, file=txtfile, quote=F, row.names=F, sep=",")
	} else if(output.format == 'txt'){
		txtfile = file.path(island.folder, paste(project.name,"_CpG_island_filtered-Avg_by_Site.txt",sep=""))
		write.table(island.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	return(island.table)
}#end def RNA.deg
