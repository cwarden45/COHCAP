cor.test.est = function(x)
	{
		index = 1:length(x)
		methyl = x[(index %% 2) == 1]
		expr = x[(index %% 2) == 0]
		
		paired.table = data.frame(methyl,expr)
		paired.table = na.omit(paired.table)
		if(nrow(paired.table) < 3){
			return(0)
		}else{
			result = cor.test(methyl, expr, na.action=na.rm)
			if(is.na(result$estimate)){
				return(0)
			}else{
				return(result$estimate)
			}
		}#end else
	}#end def cor.test.est
	
cor.test.pvalue = function(x)
	{
		index = 1:length(x)
		methyl = x[(index %% 2) == 1]
		expr = x[(index %% 2) == 0]
		
		paired.table = data.frame(methyl,expr)
		paired.table = na.omit(paired.table)
		if(nrow(paired.table) < 3){
			return(0)
		}else{
			result = cor.test(methyl, expr, na.action=na.rm)
			if(is.na(result$p.value)){
				return(1)
			}else{
				return(result$p.value)
			}
		}#end else
	}#end def cor.test.est

`COHCAP.integrate.avg.by.island` = function (island.list, project.name, project.folder, expr.file, sample.file, cor.pvalue.cutoff=0.05, cor.fdr.cutoff = 0.05, cor.cutoff = -0.2, plot.scatter=TRUE, output.format = "txt", ref="none")
{
	fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",colors())

	beta.table = island.list$beta.table
	print(dim(beta.table))
	filtered.island.stats = island.list$filtered.island.stats
	print(dim(filtered.island.stats))
	
	integrate.folder=file.path(project.folder,"Integrate")
	dir.create(integrate.folder, showWarnings=FALSE)
	
	data.folder=file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	tempID = sample(100000:999999,1)
	temp.paired = file.path(data.folder,paste("temp_paired_",tempID,".txt",sep=""))
	temp.methyl = file.path(data.folder,paste("temp_methyl_",tempID,".txt",sep=""))
	write.table(beta.table, temp.methyl, quote=F, row.names=F, sep="\t")
	
	Perl.Path = file.path(path.package("COHCAP"), "Perl")
	perl.script = file.path(Perl.Path , "integrate_island.pl")
	cmd = paste("perl \"",perl.script,"\" \"", expr.file,"\" \"", temp.methyl,"\" \"", temp.paired,"\"", sep="")
	res = system(cmd, intern=TRUE, wait=TRUE)
	message(res)

	input.table = read.table(temp.paired, header=T, sep = "\t", stringsAsFactors=TRUE)
	#print(head(input.table))
	islands = as.character(input.table[[1]])
	print(length(islands))
	genes = input.table[[2]]
	intensity.table = input.table[,3:ncol(input.table)]
	cor.coef.values = apply(intensity.table, 1, cor.test.est)
	cor.pvalues = apply(intensity.table, 1, cor.test.pvalue)
	cor.fdr.values = p.adjust(cor.pvalues, method="fdr")
	print(cor.pvalues)

	#may produce NA column with other than n=2 comparisons
	#I'll update code when I encounter a problem
	island.direction = rep("NA", times = length(islands))
	delta.beta = filtered.island.stats[[5]]
	print(length(delta.beta))
	delta.beta = delta.beta[match(islands, as.character(filtered.island.stats$island), nomatch=0)]
	print(length(delta.beta))
	island.direction[delta.beta > 0] = "Increased Methylation"
	island.direction[delta.beta < 0] = "Decreased Methylation"
	
	cor.table= data.frame(island=islands, island.direction = island.direction, gene=genes, cor=cor.coef.values, p.value=cor.pvalues, fdr=cor.fdr.values)
	if(output.format == 'xls'){
		xlsfile = file.path(data.folder, paste(project.name,"_integration_cor_stats.xlsx",sep=""))
		WriteXLS("cor.table", ExcelFileName = xlsfile)
	} else if(output.format == 'csv'){
		txtfile = file.path(data.folder, paste(project.name,"_integration_cor_stats.csv",sep=""))
		write.table(cor.table, file=txtfile, quote=F, row.names=F, sep=",")
	} else if(output.format == 'txt'){
		txtfile = file.path(data.folder, paste(project.name,"_integration_cor_stats.txt",sep=""))
		write.table(cor.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	sig.islands = cor.table$island[(cor.table$cor < cor.cutoff) & (cor.table$p.value < cor.pvalue.cutoff) & (cor.table$fdr < cor.fdr.cutoff)]
	print(paste(length(sig.islands)," significant correlations",sep=""))
	
	if(length(sig.islands) > 0){
		sig.table = cor.table[(cor.table$cor < cor.cutoff) & (cor.table$p.value < cor.pvalue.cutoff) & (cor.table$fdr < cor.fdr.cutoff),]
		if(output.format == 'xls'){
			xlsfile = file.path(integrate.folder, paste(project.name,"_integration_cor_stats.xlsx",sep=""))
			WriteXLS("sig.table", ExcelFileName = xlsfile)
		} else if(output.format == 'csv'){
			txtfile = file.path(integrate.folder, paste(project.name,"_integration_cor_stats.csv",sep=""))
			write.table(sig.table, file=txtfile, quote=F, row.names=F, sep=",")
		} else if(output.format == 'txt'){
			txtfile = file.path(integrate.folder, paste(project.name,"_integration_cor_stats.txt",sep=""))
			write.table(sig.table, file=txtfile, quote=F, row.names=F, sep="\t")
		} else {
			warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
		}
	}#end if(length(sig.islands) > 0)
	
	if(length(sig.islands) > 0)
		{	
			if(plot.scatter)
				{
					scatter.folder=file.path(integrate.folder,paste(project.name,"_Cor",sep=""))
					dir.create(scatter.folder, showWarnings=FALSE)
					
					sig.table = input.table[(cor.table$cor < cor.cutoff) & (cor.table$p.value < cor.pvalue.cutoff) & (cor.table$fdr < cor.fdr.cutoff),]
					if(output.format == 'xls'){
						xlsfile = file.path(data.folder, paste(project.name,"_cor_filter.xlsx",sep=""))
						WriteXLS("sig.table", ExcelFileName = xlsfile)
					} else if(output.format == 'txt'){
						txtfile = file.path(data.folder, paste(project.name,"_cor_filter.txt",sep=""))
						write.table(sig.table, file=txtfile, quote=F, row.names=F, sep="\t")
					} else {
						warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
					}
					rm(input.table)
				
					print("Plotting Correlated Genes....")
					print(dim(sig.table))
				
					sample.names = names(sig.table)[3:ncol(sig.table)]
					index = 1:length(sample.names)
					sample.names = sample.names[index %% 2 == 1]
					sample.names = gsub(".Methyl","",sample.names)
					#print(sample.names)
				
					sample.table = read.table(sample.file, header=F, sep = "\t", stringsAsFactors=TRUE)
					samples = as.character(sample.table[[1]])
					samples[grep("^\\d",samples, perl = TRUE)] = paste("X",samples[grep("^\\d",samples, perl = TRUE)],sep="")
					samples = gsub("-",".",samples)
					#print(samples)
					sample.group = sample.table[[2]]
					#print(sample.group)
					sample.group = sample.group[match(sample.names,samples,nomatch=0)]
					#print(sample.group)
					labelColors =  as.character(sample.group)
					if(ref == "continuous"){
						continous.var = sample.group
						continuous.color.breaks = 10	
						plot.var = as.numeric(continous.var)
						plot.var.min = min(plot.var, na.rm=T)
						plot.var.max = max(plot.var, na.rm=T)
						
						plot.var.range = plot.var.max - plot.var.min
						plot.var.interval = plot.var.range / continuous.color.breaks
						
						color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
						plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
						for (j in 1:continuous.color.breaks){
							labelColors[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
						}#end for (j in 1:continuous.color.breaks)
					}else{
						groups = levels(sample.group)
						print(groups)
						color.palette = fixed.color.palatte[1:length(groups)]
						for (j in 1:length(groups)){
							labelColors[sample.group == as.character(groups[j])] = color.palette[j]
						}#end for (j in 1:length(groups))
					}
				
					for (i in 1:length(sig.islands))
						{
							island = as.character(sig.table[i,1])
							gene = as.character(sig.table[i,2])
						
							island = gsub(":","_",island)
							island = gsub("-","_",island)
							plot.file = file.path(scatter.folder, paste(gene,"_",island,".pdf", sep=""))
							data = sig.table[i,3:ncol(sig.table)]
							index = 1:length(data)
							methyl = t(data[(index %% 2) == 1])
							expr = t(data[(index %% 2) == 0])
							pdf(file=plot.file)
							plot(expr, methyl, pch=20, col=labelColors,
								main=paste(gene,
									" (p=",
									round(cor.test(expr,methyl, na.action=na.rm)$p.value, digits=2),
									" ,r=",
									round(cor.test(expr,methyl, na.action=na.rm)$estimate, digits=2),
									")",sep=""))
							if(ref == "continuous"){
								legend("topright",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
										col=rev(color.range), pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
							}else{
								legend("topright",legend=groups,col=color.palette,  pch=15)
							}
							dev.off()
						}#end for (1 in 1:length(sig.islands))
				}	
		}#end if(length(sig.islands) > 0)

	unlink(temp.paired)
	unlink(temp.methyl)
}#end def RNA.deg
