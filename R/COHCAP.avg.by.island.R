custom.mean <- function(arr)
{
	return(mean(as.numeric(arr), na.rm = T))
}#end def ttest2

ttest2 <- function(arr, grp1, grp2)
{
	group1 <- as.numeric(arr[grp1])
	group2 <- as.numeric(arr[grp2])
	#print(grp1)
	#print(grp2)
	#stop()
	#require at least 2 replicates
	if((length(group1[!is.na(group1)]) >=3) & (length(group2[!is.na(group2)]) >=3))
	{
		result <- 1
		tryCatch(result <- t.test(group1, group2)$p.value, error = function(e) {})
		return(result)
	}
	else
	{
		return(1)
	}
}#end def ttest2

annova.pvalue <- function(arr, grp.levels)
{
	#print(arr)
	#print(grp.levels)
	
	grp.no.na <- as.factor(as.character(grp.levels[!is.na(arr)]))
	if(length(levels(grp.no.na)) >= 2)
		{
			test <- aov(arr ~ grp.levels) 
			result <- summary(test)[[1]][['Pr(>F)']][1]
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

lm.pvalue = function(arr, var1)
{
	fit = lm(as.numeric(arr)~var1)
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	return(pvalue)
}#end def lm.pvalue


lm.pvalue2 = function(arr, var1, var2)
{
	fit = lm(as.numeric(arr)~var1 + as.numeric(as.factor(var2)))
	result = summary(fit)
	pvalue = result$coefficients[2,4]
	return(pvalue)
}#end def lm.pvalue

annova.2way.pvalue <- function(arr, grp.levels, pairing.levels)
{
	#print(arr)
	#print(grp.levels)
	#print(pairing.levels)
	
	grp.no.na <- as.factor(as.character(grp.levels[!is.na(arr)]))
	
	rep.flag <- 1
	if((length(levels(grp.no.na)) >= 2) && (rep.flag == 1))
		{
			test <- aov(arr ~ grp.levels + pairing.levels) 
			result <- summary(test)[[1]][['Pr(>F)']][1]
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

`COHCAP.avg.by.island` <-function (sample.file, site.table, beta.table, project.name, project.folder, methyl.cutoff=0.7, unmethyl.cutoff = 0.3, delta.beta.cutoff = 0.2, pvalue.cutoff=0.05, fdr.cutoff=0.05, num.groups=2, num.sites=4, plot.box=TRUE, paired=FALSE, ref="none", output.format = "xls", gene.centric=TRUE)
{
	island.folder<-file.path(project.folder,"CpG_Island")
	dir.create(island.folder, showWarnings=FALSE)
	
	data.folder<-file.path(project.folder,"Raw_Data")
	dir.create(data.folder, showWarnings=FALSE)
	
	print("Reading Sample Description File....")
	sample.table <- read.table(sample.file, header=F, sep = "\t")
	samples <- as.character(sample.table[[1]])
	for (i in 1:length(samples))
		{
			if(length(grep("^[0-9]",samples[i])) > 0)
				{
					samples[i] <- paste("X",samples[i],sep="")
				}#end if(length(grep("^[0-9]",samples[i])) > 0)
		}#end def for (i in 1:length(samples))
	sample.group <- sample.table[[2]]
	sample.group <- as.factor(gsub(" ",".",sample.group))
	pairing.group <- NA
	if(paired)
		{
			pairing.group <- sample.table[[3]]
			pairing.group <- as.factor(gsub(" ",".",pairing.group))	
		}
	ref <- gsub(" ",".",ref)

	sample.names <- names(beta.table)[6:ncol(beta.table)]
	beta.values <- beta.table[,6:ncol(beta.table)]
	annotation.names <- names(beta.table)[1:5]
	annotations <- beta.table[,1:5]
	print(dim(beta.values))

	if(length(samples) != length(sample.names[match(samples, sample.names, nomatch=0)]))
		{
			warning("Some samples in sample description file are not present in the beta file!")
			warning(paste(length(samples),"items in sample description file",sep=" "))
			warning(paste(length(sample.names),"items in gene beta file",sep=" "))
			warning(paste(length(sample.names[match(samples, sample.names, nomatch=0)]),"matching items in gene beta file",sep=" "))
			#warning(sample.names[match(samples, sample.names, nomatch=0)])
			stop()
		}

	if(length(samples)>1)
		{
			beta.values <- beta.values[,match(samples, sample.names, nomatch=0)]
			colnames(beta.values)=samples
		}
		
	print(dim(beta.values))
	#print(samples)
	#print(sample.names)
	#print(colnames(beta.values))
	#print(beta.values[1,])
	rm(beta.table)
	
	groups <- levels(sample.group)
	print(paste("Group:",groups, sep=" "))

	if((length(groups) != num.groups) & (ref != "continuous"))
		{
			print(paste("Group:",groups, sep=" "))

			warning(paste("There are ",length(groups)," but user specified algorithm for ",num.groups," groups.",sep=""))
			warning("Please restart algorithm and specify the correct number of groups")
			stop()
		}
	
	if (ref == "continuous"){
		print("Continous Variable :")
		print(sample.table[[2]])
	}
	
	print("Checking CpG Site Stats Table")
	print(dim(site.table))
	site.table <- site.table[!is.na(site.table[[5]]), ]
	print(dim(site.table))
	
	cpg.islands <- levels(as.factor(as.character(site.table[[5]])))
	print(length(cpg.islands))

	genes <- array(dim=length(cpg.islands))
	island.values<- matrix(nrow=length(cpg.islands), ncol=ncol(beta.values))
	colnames(island.values) <- samples

	#average CpG sites per CpG Island
	print("Average CpG Sites per CpG Island")
	for (i in 1:length(cpg.islands))
		{
			island <- cpg.islands[i]
			#print(island)
			cpg.sites <- site.table[site.table[[5]] == island,1]
			#print(cpg.sites)
			data.mat <- beta.values[match(cpg.sites, annotations[[1]],nomatch=0),]
			#print(data.mat)
			if(nrow(data.mat) >= num.sites)
				{
					info.mat <- site.table[site.table[5] == island,]
					genes[i] <- NA
					if(length(levels(as.factor(as.character(info.mat[[4]])))) == 1)
						{
							genes[i] <- as.character(info.mat[1,4])
						}
					island.values[i,]=apply(data.mat,2,mean, na.rm=T)
				}#end if(nrow(data.mat) >= num.sites)
		}#end for (i in 1:length(cpg.islands))
	island.avg.table <- data.frame(island=cpg.islands, gene=genes, island.values)
	if(gene.centric)
		{
			island.avg.table <- island.avg.table[!is.na(island.avg.table[[2]]), ]
		}
	annotations <- island.avg.table[,1:2]
	annotation.names <- c("island","gene")
	beta.values <- island.avg.table[,3:ncol(island.avg.table)]
	#print(beta.values)
	
	rm(site.table)
	rm(cpg.islands)
	rm(genes)
	rm(island.values)

	#CpG Island statistical analysis
	island.table <- annotations 
	rm(annotations)
	
	if(ref == "continuous"){
			continous.var = sample.table[[2]]
			if (paired == TRUE){
				lm.pvalue = apply(beta.values, 1, lm.pvalue2, continous.var, pairing.group)
			} else{
				lm.pvalue = apply(beta.values, 1, lm.pvalue, continous.var)
			}
			lm.fdr = p.adjust(lm.pvalue, "fdr")
			beta.cor = apply(beta.values, 1, custom.cor, var1=continous.var)
			#upper.beta is max if positive correlation, min if negative correlation
			#lower.beta is min if positive correlation, max if negative correlation
			#in both cases, delta.beta is upper.beta - lower.beta
			beta.min= apply(beta.values, 1, min, na.rm=TRUE)
			beta.max= apply(beta.values, 1, max, na.rm=TRUE)
			
			upper.beta = beta.max
			upper.beta[beta.cor < 0] = beta.min[beta.cor < 0]
			lower.beta = beta.min
			lower.beta[beta.cor < 0] = beta.max[beta.cor < 0]
			
			rm(beta.min)
			rm(beta.max)
			
			delta.beta = upper.beta - lower.beta

			#make format compatible with 2-group code
			island.table = data.frame(island.table, upper.beta = upper.beta, lower.beta = lower.beta,
									delta.beta = delta.beta, island.pvalue = lm.pvalue, island.fdr = lm.fdr,
									cor.coef = beta.cor)
	} else if(length(groups) == 1) {
	col.names <- c(annotation.names)
	print("Averaging Beta Values for All Samples")

	if(length(sample.names) > 1)
		{
			avg.beta <- apply(beta.values, 1, mean, na.rm = T)
		}
	else
		{
			avg.beta<- beta.values
		}
	island.table <- data.frame(island.table, avg.beta=avg.beta)
	} else if(length(groups) == 2) {
	print("Differential Methylation Stats for 2 Groups with Reference")
	trt <- groups[groups != ref]
	#print(ref)
	#print(trt)
	#print(samples)
	#print(sample.group)
	trt.samples <- samples[which(sample.group == trt)]
	#print(trt.samples)
	ref.samples <- samples[which(sample.group == ref)]
	#print(ref.samples)

	all.indices <- 1:length(sample.names)
	trt.indices <- all.indices[which(sample.group == trt)]
	#print(trt.indices)
	ref.indices <- all.indices[which(sample.group == ref)]
	#print(ref.indices)

	trt.beta.values <- beta.values[, trt.indices]
	ref.beta.values <- beta.values[, ref.indices]

	if(length(trt.indices) > 1)
		{
			trt.avg.beta <- apply(trt.beta.values, 1, mean, na.rm = T)
		}
	else
		{
			trt.avg.beta <- trt.beta.values
		}
		
	if(length(ref.indices) > 1)
		{
			ref.avg.beta <- apply(ref.beta.values, 1, mean, na.rm = T)
		}
	else
		{
			ref.avg.beta<- ref.beta.values
		}


	delta.beta <- trt.avg.beta - ref.avg.beta

	if(paired)
		{
			print("Factor in Paired Samples")
			beta.pvalue <- apply(beta.values, 1, annova.2way.pvalue, grp.levels=sample.group, pairing.levels=pairing.group)
		}#end if(pair.flag == 1)
	else
		{
			beta.pvalue <- apply(beta.values, 1, ttest2, grp1=trt.indices, grp2=ref.indices)
		}#end else
	beta.fdr <- p.adjust(beta.pvalue, method="fdr")

	#print(dim(island.table))
	#print(length(trt.avg.beta))
	#print(length(ref.avg.beta))
	#print(length(delta.beta))
	#print(length(beta.pvalue))
	#print(length(beta.fdr))

	island.table <- data.frame(island.table, trt.beta = trt.avg.beta, ref.beta = ref.avg.beta, delta.beta = delta.beta, pvalue = beta.pvalue, fdr = beta.fdr)
	colnames(island.table) <- c(annotation.names,
								paste(trt,"avg.beta",sep="."),
								paste(ref,"avg.beta",sep="."),
								paste(trt,"vs",ref,"delta.beta",sep="."),
								"island.pvalue",
								"island.fdr")
	} else if(length(groups) > 2) {
	print("Calculating Average Beta and ANOVA p-values")
	col.names <- c(annotation.names)

	for (i in 1:length(groups))
	{
		temp.grp <- groups[i]
		all.indices <- 1:length(sample.names)
		grp.indices <- all.indices[which(sample.group == temp.grp)]
		grp.beta <- beta.values[,which(sample.group == temp.grp)]
		if(length(grp.indices) > 1)
			{
				avg.beta <- apply(grp.beta, 1, mean, na.rm = T)
			}
		else
			{
				avg.beta<- grp.beta
			}
		col.names <- c(col.names, paste(temp.grp,"avg.beta",sep="."))
		island.table <- data.frame(island.table, avg.beta=avg.beta)
		colnames(island.table) <- col.names
	}#end for (i in 1:length(groups)

	if(paired)
		{
			print("Factor in Paired Samples")
			pvalue <- apply(beta.values, 1, annova.2way.pvalue, grp.levels=sample.group, pairing.levels=pairing.group)
		}#end if(pair.flag == 1)
	else
		{
			pvalue <- apply(beta.values, 1, annova.pvalue, grp.levels = sample.group)
		}#end else
	anova.fdr <- p.adjust(pvalue, method="fdr")
	col.names <- c(col.names, "island.pvalue", "island.fdr")
	island.table <- data.frame(island.table, anova.pvalue = pvalue, anova.fdr = anova.fdr)
	colnames(island.table) <- col.names
	} else {
	warning("Invalid groups specified in sample description file!")
	}
	

	if(output.format == 'xls'){
		xlsfile <- file.path(data.folder, paste(project.name,"_CpG_island_stats-Avg_by_Island.xlsx",sep=""))
		WriteXLS("island.table", ExcelFileName = xlsfile)
	} else if(output.format == 'txt'){
		txtfile <- file.path(data.folder, paste(project.name,"_CpG_island_stats-Avg_by_Island.txt",sep=""))
		write.table(island.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	#filter CpG islands
	print(dim(island.table))
	filter.table = island.table
	if((length(groups) == 2)|(ref == "continuous")){
		temp.trt.beta <- island.table[[3]]
		temp.ref.beta <- island.table[[4]]
		temp.delta.beta <- abs(island.table[[5]])
		temp.pvalue <- island.table[[6]]
		temp.fdr <- island.table[[7]]
		
		if(unmethyl.cutoff > methyl.cutoff)
			{
				print("unmethyl.cutoff > methyl.cutoff ...")
				print("so, delta-beta decides which group should be subject to which cutoff")
				temp.delta.beta <- island.table[[5]]
				temp.methyl.up <- filter.table[(temp.trt.beta >= methyl.cutoff) & (temp.ref.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta >= delta.beta.cutoff),]
				temp.methyl.down <- filter.table[(temp.ref.beta >= methyl.cutoff) & (temp.trt.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta <= -delta.beta.cutoff),]
				filter.table <- rbind(temp.methyl.up, temp.methyl.down)	}
		else
			{
				temp.methyl.up <- filter.table[(temp.trt.beta >= methyl.cutoff) & (temp.ref.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta >= delta.beta.cutoff),]
				temp.methyl.down <- filter.table[(temp.ref.beta >= methyl.cutoff) & (temp.trt.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta >= delta.beta.cutoff),]
				filter.table <- rbind(temp.methyl.up, temp.methyl.down)
			}
	} else if(length(groups) == 1) {
		temp.avg.beta <- island.table$avg.beta
		filter.table <- filter.table[(temp.avg.beta >= methyl.cutoff) | (temp.avg.beta <=unmethyl.cutoff),]
	} else if(length(groups) > 2) {
		temp.pvalue <- island.table[[ncol(island.table) -1]]
		temp.fdr <- island.table[[ncol(island.table)]]
		filter.table <- filter.table[(temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff),]
	} else {
	warning("Invalid groups specified in sample description file!")
	}
	print(dim(filter.table))
	
	sig.islands <- filter.table[[1]]
	print(paste("There are ",length(sig.islands)," differentially methylated islands",sep=""))
	
	if(output.format == 'xls'){
		xlsfile <- file.path(island.folder, paste(project.name,"_CpG_island_filtered-Avg_by_Island.xlsx",sep=""))
		WriteXLS("filter.table", ExcelFileName = xlsfile)
	} else if(output.format == 'txt'){
		txtfile <- file.path(island.folder, paste(project.name,"_CpG_island_filtered-Avg_by_Island.txt",sep=""))
		write.table(filter.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	print(dim(island.avg.table))
	island.avg.table <- island.avg.table[match(sig.islands,island.avg.table[[1]],nomatch=0),]
	print(dim(island.avg.table))


	if(output.format == 'xls'){
		xlsfile <- file.path(data.folder, paste(project.name,"_CpG_island_filtered_beta_values-Avg_by_Island.xlsx",sep=""))
		WriteXLS("island.avg.table", ExcelFileName = xlsfile)
	} else if(output.format == 'txt'){
		txtfile <- file.path(data.folder, paste(project.name,"_CpG_island_filtered_beta_values-Avg_by_Island.txt",sep=""))
		write.table(island.avg.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
if((plot.box) & (nrow(island.avg.table) > 0))
	{
		if (ref == "continuous"){
			print("Plotting Significant Islands for Continuous Variable..")
			continous.var = sample.table[[2]]
			
			plot.folder<-file.path(island.folder,paste(project.name,"_Line_Plots",sep=""))
			dir.create(plot.folder, showWarnings=FALSE)
			
			for (i in 1:length(sig.islands))
				{
					island = as.character(island.avg.table[i,1])
					gene = as.character(island.avg.table[i,2])
					
					gene = gsub(";","_",gene)
					
					island = gsub(":","_",island)
					island = gsub("-","_",island)
					plot.file = file.path(plot.folder, paste(gene,"_",island,".pdf", sep=""))
					methyl.level = as.numeric(island.avg.table[i,3:ncol(island.avg.table)])
					
					pdf(file=plot.file)
					plot(continous.var, methyl.level, pch=16, col="black",
						main=paste(gene,sep=""), ylab="Methylation (Beta)", xlab="")
					abline(lm(methyl.level~continous.var),lwd=2, col="red")
					dev.off()
				}#end for (1 in 1:length(sig.islands))
		}else{
			print("Plotting Significant Islands Box-Plots..")
			plot.folder<-file.path(island.folder,paste(project.name,"_Box_Plots",sep=""))
			dir.create(plot.folder, showWarnings=FALSE)
			
			labelColors <- as.character(sample.group)
			
			colors <- c("red","blue","green","orange","purple","cyan","pink","maroon","yellow","grey","black",colors())
			colors <- colors[1:length(groups)]
			
			#print(samples)
			#print(sample.names)
			
			for (i in 1:length(sig.islands))
				{
					island = as.character(island.avg.table[i,1])
					gene = as.character(island.avg.table[i,2])
					
					gene = gsub(";","_",gene)
					
					island = gsub(":","_",island)
					island = gsub("-","_",island)
					plot.file = file.path(plot.folder, paste(gene,"_",island,".pdf", sep=""))
					data = t(island.avg.table[i,3:ncol(island.avg.table)])
					
					pdf(file=plot.file)
					plot(sample.group, data, pch=20, col=colors,
						main=paste(gene,sep=""), ylab="Methylation (Beta)")
					dev.off()
				}#end for (1 in 1:length(sig.islands))
		}#end else
	}#end if((plot.box) & (nrow(island.avg.table) > 0))
	
	integrate.tables = list(beta.table=island.avg.table, filtered.island.stats=filter.table)
	return(integrate.tables)
}#end def RNA.deg
