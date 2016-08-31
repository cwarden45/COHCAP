custom.mean = function(arr)
{
	return(mean(as.numeric(arr), na.rm = T))
}#end def custom.mean

custom.cor = function(arr, var1)
{
	return(cor(as.numeric(arr), var1, use="complete.obs"))
}#end def custom.cor

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

`COHCAP.site` <-function (sample.file, beta.table, project.name, project.folder, methyl.cutoff=0.7, unmethyl.cutoff = 0.3, delta.beta.cutoff = 0.2, pvalue.cutoff=0.05, fdr.cutoff=0.05, ref="none", num.groups=2, create.wig = "avg", paired=FALSE, output.format='xls')
{
	site.folder<-file.path(project.folder,"CpG_Site")
	dir.create(site.folder, showWarnings=FALSE)
	
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

	#print(samples)
	#print(sample.names)

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

	groups = levels(sample.group)
	
	if((length(groups) != num.groups) & (ref != "continuous"))
		{
			print(paste("Group:",groups, sep=" "))

			warning(paste("There are ",length(groups)," but user specified algorithm for ",num.groups," groups.",sep=""))
			warning("Please restart algorithm and specify the correct number of groups")
			stop()
		}
	
	if (ref == "continuous"){
		print(paste("Continous Variable : ",paste(sample.table[[2]],collapse=","),sep=""))
	}

	stat.table <- annotations

	if((create.wig == "sample")|(create.wig == "avg.and.sample")){
			#create .wig files for each sample
			wig.folder<-file.path(site.folder,paste(project.name,"wig",sep="_"))
			dir.create(wig.folder, showWarnings=FALSE)

			temp.stat.file <- file.path(wig.folder, "temp.txt")
			Perl.Path <- file.path(path.package("COHCAP"), "Perl")
			perl.script <- file.path(Perl.Path , "create_wig_files.pl")
			
			sample.beta = beta.table
			sample.beta = sample.beta[order(sample.beta$Chr, sample.beta$Loc), ]
			colnames(sample.beta) = paste(colnames(beta.table),"beta",sep=".")
			write.table(sample.beta, temp.stat.file, quote=F, row.names=F, sep="\t")
			cmd <- paste("perl \"",perl.script,"\" \"", temp.stat.file,"\" \"", wig.folder,"\"", sep="")
			res <- system(cmd, intern=TRUE, wait=TRUE)
			message(res)

			shortcut.beta.file = file.path(wig.folder, "SiteID.beta.wig")
			unlink(shortcut.beta.file)
			
			shortcut.beta.file = file.path(wig.folder, "Chr.beta.wig")
			unlink(shortcut.beta.file)

			shortcut.beta.file = file.path(wig.folder, "Gene.beta.wig")
			unlink(shortcut.beta.file)
			
			shortcut.beta.file = file.path(wig.folder, "Loc.beta.wig")
			unlink(shortcut.beta.file)

			shortcut.beta.file = file.path(wig.folder, "Island.beta.wig")
			unlink(shortcut.beta.file)
		}#end if(create.wig)
	
	rm(beta.table)
	
	if(ref == "continuous"){
			continous.var = sample.table[[2]]
			lm.pvalue = apply(beta.values, 1, lm.pvalue, continous.var)
			lm.fdr = p.adjust(lm.pvalue, "fdr")
			beta.cor = apply(beta.values, 1, custom.cor, var1=continous.var)
			#upper.beta is max if positive correlation, min if negative correlation
			#lower.beta is min if positive correlation, max if negative correlation
			#in both cases, delta.beta is upper.beta - lower.beta
			beta.min= apply(beta.values, 1, min, na.rm=TRUE)
			beta.max= apply(beta.values, 1, min, na.rm=TRUE)
			
			upper.beta = beta.max
			upper.beta[beta.cor < 0] = beta.min[beta.cor < 0]
			lower.beta = beta.min
			lower.beta[beta.cor < 0] = beta.max[beta.cor < 0]
			
			rm(beta.min)
			rm(beta.max)
			
			delta.beta = upper.beta - lower.beta

			#make format compatible with 2-group code
			stat.table <- data.frame(stat.table, upper.beta = upper.beta, lower.beta = lower.beta,
									delta.beta = delta.beta, pvalue = lm.pvalue, fdr = lm.fdr,
									cor.coef = beta.cor)
	} else if(length(groups) == 1) {
	col.names <- c(annotation.names)
	print("Averaging Beta Values for All Samples")

	#print(beta.values[1:3,])
	if(length(samples) > 1)
		{
			avg.beta <- apply(beta.values, 1, custom.mean)
		}
	else
		{
			avg.beta<- beta.values
		}
	#print(avg.beta)

	stat.table <- data.frame(stat.table, avg.beta=avg.beta)
	} else if(length(groups) == 2) {
	print("Differential Methylation Stats for 2 Groups with Reference")
	trt <- groups[groups != ref]
	#print(trt)
	#print(ref)
	trt.samples <- samples[which(sample.group == trt)]
	#print(trt.samples)
	ref.samples <- samples[which(sample.group == ref)]
	#print(ref.samples)

	all.indices <- 1:length(samples)
	trt.indices <- all.indices[which(sample.group == trt)]
	ref.indices <- all.indices[which(sample.group == ref)]

	if (length(trt.indices) == 1){
		trt.avg.beta = beta.values[, trt.indices]
	} else {
		trt.beta.values <- beta.values[, trt.indices]	
		trt.avg.beta <- apply(trt.beta.values, 1, custom.mean)
	}
	
	if (length(ref.indices) == 1){
		ref.avg.beta = beta.values[, ref.indices]
	} else {
		ref.beta.values <- beta.values[, ref.indices]
		ref.avg.beta <- apply(ref.beta.values, 1, custom.mean)
	}

	delta.beta <- trt.avg.beta - ref.avg.beta

	if(paired)
		{
			print("Factor in Paired Samples")
			beta.pvalue <- unlist(apply(beta.values, 1, annova.2way.pvalue, grp.levels=sample.group, pairing.levels=pairing.group))
		}#end if(paired)
	else
		{
			beta.pvalue <- apply(beta.values, 1, ttest2, grp1=trt.indices, grp2=ref.indices)
		}#end else

	#print(beta.pvalue)
		
	beta.fdr <- p.adjust(beta.pvalue, method="fdr")

	#print(names(beta.values))
	print(dim(stat.table))
	#print(length(annotation.names))
	#print(length(trt.avg.beta))
	#print(length(ref.avg.beta))
	#print(length(delta.beta))
	#print(length(beta.pvalue))
	#print(length(beta.fdr))

	stat.table <- data.frame(stat.table, trt.beta = trt.avg.beta, ref.beta = ref.avg.beta, delta.beta = delta.beta, pvalue = beta.pvalue, fdr = beta.fdr)
	print(dim(stat.table))
	#print(names(stat.table))
	#test <- c(annotation.names,
	#							paste(trt,"avg.beta",sep="."),
	#							paste(ref,"avg.beta",sep="."),
	#							paste(trt,"vs",ref,"delta.beta",sep="."),
	#							paste(trt,"vs",ref,"pvalue",sep="."),
	#							paste(trt,"vs",ref,"fdr",sep="."))
	#print(test)
	#print(names(stat.table))
	colnames(stat.table) <- c(annotation.names,
								paste(trt,"avg.beta",sep="."),
								paste(ref,"avg.beta",sep="."),
								paste(trt,"vs",ref,"delta.beta",sep="."),
								paste(trt,"vs",ref,"pvalue",sep="."),
								paste(trt,"vs",ref,"fdr",sep="."))
	#print(names(stat.table))
	} else if(length(groups) > 2) {
	print("Calculating Average Beta and ANOVA p-values")
	col.names <- c(annotation.names)

	for (i in 1:length(groups))
	{
		temp.grp <- groups[i]
		grp.beta <- beta.values[,which(sample.group == temp.grp)]
		avg.beta <- apply(grp.beta, 1, custom.mean)
		col.names <- col.names <- c(col.names, paste(temp.grp,"avg.beta",sep="."))
		stat.table <- data.frame(stat.table, avg.beta=avg.beta)
		colnames(stat.table) <- col.names
	}#end for (i in 1:length(groups)

	if(paired)
		{
			print("Factor in Paired Samples")
			pvalue <- apply(beta.values, 1, annova.2way.pvalue, grp.levels=sample.group, pairing.levels=pairing.group)
		}#end if(paired)
	else
		{
			pvalue <- apply(beta.values, 1, annova.pvalue, grp.levels = sample.group)
		}#end else

	anova.fdr <- p.adjust(pvalue, method="fdr")
	col.names <- c(col.names, "annova.pvalue", "annova.fdr")
	stat.table <- data.frame(stat.table, anova.pvalue = pvalue, anova.fdr = anova.fdr)
	colnames(stat.table) <- col.names
	} else {
	print("Invalid groups specified in sample description file!")
	}

	stat.table <- stat.table[order(stat.table$Chr, stat.table$Loc), ]
	if(output.format == 'xls'){
		xlsfile <- file.path(data.folder, paste(project.name,"_CpG_site_stats.xlsx",sep=""))
		WriteXLS("stat.table", ExcelFileName = xlsfile)
	} else if(output.format == 'txt') {
		txtfile <- file.path(data.folder, paste(project.name,"_CpG_site_stats.txt",sep=""))
		write.table(stat.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		print(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}

	if((create.wig == "avg")|(create.wig == "avg.and.sample")){
			#create .wig file based upon average methylation values
			if(ref == "continuous"){
				print("Not creating average .wig files for continuous variable")
			} else{
				wig.folder<-file.path(site.folder,paste(project.name,"wig",sep="_"))
				dir.create(wig.folder, showWarnings=FALSE)

				temp.stat.file <- file.path(wig.folder, "temp.txt")
				Perl.Path <- file.path(path.package("COHCAP"), "Perl")
				perl.script <- file.path(Perl.Path , "create_wig_files.pl")
				write.table(stat.table, temp.stat.file, quote=F, row.names=F, sep="\t")
				cmd <- paste("perl \"",perl.script,"\" \"", temp.stat.file,"\" \"", wig.folder,"\"", sep="")
				res <- system(cmd, intern=TRUE, wait=TRUE)
				message(res)
			}
	}#end if(create.wig)
	
	#filter sites
	print(dim(stat.table))
	filter.table <- stat.table
	if(length(groups) == 1) {
		temp.avg.beta <- stat.table$avg.beta
		filter.table <- filter.table[(temp.avg.beta >= methyl.cutoff) | (temp.avg.beta <=unmethyl.cutoff),]
	} else if((length(groups) == 2)|(ref == "continuous")) {
		temp.trt.beta <- stat.table[[6]]
		temp.ref.beta <- stat.table[[7]]
		temp.delta.beta <- abs(stat.table[[8]])
		temp.pvalue <- stat.table[[9]]
		temp.fdr <- stat.table[[10]]
		
		if(unmethyl.cutoff > methyl.cutoff)
			{
				print("unmethyl.cutoff > methyl.cutoff ...")
				print("so, delta-beta decides which group should be subject to which cutoff")
				temp.delta.beta <- stat.table[[8]]
				temp.methyl.up <- filter.table[(temp.trt.beta >= methyl.cutoff) & (temp.ref.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta >= delta.beta.cutoff),]
				temp.methyl.down <- filter.table[(temp.ref.beta >= methyl.cutoff) & (temp.trt.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta <= -delta.beta.cutoff),]
				filter.table <- rbind(temp.methyl.up, temp.methyl.down)
			}
		else
			{
				temp.methyl.up <- filter.table[(temp.trt.beta >= methyl.cutoff) & (temp.ref.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta >= delta.beta.cutoff),]
				temp.methyl.down <- filter.table[(temp.ref.beta >= methyl.cutoff) & (temp.trt.beta<=unmethyl.cutoff) & (temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff) & (temp.delta.beta >= delta.beta.cutoff),]
				filter.table <- rbind(temp.methyl.up, temp.methyl.down)
			}
	} else if(length(groups) > 2) {
		temp.pvalue <- stat.table[[ncol(stat.table) -1]]
		temp.fdr <- stat.table[[ncol(stat.table)]]
		filter.table <- filter.table[(temp.pvalue <= pvalue.cutoff) & (temp.fdr <= fdr.cutoff),]
	} else {
	print("Invalid groups specified in sample description file!")
	}
	print(dim(filter.table))
	filter.table = filter.table[!is.na(filter.table$SiteID),]
	print(dim(filter.table))

	if(output.format == 'xls'){
		xlsfile <- file.path(site.folder, paste(project.name,"_CpG_site_filter.xlsx",sep=""))
		WriteXLS("filter.table", ExcelFileName = xlsfile)
	} else if(output.format == 'txt') {
		txtfile <- file.path(site.folder, paste(project.name,"_CpG_site_filter.txt",sep=""))
		write.table(filter.table, file=txtfile, quote=F, row.names=F, sep="\t")
	} else {
		warning(paste(output.format," is not a valid output format!  Please use 'txt' or 'xls'.",sep=""))
	}
	
	return(filter.table)
}#end def COHCAP.site
