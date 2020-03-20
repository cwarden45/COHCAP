colLab <- function(n, labelColors, clusMember) { 
   if(is.leaf(n)) { 
       a <- attributes(n) 
	   #print(a)
       # clusMember - vector of sample names (ordered to match label color.palette)
       # labelColors - a vector of color.palette for the above grouping 
       labCol <- labelColors[clusMember == a$label]
	   #print(labCol)
       attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol) 
   } 
   n 
}

`COHCAP.qc` <-function (sample.file, beta.table, project.name, project.folder, plot.legend=TRUE, color.palette = c("red","blue","green","orange","purple","cyan","pink","maroon","yellow","grey","black",colors()))
{
qc.folder<-file.path(project.folder,"QC")
dir.create(qc.folder, showWarnings=FALSE)

sample.table <- read.table(sample.file, header=F, sep = "\t", stringsAsFactors=TRUE)
samples <- as.character(sample.table[[1]])
for (i in 1:length(samples))
	{
		if(length(grep("^[0-9]",samples[i])) > 0)
			{
				samples[i] <- paste("X",samples[i],sep="")
			}#end if(length(grep("^[0-9]",samples[i])) > 0)
	}#end def for (i in 1:length(samples))
sample.group <- sample.table[[2]]

sample.names <- names(beta.table)[6:ncol(beta.table)]
beta.values <- beta.table[,6:ncol(beta.table)]
print(dim(beta.values))
methyl.max <- ceiling(max(beta.values))
methyl.min <- ceiling(min(beta.values))

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
rm(beta.table)

groups <- levels(sample.group)
print(paste("Group: ",groups,sep=""))
color.palette <- color.palette[1:length(groups)]
print(paste("Color: ",color.palette[1:length(groups)], sep=""))

#sample histogram
print("Calculating Sample Statistics...")
q0 <- array(dim=length(samples))
q25 <- array(dim=length(samples))
q50 <- array(dim=length(samples))
q75 <- array(dim=length(samples))
q100 <- array(dim=length(samples))

print("Creating Sample Histogram...")
hist.file <- file.path(qc.folder, paste(project.name,"_hist.pdf",sep=""))
pdf(file = hist.file)
#print(samples)
#print(length(samples))
for (i in 1:length(samples))
	{
		#print(i)
		#print(paste("Working on Density Distribution for ",samples[i],sep=""))
		data <- -1
		if(length(samples)>1)
			{
				data <- as.numeric(t(beta.values[i]))
			}
		else
			{
				data <- as.numeric(beta.values)
			}
		#print(dim(data))
		quant <- quantile(data, na.rm=T)
		q0[i] <- quant[1]
		q25[i] <- quant[2]
		q50[i] <- quant[3]
		q75[i] <- quant[4]
		q100[i] <- quant[5]
		
		col <- "black"
		if(typeof(sample.group) != "double")
			{
				expression.group <- sample.group[i]
				#print(expression.group)
				for (j in 1:length(groups))
					{
						if(expression.group == groups[j])
							{
							col = color.palette[j]
							}
					}
			}#end if(typeof(sample.group) != "double")
		
		if(i == 1)
			{
				den <- density(data, na.rm=T)
				expr <- den$x
				freq <- den$y
				plot(expr, freq, type="l", xlab = "Beta / Percentage Methylation", ylab = "Density", ylim=c(0,8), col=col)
				if(plot.legend)
					{
						legend("topright",legend=groups,col=color.palette, lwd=3)
					}
			}#end if(i == 1)
		else
			{
				den <- density(data, na.rm=T)
				expr <- den$x
				freq <- den$y
				lines(expr, freq, type = "l", col=col)
			}#end else
	}#end for (i in 1:length(bed.indices))
dev.off()
#print(samples)
#print(q50)
hist.table <- data.frame(sample = samples, min=q0, bottom25=q25, median=q50, top25=q75, max=q100)
hist.text.file <- file.path(qc.folder,paste(project.name,"_descriptive_statistics.txt",sep=""))
write.table(hist.table, hist.text.file, quote=F, row.names=F, sep="\t")
rm(hist.table)

if(length(samples) > 1)
	{
		dist1 <- dist(as.matrix(t(beta.values)))
		clusMember <- sample.group
		labelColors <- as.character(clusMember)
		for (i in 1:length(groups))
			{
			labelColors[clusMember == groups[i]] = color.palette[i]
			}
		hc <- hclust(dist1)
		rm(dist1)
		dend1 <- as.dendrogram(hc)
		rm(hc)
		cluster.file <- file.path(qc.folder,paste(project.name,"_cluster.pdf",sep=""))
		print("Creating Dendrogram...")
		pdf(file = cluster.file)
		#print(labelColors)
		#print(clusMember)
		dend1 <- dendrapply(dend1, colLab, labelColors=labelColors, clusMember=samples) 
		a <- attributes(dend1) 
		attr(dend1, "nodePar") <- c(a$nodePar, lab.col = labelColors) 
		op <- par(mar = par("mar") + c(0,0,0,10)) 
		plot(dend1, horiz=T)
		par(op) 
		dev.off()
		rm(dend1)
		
		#PCA
		print("Creating PCA Plot...")
		#print(dim(beta.values))
		#print(dim(na.omit(data.matrix(beta.values))))
		#print(data.matrix(beta.values)[1,])
		#print(na.omit(data.matrix(beta.values)[1,]))
		pca.values <- prcomp(na.omit(data.matrix(beta.values)))
		#print(attributes(pca.values))
		#print(pca.values$rotation)
		pc.values <- data.frame(pca.values$rotation)
		variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
		pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
		pca.text.file <- file.path(qc.folder,paste(project.name,"_pca.txt",sep=""))
		write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")
		pca.file <- file.path(qc.folder,paste(project.name,"_pca.pdf",sep=""))
		pdf(file=pca.file)
		plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19)
		if(plot.legend)
			{
				legend("topright",legend=groups,col=color.palette,  pch=19)
			}
		dev.off()
	}#end 

}#end def RNA.qc