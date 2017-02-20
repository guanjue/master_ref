library(gplots)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
print(length(args))
for (i in seq(1,length(args))) {
	print(i)
	print(args[i])
	data1=read.table(args[i],header=F)
	data1=as.matrix(data1)
	png(paste(args[i],'hist.png',sep='.'))
	plot(density(data1),xlim=c(-250,250))#,ylim=c(0,0.01))
	#heatmap.2(data1, Rowv=NULL, Colv=NULL, col=my_palette1, margins=c(0.01,0.01),symm=F,symkey=F,symbreaks=F, key=F, labRow=F, labCol=F, keysize=0.001, scale="none", density.info="none", trace="none", dendrogram='none')
	#heatmap(data1, Rowv=NULL, Colv=NULL, col=my_palette1, margins=c(0.01,0.01), dendrogram='none')
	dev.off()
}

