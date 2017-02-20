library(gplots)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
data1=read.table(args[1],header=F)


my_palette1 = colorRampPalette(c('cornsilk','darkslategray1','darkorchid1','gray','sienna1'))

data1=as.matrix(cbind(data1,data1))


png(paste(args[3],'color_split.png',sep='.'),res=300,width=20,height=800)
heatmap.2(data1, Rowv=NULL, Colv=NULL, col=my_palette1, margins=c(0.01,0.01),symm=F,symkey=F,symbreaks=F, key=F, labRow=F, labCol=F, keysize=0.001, scale="none", density.info="none", trace="none", dendrogram='none')
#heatmap(data1, Rowv=NULL, Colv=NULL, col=my_palette1, margins=c(0.01,0.01), dendrogram='none')
dev.off()
