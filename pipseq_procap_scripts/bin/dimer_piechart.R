args <- commandArgs(trailingOnly = TRUE)
data=read.table(args[1],header=F)
data1=data[seq(2, dim(data)[1], 2),1]
pie_data=table(as.matrix(data1))
pie_data=pie_data/sum(pie_data)
pie_data=pie_data[order(-pie_data)]

print(rownames(pie_data))

lbls_num=paste(rownames(pie_data),round(pie_data,digits=2),sep=':')

png(args[2])
par(mar=c(1.5,1.5,2.3,1.5))
pie(pie_data, cex=1.2, cex.main=2.0,main=paste("Total number of Genes:", toString(length(data1)),sep=' '), labels = lbls_num )
dev.off()


