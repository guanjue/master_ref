library(DNAshapeR)
plotDNAshape <- function(input,output1,output2,seqext)
{
	fa_p_shape=getShape(input, shapeType = 'All', parse = TRUE) ## calculate DNA shape score
	pdf(output1)	
	par(mfrow=c(2,2))
	### plot DNA shape scores
		plot(colMeans(fa_p_shape$MGW),main='MGW shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(4.5,5.5)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$MGW),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp
		
		plot(colMeans(fa_p_shape$HelT),main='HelT shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(33,35)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$HelT),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp

		plot(colMeans(fa_p_shape$ProT),main='ProT shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(-9,-7)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$ProT),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp

		plot(colMeans(fa_p_shape$Roll),main='Roll shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(-1.5,-0.5)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$Roll),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp
	dev.off()
	
	png(output2)	
	par(mfrow=c(2,2))
	### plot DNA shape scores
		plot(colMeans(fa_p_shape$MGW),main='MGW shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(4.5,5.5)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$MGW),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp

		plot(colMeans(fa_p_shape$HelT),main='HelT shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(33,35)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$HelT),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp

		plot(colMeans(fa_p_shape$ProT),main='ProT shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(-9,-7)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$ProT),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp

		plot(colMeans(fa_p_shape$Roll),main='Roll shape score',ylab="DNA shape score",xlab='Distance from MEME motif (bp)',xaxt = "n",ylim=c(-1.5,-0.5)) ## plot DNA shape score mean
		lines(colMeans(fa_p_shape$Roll),col='blue',lwd=2) ## draw a line connect DNA shape score means
		axis(1, at=c(seqext-1-50,seqext-1-25,seqext-1,seqext-1+25,seqext-1+50), labels=c(-50,-25,0,25,50)) ## change x axis label to 101 bp
	dev.off()
}
args<-commandArgs(TRUE)
plotDNAshape(toString(args[1]),toString(args[2]),toString(args[3]),as.numeric(args[4]))

### e.g. Rscript plotDNAshape.R test.fa plot.pdf plot.png 52
