#### Script created by Dr. Reema Singh
#### Contact: Reema Singh (Email: res498@usask.ca)

library(DESeq2)
library(vsn)
library(ggplot2)

#setwd("/home/reemas/scratch/EBOLA_DE/DE")

Counts <- read.table(file="readCount.txt",sep="\t", header = TRUE, check.names = FALSE)
names(Counts) <- gsub("Aligned.sortedByCoord.out.bam", "", names(Counts))
StudyDes <- read.table(file="StudyDesign.txt",sep="\t", header = TRUE)

ddsMat <- DESeqDataSetFromMatrix(countData = Counts, colData = StudyDes, design = ~Treatment)
nrow(ddsMat)
keep <- rowSums(counts(ddsMat) >= 15) >= 6
ddsMat <- ddsMat[keep,]
nrow(ddsMat)

vsd <- vst(ddsMat, blind = FALSE)
head(assay(vsd))
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(StudyDes$Treatment, StudyDes$Lane, sep="-")


mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, colData(vsd))
pdf("Figure1_MDSPlot_TreatmentTimeLane.pdf")
p <- qplot(X1,X2,color=Treatment,shape=Lane,data=as.data.frame(mds),xlab = "Dimension1", ylab = "Dimension2")+geom_point(size =3)
p + theme_classic()+theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=9), axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
dev.off()

ddsMat <- DESeq(ddsMat)
write.csv(as.data.frame(counts(ddsMat)),file="EbolaRawCount.csv")

res <- results(ddsMat, alpha = 0.05,contrast=c("Treatment", "ZEBOV.24hr","MOCK.24hr"))
mcols(res, use.names = TRUE)
summary(res)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file="ZEBOV24hrvsMOCK24hr.csv")

res1 <- results(ddsMat, alpha = 0.05, contrast=c("Treatment", "RESTV.24hr", "MOCK.24hr"))
mcols(res1, use.names = TRUE)
summary(res1)
res1Ordered <- res1[order(res1$pvalue),]
write.csv(as.data.frame(res1Ordered), file="RESTV24hrvsMOCK24hr.csv")

res2 <- results(ddsMat, alpha=0.05, contrast=c("Treatment", "ZEBOV.48hr", "MOCK.48hr"))
mcols(res2, use.names = TRUE)
summary(res2)
res2Ordered <- res2[order(res2$pvalue),]
write.csv(as.data.frame(res2Ordered), file="ZEBOV48hrvsMOCK48hr.csv")

res3 <- results(ddsMat, alpha=0.05, contrast=c("Treatment", "RESTV.48hr", "MOCK.48hr"))
mcols(res3, use.names = TRUE)
summary(res3)
res3Ordered <- res3[order(res3$pvalue),]
write.csv(as.data.frame(res3Ordered), file="RESTV48hrvsMOCK48hr.csv")

pdf("Figure2_MA-plot.pdf")
par(mfrow=c(2,2))
plotMA(res, ylim=c(-2,2), cex=.6, alpha=0.05, colSig = "red", main ="A) ZEBOV.24hr vs MOCK.24hr")
abline(h=c(-1,1), col="dodgerblue", lwd=2)
plotMA(res1, ylim=c(-2,2), cex=.6, alpha=0.05, colSig = "red", main ="B) RESTV.24hr vs MOCK.24hr")
abline(h=c(-1,1), col="dodgerblue", lwd=2)
plotMA(res2, ylim=c(-2,2), cex=.6, alpha=0.05, colSig = "red", main ="C) ZEBOV.48hr vs MOCK.48hr")
abline(h=c(-1,1), col="dodgerblue", lwd=2)
plotMA(res3, ylim=c(-2,2), cex=.6, alpha=0.05, colSig = "red", main ="D) RESTV.48hr vs MOCK.48hr")
abline(h=c(-1,1), col="dodgerblue", lwd=2)
dev.off()

