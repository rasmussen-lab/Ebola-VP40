#### Script created by Dr. Reema Singh
##Contact: Reema Singh (Email: res498@usask.ca)

library(DESeq2)
library(vsn)
library(ggplot2)
library(tidyverse)

Counts <- read.table(file="readCount.txt",sep="\t", header = TRUE, check.names = FALSE)
names(Counts) <- gsub("Aligned.sortedByCoord.out.bam", "", names(Counts))
StudyDes <- read.table(file="StudyDesign.txt",sep="\t", header = TRUE)
ddsMat <- DESeqDataSetFromMatrix(countData = Counts, colData = StudyDes, design = ~Treatment)
keep <- rowSums(counts(ddsMat) >= 15) >= 6
ddsMat <- ddsMat[keep,]

###### Differential Expression

ddsMat <- DESeq(ddsMat)
res <- results(ddsMat, alpha = 0.05,contrast=c("Treatment", "ZEBOV.24hr","MOCK.24hr"))
mcols(res, use.names = TRUE)
summary(res)
resSigUp <- subset(res, res$padj < 0.01 & res$log2FoldChange > 1.5)
up <- rownames(resSigUp)
resSigDown <- subset(res, res$padj < 0.01 & res$log2FoldChange > -1.5)
down <- rownames(resSigDown)
DEG <- c(down,up)

res1 <- results(ddsMat, alpha = 0.05, contrast=c("Treatment", "RESTV.24hr", "MOCK.24hr"))
mcols(res1, use.names = TRUE)
summary(res1)
resSigUp1 <- subset(res1, res1$padj < 0.01 & res1$log2FoldChange > 1.5)
up1 <- rownames(resSigUp1)
resSigDown1 <- subset(res1, res1$padj < 0.01 & res1$log2FoldChange > -1.5) 
down1 <- rownames(resSigDown1)
DEG1 <- c(down1,up1)

res2 <- results(ddsMat, alpha=0.05, contrast=c("Treatment", "ZEBOV.48hr", "MOCK.48hr"))
mcols(res2, use.names = TRUE)
summary(res2)
resSigUp2 <- subset(res2, res2$padj < 0.01 & res2$log2FoldChange > 1.5)
resSigDown2 <- subset(res2, res2$padj < 0.01 & res2$log2FoldChange > -1.5)
up2 <- rownames(resSigUp2)
down2 <- rownames(resSigDown2)
DEG2 <- c(down2,up2)

res3 <- results(ddsMat, alpha=0.05, contrast=c("Treatment", "RESTV.48hr", "MOCK.48hr"))
mcols(res3, use.names = TRUE)
summary(res3)
resSigUp3 <- subset(res3, res3$padj < 0.01 & res3$log2FoldChange > 1.5)
resSigDown3 <- subset(res3, res3$padj < 0.01 & res3$log2FoldChange > -1.5)
up3 <- rownames(resSigUp3)
down3 <- rownames(resSigDown3)
DEG3 <- c(down3,up3)

#### Union of DEGs

T24h <- union(DEG,DEG1)
T48h <- union(DEG2,DEG3)
DEGs <- union(T24h,T48h)

#### MDS plot fo DEGs

ddsMatDEGs <- ddsMat[DEGs,]

vsd <- vst(ddsMatDEGs, blind = FALSE)
head(assay(vsd))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(StudyDes$Treatment, StudyDes$Lane, sep="-")

mds <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mds, colData(vsd))

### after excluding Mocks

mdsFinal <- mds[ ! mds$Treatment == "MOCK.24hr", ]
mdsFinal <- mdsFinal[ ! mdsFinal$Treatment == "MOCK.48hr", ]
MDS <- mdsFinal %>% group_by(Treatment) %>% mutate(chull_val = list(chull(X1,X2)))

pdf("Figure3_MDS_DEGs_WithoutMock.pdf")

ggplot(MDS, aes(X1,X2))+ geom_path(aes(colour = Treatment))+ geom_point(aes(colour = Treatment))+geom_polygon(data =mdsFinal, alpha = 0.2,fill=NA)+ theme_classic()+theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),legend.title = element_blank(),legend.position="top", strip.text = element_text(face="bold", size=9), axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=10))+xlab("Dimension1")+ylab("Dimension2")+scale_color_manual(values = c("darkred","darkgreen","darkblue","purple"))

dev.off()
