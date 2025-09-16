#### Script created by Dr. Reema Singh
## Contact: Reema Singh (Email: res498@usask.ca)

library(Rsubread)
setwd("/home/reemas/scratch/EBOLA_DE/Alignment")
file_list <- list.files(pattern="*.bam$")

Sample1 <- featureCounts(files=file_list,annot.ext="/home/reemas/scratch/EBOLA_DE/HumanGenome/Homo_sapiens.GRCh38.109.gtf", isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",isPairedEnd=TRUE)

write.table(Sample1$counts, file="readCount.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(Sample1$stat, file="readCountStat.txt",sep="\t",row.names=TRUE,col.names=TRUE,quote=FALSE)


