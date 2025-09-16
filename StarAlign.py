#### Script created by Dr. Reema Singh
#### Contact: Reema Singh (Email: res498@usask.ca)

import os
import sys

##### Use following argument to run the program from command line

## python StarAlign.py Input.txt HumanGenome_STARIndex FASTQ  Alignment

Input = sys.argv[1]
GenomeIndex = sys.argv[2]
Fastq = sys.argv[3]
Output = sys.argv[4]

with open(Input) as acces:
    for text in acces:
        numb = text.split("\t")[0]
        numb1 = text.split("\t")[1]
        numb2 = text.split("\t")[2]
        
        out1 = numb.rstrip() 
        out2 = out1+'Aligned.sortedByCoord.out.bam'

        cmd = "STAR --runThreadN 20 --genomeDir" + ' ' + GenomeIndex + ' '+"--readFilesIn" + ' ' + Fastq + "/"+ numb1 + ' ' + Fastq + "/" +numb2.rstrip()+ ' ' +  "--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000 --outBAMsortingThreadN 20 --outFileNamePrefix" + ' ' +Output+"/"+ out1+' '+"--genomeLoad NoSharedMemory --outSAMunmapped Within --outSAMstrandField intronMotif --outSJtype Standard"
        #print(cmd)
        os.system(cmd)
        cmd1 = "samtools index"+ ' '+Output+"/"+out2
        #print(cmd1)
        os.system(cmd1)
        
