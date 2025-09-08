#Code Task 1
#Move to your working directory using cd Oluwaseyi_Ashaka
Use one line of command to get the sequence file
wget wget https://github.com/rieseberglab/fastq-examples/raw/refs/heads/master/data/HI.4019.002.index_7.ANN0831_R1.fastq.gz 
#Use ls to check the downloaded file
#Use fastqc command to do the quality check on the downloaded reads
fastqc HI.4019.002.index_7.ANN0831_R1.fastq.gz
#Use ls to locate the generated .html and .zip report files
ls
#Open another terminal and use the code below to transfer the .html file
scp -r pasteur@135.181.163.242:/home/pasteur/Oluwaseyi_Ashaka/HI.4019.002.index_7.ANN0831_R1_fastqc.html ./
#Download the html file from the terminal and view on the browser

#Code Task 2
#Open the FastQC .html report from Task 1.
#Identify and note:
#Bases with low quality in “Per base sequence quality.”
Base 1-9 in the reads are of low quality
#Any warnings about adapter sequences in “Adapter Content.”
There are no warnings but on the graph there are some adapter sequences after 60 reads
#Any “Overrepresented sequences” flagged.
There are three (3) over-represented sequences 
Two were from Illumina Multiplexing PCR Primer, which are 96% over 25bp and 27bp, while the most over-represented sequence was the TruSeq Adapter, Index 7 which are 100% over 50bp)
#Write a short summary (2–3 sentences) describing the quality of the sequencing data and potential issues.
The Per base sequence quality shows that after 100bp read quality was reduced. The Per sequence quality scores showed that sequence reads 40 bp had good quality score. This means that the reads should be between 40 and 100 bp. The per-base sequence content showed that sequence lower than 10bp had inconsistencies. The overrepresented sequences were sequencing primers and adapters. 
Conclusion: The read should be trimmed and sequencing may be repeated with reads between 40 and 50bp
