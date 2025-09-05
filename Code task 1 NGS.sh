#Move to your working directory using cd Oluwaseyi_Ashaka
Use one line of command to get the sequence file
wget wget https://github.com/rieseberglab/fastq-examples/raw/refs/heads/master/data/HI.4019.002.index_7.ANN0831_R1.fastq.gz 
#Use ls to check the downloaded file
#Use the fastqc command to do the quality check on the downloaded reads
fastqc HI.4019.002.index_7.ANN0831_R1.fastq.gz
#Use ls to locate the generated .html and .zip report files
ls
#Open another terminal and use the code below to transfer the .html file
scp -r pasteur@135.181.163.242:/home/pasteur/Oluwaseyi_Ashaka/HI.4019.002.index_7.ANN0831_R1_fastqc.html ./
#Download the .html file from the terminal and view it in the browser
