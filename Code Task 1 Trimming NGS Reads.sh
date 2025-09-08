Code Task 1 Trimming NGS Reads
#Create a directory called qc_reads in the current working directory where the sequence reads are present
mkdir qc_reads
#Define an array of sample names: ACBarrie, Alsen, Baxter, Chara, Drysdale.
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

#Use nano trim.sh to create the bash script
#!/bin/bash
#This tells the system that the script should be run with the bash shell.
#Defines a bash array called SAMPLES with the sample names.
#Note that each sample has corresponding paired-end sequencing reads (_R1.fastq.gz and _R2.fastq.gz).
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

for SAMPLE in "${SAMPLES[@]}"; do
#Starts a loop to go through each sample name in the array.
  fastp \
    -i "$PWD/${SAMPLE}_R1.fastq.gz" \
    -I "$PWD/${SAMPLE}_R2.fastq.gz" \
    -o "qc_reads/${SAMPLE}_R1.fastq.gz" \
    -O "qc_reads/${SAMPLE}_R2.fastq.gz" \
    --html "qc_reads/${SAMPLE}_fastp.html" 
#-i → specifies the input file for read 1 (forward reads).
#-I (capital i) → specifies the input file for read 2 (reverse reads) in paired-end sequencing.
#-o → specifies the output file for cleaned read 1.
#-O (capital o) → specifies the output file for cleaned read 2.
#--html → tells fastp to generate an HTML report.
done
#Ends the loop.

#Exit the nano script and use cd to move the terminal to the current working directory where the raw sequences can be found
#Also copy the trim.sh to the current working directory
cd raw_reads
cp trim.sh raw_reads/

bash trim.sh
#This will execute the command 
#Use ls qc_reads to check the html files created
#Note that the output of fastp can be subjected to fastqc, and the report can be checked with multiqc 
