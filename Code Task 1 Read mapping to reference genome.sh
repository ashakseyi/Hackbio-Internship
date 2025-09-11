
#**Map Trimmed Reads to a Reference Genome**
#Download the reference genome (reference.fasts) into a folder named references
#Use nano to build an aligner script using nano aligner.sh
#Note that the aligner,sh must be outside with the other folder references and trimmed reads to be used for this alignment

#!/bin/bash
# Call out SAMPLES which contains all the names of the sequences using ‘ ‘ in a bracket ()
SAMPLES=(
  "ACBarrie"
  "Alsen"
  "Baxter"
  "Chara"
  "Drysdale"
)

# Index the reference genome from using the path in the current working directory
bwa index references/reference.fasta
# Create directories: repaired and alignment_map
mkdir repaired
mkdir alignment_map

#Use a for loop to call out SAMPLE in SAMPLES where [@] call out all the SAMPLE in {SAMPLES}
for SAMPLE in "${SAMPLES[@]}"; do
#Repair the trimmed sequences for each SAMPLE specifying input and output reads for both paired reads with a single optional output all into the repaired folder
    repair.sh in1="trimmed_reads/${SAMPLE}_R1.fastq.gz" in2="trimmed_reads/${SAMPLE}_R2.fastq.gz" out1="repaired/${SAMPLE}_R1_rep.fastq.gz" out2="repaired/${SAMPLE}_R2_rep.fastq.gz" outsingle="repaired/${SAMPLE}_single.fq"
#Call out present working directory then perform alignment_mapping using repaired sequences in the repaired folder and converting it into bam files
    echo $PWD
    bwa mem -t 1 \
    references/reference.fasta \
    "repaired/${SAMPLE}_R1_rep.fastq.gz" "repaired/${SAMPLE}_R2_rep.fastq.gz" \
  | samtools view -b \
  > "alignment_map/${SAMPLE}.bam"
#This ends the loop
done

#Run the alignment of the multiple sequences using bash aligner.sh
bash aligner.sh
#Confirm the file sizes using ls -lh for the repaired folder and the alignment_map folder
