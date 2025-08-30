#!/bin/bash
# Author: Oluwaseyi Ashaka
# Description: Linux & Bash Fundamentals Project 1 and 2
# Note: Run this script with: bash scriptname.sh

#######################################
# Project 1
#######################################

# Print your name
echo "Oluwaseyi_Ashaka"

# Create a folder titled your name
mkdir -p ~/Oluwaseyi_Ashaka

# Create another new directory titled biocomputing and change to that directory
mkdir -p ~/biocomputing && cd ~/biocomputing || exit

# Download files
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.fna
wget https://raw.githubusercontent.com/josoga2/dataset-repos/main/wildtype.gbk

# Move the .fna file to the folder titled your name
mv wildtype.fna ~/Oluwaseyi_Ashaka/

# Delete duplicate gbk file if it exists
rm -f wildtype.gbk.1

# Confirm if the .fna file is mutant or wild type (tatatata vs tata)
grep "tatatata" wildtype.fna
grep "tata" wildtype.fna

# If mutant, print all matching lines into a new file
grep "tatatata" wildtype.fna > mutant_1.txt

# Count number of lines (excluding header) in the .gbk file
grep -v LOCUS wildtype.gbk | wc -l

# Print the sequence length of the .gbk file
grep "LOCUS" wildtype.gbk

# Print the source organism of the .gbk file
grep "SOURCE" wildtype.gbk

# List all the gene names of the .gbk file
grep "/gene=" wildtype.gbk

# List contents of both directories
ls ~/biocomputing ~/Oluwaseyi_Ashaka


#######################################
# Project 2
#######################################

# Activate your base conda environment (assumes conda is installed)
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

# Create a conda environment named funtools
conda create -n funtools python=3.10 -y

# Activate the funtools environment
conda activate funtools

# Add channels
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install Figlet
conda install -y pyfiglet || sudo apt-get install -y figlet

# Run figlet
figlet "Oluwaseyi Ashaka"

# Install bioinformatics tools
conda install -y bwa
conda install -y blast
conda install -y samtools
conda install -y bedtools
conda install -y spades
conda install -y bcftools
conda install -y fastp

# Create and activate new environment for multiqc
conda create -n myenv python=3.11 -y
conda activate myenv

# Add channels again for bioinformatics
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install multiqc
conda install -y multiqc

echo "All tasks completed successfully!"
