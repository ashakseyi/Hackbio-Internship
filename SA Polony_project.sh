#Project Code Task
#!/usr/bin/env bash
set -euo pipefail
# sa_polony_amr_pipeline.sh
# Pipeline to:
# 1. download sample list (100 isolates),
# 2. QC, trim, (optionally) assemble,
# 3. species id via BLAST,
# 4. AMR detection via abricate (resfinder, card, ncbi) and virulence (vfdb),
# 5. summarize results for a linkedIn post.
#
# Usage: bash sa_polony_amr_pipeline.sh
# -----------------------
# USER CONFIG
# -----------------------
#Make a directory named polony__raw_reads
mkdir polony_raw_reads
#change working directory to polony_raw_reads
#Transfer data from a server and follow redirection using -L then give an output script
curl -L 'https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh' -o Polony_100_download.sh
#Use screen to run download in background or a session so it survives logout
screen -S myproject
#Use bash to run the script
bash Polony_100_download.sh
#Change working directory
cd polony_raw_reads
#Confirm download using ls -1 *.gz | wc -l 
ls -1 *.gz | wc -l

# -----------------------
# Quality Control
# -----------------------
#Make a directory polony_qc_raw_reads
mkdir -p polony_qc_raw_reads
#Use Screen
screen -S myproject
#This code selects all FASTQ files using *.fastq.gz and outputs it to polony_qc_raw_reads using 8 CPU threads to run multiple jobs in parallel
fastqc polony_raw_reads/*.fastq.gz -o polony_qc_raw_reads -t 8
#Change working directory
cd polony_qc_raw_reads
#Confirm using ls -1 *.html | wc -l to count
ls -1 *.html | wc -l
# -----------------------
# Trimming using fastp
# -----------------------
#Make a directory trimmed_reads 
mkdir -p trimmed_reads
#Write a script for the trimming using fastp
nano trim.sh
#The scriptgoes through every file in polony_raw_reads/ ending with _1.fastq.gz (the forward read) and define matching reverse reads
#The code use basename to call out read name without directory name after which fastp is executed
for r1 in polony_raw_reads/*_1.fastq.gz
do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}   # find matching R2 file
    base=$(basename "$r1" _1.fastq.gz) # sample name without _1.fastq.gz

    fastp \
      -i "$r1" \
      -I "$r2" \
      -o trimmed_reads/${base}_1.trimmed.fastq.gz \
      -O trimmed_reads/${base}_2.trimmed.fastq.gz \
      -h trimmed_reads/${base}_fastp.html \
      -j trimmed_reads/${base}_fastp.json \
      -w 8
done 

#Use ls -1 *fastp.html | wc -l to count the number of fastp.html file

#....................................
# QC trimmed reads
#....................................
#Make a directory qc_trimmed_reads
mkdir qc_trimmed_reads
#Use Screen
screen -S qc_trimmed_reads
#This code selects all FASTQ files using *.trimmed.fastq.gz and outputs it to qc_trimmed_reads using 8 CPU threads to run multiple jobs in parallel
fastqc trimmed_reads/*.trimmed.fastq.gz -o qc_trimmed_reads -t 8
#Change working directory
cd qc_trimmed_reads
#Confirm using ls -1 *.html | wc -l to count
ls -1 *.html | wc -l

# -----------------------
# Genome assembly with SPAdes
# -----------------------
echo "[`date`] assembling with SPAdes (may take long). Will run only if spades.py is available."
if command -v spades.py >/dev/null 2>&1; then
  mkdir -p assembly
  # assemble paired samples if present
  for r1 in trim/*_1.trim.fastq.gz; do
    [[ -e "$r1" ]] || continue
    sample=$(basename "$r1" | sed -E 's/_1.trim.fastq.gz$//')
    r2="trim/${sample}_2.trim.fastq.gz"
    outdir="assembly/${sample}"
    if [[ ! -d "${outdir}" ]]; then
      spades.py -1 "${r1}" -2 "${r2}" -o "${outdir}" --threads "${THREADS}" --careful > logs/spades_${sample}.log 2>&1 || true
    fi
  done
  # single-end assemblies (if any)
  for se in trim/*.trim.fastq.gz; do
    [[ -e "$se" ]] || continue
    sample=$(basename "$se" | sed -E 's/.trim.fastq.gz$//')
    outdir="assembly/${sample}"
    if [[ ! -d "${outdir}" ]]; then
      spades.py -s "${se}" -o "${outdir}" --threads "${THREADS}" --careful > logs/spades_${sample}.log 2>&1 || true
    fi
  done
else
  echo "spades.py not found: skipping assembly step."
fi

# -----------------------
# 5) Species identification using BLAST
#  	Sequence ID					Blast result  
1 	SRR27013147					Listeria monocytogenes
2	SRR27013245					Brucella anthropic
3	SRR27013246					Listeria monocytogenes
4	SRR27013247					Listeria monocytogenes
5	SRR27013248					Listeria monocytogenes
6	SRR27013249					Listeria monocytogenes
7	SRR27013250					Listeria monocytogenes
8	SRR27013251					Listeria monocytogenes
9	SRR27013252					Listeria monocytogenes
10	SRR27013253					Listeria monocytogenes
11	SRR27013254					Listeria monocytogenes
12	SRR27013256					Listeria monocytogenes
13	SRR27013257					Listeria monocytogenes
14	SRR27013258					Listeria monocytogenes
15	SRR27013259					Listeria monocytogenes
16	SRR27013260					Listeria monocytogenes
17	SRR27013261					Listeria monocytogenes
18	SRR27013262					Listeria monocytogenes
19	SRR27013264					Listeria monocytogenes
20	SRR27013265					Listeria monocytogenes
21	SRR27013266					Listeria monocytogenes
22	SRR27013267					Listeria monocytogenes
23	SRR27013268					Listeria monocytogenes						
	
								
# -----------------------
# 6) AMR detection with abricate
# -----------------------
#Use nano to create ABRicate.sh that contains the script below
nano ABRicate.sh
#!/bin/bash
# Run ABRicate on all SPAdes assemblies
# Input: assemblies/*_spades/contigs.fasta
# Output: abricate_results/SAMPLE.tab

ASSEMBLY_DIR="assemblies"
OUT_DIR="abricate_results"

mkdir -p "$OUT_DIR"

for contig in ${ASSEMBLY_DIR}/*_spades/contigs.fasta
do
    SAMPLE=$(basename $(dirname "$contig") _spades)

    echo "Running ABRicate for $SAMPLE ..."
    abricate "$contig" > "${OUT_DIR}/${SAMPLE}_abricate.tab"
done

echo "All ABRicate analyses completed."
# -----------------------
# Resistance Information extraction
# ---------------------------------
#Use nano to create resistance_extract.sh that contains the script below
nano resistance_extract.sh 
#!/bin/bash
# Extract resistance info from ABRicate .tab outputs
# Usage: ./extract_resistance.sh abricate_results/ output_folder/

IN_DIR=${1:-abricate_results}
OUT_DIR=${2:-resistance_reports}

mkdir -p "$OUT_DIR"

# 1. Full report: Sample, Gene, Product, Resistance
echo -e "Sample\tGene\tProduct\tResistance" > "$OUT_DIR/full_resistance_report.txt"
for f in ${IN_DIR}/*.tab
do
    sample=$(basename "$f" _abricate.tab)
    awk -F"\t" -v S="$sample" 'NR > 1 && $6 != "" {print S "\t" $6 "\t" $14 "\t" $15}' "$f" \
        >> "$OUT_DIR/full_resistance_report.txt"
done

# 2. Count occurrences of each resistance class
awk -F"\t" 'NR > 1 && $15 != "" {count[$15]++} END {for (c in count) print c, count[c]}' ${IN_DIR}/*.tab \
    | sort > "$OUT_DIR/resistance_class_counts.txt"

echo "✅ Resistance reports generated in $OUT_DIR/"
echo " - full_resistance_report.txt (sample, gene, product, resistance)"
echo " - resistance_class_counts.txt"

# --------------------------------------------------
#Summarizes AMR prevalence and resistance implications
# --------------------------------------------------

# -----------------------
# 7) Toxin gene detection (hly plcA plcB)
# -----------------------
# The literature documents the presence of hly gene responsible for listeriolysin O (LLO), a crucial virulence factor that forms pores in host cell membranes, enabling the bacteria to escape from phagosomes into the host cell cytoplasm and cause disease (listeriosis) 
# L. monocytogenes also contains plcA gene that produces phosphatidylinositol-specific phospholipase C (PI-PLC), and the plcB gene encodes a broad-range phospholipase C (PC-PLC). Both are crucial virulence factors, acting alongside listeriolysin O (LLO) to promote the rupture of phagocytic vacuoles and facilitate the bacteria's escape into the host cell's cytoplasm.
# Brucella anthropic does not possess any toxin gene 

# -----------------------
# 8) Summarize AMR prevalence across isolates
# -----------------------
# All the Listeria monocytogens isolates carry FOSFOMYCIN and LINCOSAMIDE resistance genes
# The Brucella anthropic isolate (SRR27013245) carries the CEPHALOSPORIN resistance gene in addition to FOSFOMYCIN and LINCOSAMIDE resistance genes
# The prevalence of FOSFOMYCIN and LINCOSAMIDE resistance genes is 100% across all isolates (Listeria monocytogenes and Brucella anthropica)
# The prevalence of CEPHALOSPORIN resistance gene was 2% (1/50) among all the isolates, but 0% among Listeria monocytogenes 
# ----------------------
# Recommended Alternative Drugs
# ----------------------
# For Listeria monocytogenes
# Aminopenicillins (Ampicillin/Amoxicillin) + Aminoglycosides (Gentamicin): This synergistic combination is the reference treatment for listeriosis and is bactericidal. It should be used at high doses and for a prolonged duration due to the intracellular nature of L. monocytogenes. 
# Trimethoprim-Sulfamethoxazole: This is a common alternative, particularly for patients with a penicillin allergy. 
# Vancomycin: This glycopeptide antibiotic is a potential alternative treatment option. 
# Rifampin: Rifampicin has also been proposed as a possible alternative treatment. 
# Linezolid: The oxazolidinone antibiotic linezolid is another alternative that may be considered.

# For Brucella anthropic
# Doxycycline + Gentamicin: This combination is considered highly effective and cost-effective for brucellosis treatment and is a strong alternative. 
# Doxycycline + Rifampicin: This is a traditional and effective option for brucellosis. 
# Doxycycline + Streptomycin: A standard treatment regimen for brucellosis.
# Tigecycline: A newer antibiotic with promising in vitro activity against Brucella and potential for monotherapy or shorter treatment durations. 
# Trimethoprim-sulfamethoxazole (TMP-SMZ): Can be used in certain cases, such as in children or with other antibiotics like gentamicin. 
# Fluoroquinolones: While showing in vitro activity against Brucella, their clinical use is limited, particularly in children and in combination with other drugs due to resistance concerns and potential for increased relapse rates. 
# -----------------------
# 9) Report with methods, results, and public health discussion
# -----------------------
# AMR & Outbreak analysis — South Africa polony outbreak (2017–2018)
## Project goal
- Obtain the sequences
- Do a quality check using fastqc on raw sequence reads
- Trim the sequencing adapters using fastp
- Do a second fastqc quality check on the trimmed reads
- Perform a denovo assembly using SPAdes.py on trimmed reads
- Identify organism(s) from the output of SPAdes from Bandage tool using BLAST.
- Characterize AMR gene profile across isolates using ABRICATE.
- Summarize results and propose treatment options considering resistant isolates.
- Identify from literature the presence of toxin genes (hly, plcA, plcB).
- Provide reproducible code.

## Methods
# Bash script: sa_polony_amr_pipeline.sh was used for this analysis which contained several other scripts
- Data download using provided script.
- QC: FastQC. Trimming: fastp.
- Assembly: SPAdes.
- Species ID: BLAST of largest contigs (blastn).
- AMR detection: abricate 
- Toxin detection: This was checked from the literature

## Key outputs (paths)
# After running the pipeline (sa_polony_amr_pipeline.sh), here are the important directories and files:
# Raw reads (downloaded):
polony_raw_reads/*.fastq.gz
# QC reports (raw):
polony_qc_raw_reads/*.html
# Trimmed reads:
trimmed_reads/*_trimmed.fastq.gz
# QC reports (trimmed):
qc_trimmed_reads/*.html
# Genome assemblies:
assembly/<sample>/contigs.fasta
(SPAdes output, input for BLAST/abricate)
# Species identification (manual/BLAST results):
Tabular list of isolate IDs vs identified species
# AMR detection (abricate outputs):
abricate_results/<sample>_abricate.tab
# Resistance summary reports:
resistance_reports/full_resistance_report.txt
resistance_reports/resistance_class_counts.txt
# Toxin/virulence results (hly, plcA, plcB):
summary of result from literature search 

## Results (fill in after pipeline finishes)
Species Identification
Majority: Listeria monocytogenes (the known outbreak agent).
Minority: Brucella anthropic (1 isolate, unexpected finding).

AMR Profile (from abricate):
Listeria monocytogenes:
100% carried fosfomycin and lincosamide resistance genes.
Brucella anthropic:
Carried fosfomycin, lincosamide, plus an additional cephalosporin resistance gene.

Prevalence summary:
FOSFOMYCIN = 100% (all isolates)
LINCOSAMIDE = 100% (all isolates)
CEPHALOSPORIN = 2% (only the Brucella isolate)
Toxin Genes (literature-supported):

Listeria monocytogenes:
hly (listeriolysin O, pore-forming toxin)
plcA (PI-PLC)
plcB (PC-PLC)

Brucella anthropic:
No hly/plc toxins documented.

## Public Health Discussion
The 2017–2018 South Africa polony outbreak was the world’s largest recorded Listeria outbreak, linked to contaminated ready-to-eat meat. This pipeline shows how genomic surveillance can confirm the pathogen and characterize resistance/virulence factors.

AMR implications:
High prevalence (100%) of fosfomycin and lincosamide resistance suggests these are not viable treatment options.
Fortunately, standard first-line therapies for listeriosis (ampicillin + gentamicin, or TMP-SMX for allergic patients) remain effective since no resistance genes for these were detected.
For Brucella anthropic, cephalosporin resistance is concerning, but traditional regimens (doxycycline + gentamicin/rifampicin) are still valid.

Virulence implications:
Detection of Listeria virulence genes (hly, plcA, plcB) reinforces their role in intracellular survival and disease severity.
This highlights why Listeria infections are particularly dangerous in pregnant women, neonates, and immunocompromised individuals.

Public health message:
This project underscores the power of bioinformatics pipelines in real-time outbreak investigation.
By integrating QC, assembly, AMR profiling, and virulence factor detection, we can rapidly assess both treatment risks and transmission dynamics.

Such genomic surveillance can guide outbreak response, antibiotic stewardship, and food safety interventions.
## Interpretation & suggested antibiotics
- If organism is *Listeria monocytogenes*: standard therapy often includes **ampicillin +/- gentamicin** for severe invasive disease (neonatal/meningitis). If beta-lactam resistance markers are absent, ampicillin remains appropriate. If aminoglycoside resistance genes are present, gentamicin addition may be reconsidered.
- If we observe genes encoding extended-spectrum beta-lactamases (ESBLs)/AmpC, avoid penicillins and cephalosporins; consider agents based on susceptibility (e.g., carbapenems for serious infections resistant to others), but clinical decisions must be guided by phenotypic AST and local guidelines.
- NOTE: This pipeline provides genomic predictions. Clinical treatment requires clinical microbiology and physician judgement.

## Next steps for submission
- Push scripts and REPORT.md to a GitHub repo.
- Create a clear README with run instructions and resource requirements 
- Create a LinkedIn long-form post describing the project and posting the GitHub link

