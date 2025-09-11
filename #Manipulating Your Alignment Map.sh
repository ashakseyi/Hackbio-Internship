**#Manipulating Your Alignment Map **
# Use for loop to sort and move sorted files to mapped_sorted
 for bam in alignment_map/*.bam; do     samtools sort "$bam" -o mapped_sorted/$(basename "$bam" .bam)_sorted.bam; done
 cd mapped_sorted/
 ls
#Copy the mapped_sorted file to IGV
 cp mapped_sorted/*_sorted.bam IGV/
 ls IGV/
 cd IGV/
#Index the aligned mapped and sorted file to get a bai file which is like the table of contents 
 samtools index ACBarrie_sorted.bam 
#Use the SCP command on another terminal to move the bai and bam files to your shell google cloud terminal and download to your computer
scp -r pasteur@135.181.163.242:/home/pasteur/Oluwaseyi/raw_reads/IGV/ACBarrie_sorted.bam.bai ./
scp -r pasteur@135.181.163.242:/home/pasteur/Oluwaseyi/raw_reads/IGV/ACBarrie_sorted.bam ./

#Note that the above will require a password
 
#open the webpage https://igv.org/app/ on a web browser
#Look for the reference genome or upload the website
#At the track tab upload both the .bam and .bai files for visualization, and check the link below for the visuals
https://github.com/ashakseyi/Hackbio-Internship/blob/084e6018069f49e6418543f0034b6d05dd0ad5b2/IGV%20Visualization%20mapping%20reads%20to%20a%20genome.png
