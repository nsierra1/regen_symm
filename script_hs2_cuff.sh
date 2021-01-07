# PIPELINE FOR PROCESSING PAIRED END RNA-seq DATA USING HISAT2(MAPPING) 
# AND CUFFLINKS(ASSEMBLY AND QUANTIFICATION)
# Noemie Sierra, July 2019	
	
	# HOME (~) FOLDER CONTAINS: 
	  # PROJECT FOLDER
		# fastq_files - containing PE .fastq.gz data - label_R1 & label_R2
		# slurm-logs - for stdout/stderr output
	  # DATA FOLDER
		# GENOMES and other reference info; all indices and 
	 	# other reference file manipulations are saved here
		
	# Each of the following was run as a separate script on the Peloton Cluster (CSE, UC Davis, 2019)
	# using the following slurm queue submission:
		# sbatch -p (med/high) -n (cores) script
		

  #########################
  # DOWNLOAD THE RAW DATA #
  #########################
  
   
#!/bin/bash
#SBATCH -D ~/
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=#enter_email_here
#SBATCH -o ~/projects/regen_symm/slurm-logs/datadown-o%j.txt
#SBATCH -e ~/projects/regen_symm/slurm-logs/datadown-e%j.txt
#SBATCH -J datadown
#SBATCH -t 04:00:00

cd ~/projects

# Make a new project directory
mkdir regen_symm
cd regen_symm
mkdir slurm-logs
mkdir raw
cd raw

wget -r --user="gec" --password="aardvark dryer rummage" https://jumpgate.caltech.edu/runfolders/volvox02/181214_SN787_0929_AHTG2TBCX2/Unaligned/

#alternative for bulk url download is to save in a textfile and use -i url_list.txt

# Get them out of the nested directory structure and make sure they were not corrupted during the copy process
mv jumpgate.caltech.edu/runfolders/volvox02/181214_SN787_0929_AHTG2TBCX2/Unaligned/Project_*/Sample*/*.fastq.gz .

shasum -a 256 jumpgate.caltech.edu/runfolders/volvox02/181214_SN787_0929_AHTG2TBCX2/Unaligned/Project_*/Sample*/*.fastq.gz >> nested_loc.txt
shasum -a 256 ./*.fastq.gz >> final_loc.txt

# Cat together lane data into a single fastq file
for dir in 210*/; do 
  cd $dir
  bn=`basename $dir` 
  cat *R1_0??.fastq.gz > ~/projects/regen_symm/fastq_files/${bn}_R1.fastq.gz
  cat *R2_0??.fastq.gz > ~/projects/regen_symm/fastq_files/${bn}_R2.fastq.gz
  cd ..
done




  #################################
  # RUN HISAT TO MAP PE FRAGMENTS #
  #################################


#!/bin/bash -l

#SBATCH -D ~/projects/regen_symm/
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=#enter_email_here 
#SBATCH -o ~/projects/regen_symm/slurm-logs/hs2map-o%j.txt
#SBATCH -e ~/projects/regen_symm/slurm-logs/hs2map-e%j.txt
#SBATCH -J hs2map
#SBATCH -t 60:00:00

PATH=$PATH:~/tools/hisat2-2.1.0/hisat2

# Index the genome
cd ~/data
mkdir hisat2_genomeindex
hisat2-build Aurelia.Genome_v1.2_11-27-18.fasta hisat2_genomeindex/Aurelia.Genome_v1.2_11-27-18.index

# Make a sample names list to reference (unique names, without _R1/R2)
cd ~/projects/regen_symm/fastq_files

for s in *.fastq.gz; do
  echo ${s/_R?.fastq.gz} >> ~/projects/regen_symm/samples_a.txt
done
sort ~/projects/regen_symm/samples_a.txt | uniq > ~/projects/regen_symm/samples.txt
rm ~/projects/regen_symm/samples_a.txt

# Run hisat2 - mapping fragments to the genome 
cd ~/projects/regen_symm/
mkdir hs2_sam
mkdir hs2_sam/summaryfiles

while read f; do
  hisat2 --no-softclip \
  -x ~/data/hisat2_genomeindex/Aurelia.Genome_v1.2_11-27-18.index \
  -1 fastq_files/"$f"_R1.fastq.gz \
  -2 fastq_files/"$f"_R2.fastq.gz \
  -S hs2_sam/$f.sam \
  --summary-file hs2_sam/summaryfiles/summaryfile-$f.txt
done < ~/projects/regen_symm/samples.txt

while read f; do
 head -n 2 fastq_files/"$f"_R1.fastq.gz
done < ~/projects/regen_symm/samples.txt


# Save individual summary files to a single file
cat hs2_sam/summaryfiles/summaryfile-*.txt >> ../hs2_sam/SummaryFile_hs2map.txt

# Save a note about new files (usage/generation)
echo -e "The files in regen_symm/hs2_sam/ have been sorted from the hisat2 generated sam files using
cufflinks sort (ID.sam.sorted) for use in cufflinks downstream per the recommendation
of the user manual: 
   http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html 

" | cat > hs2_sam/content_notes.txt



  ########################################################
  # RUN CUFFLINKS TO ASSEMBLE FRAGMENTS INTO TRANSCRIPTS #
  ########################################################


#!/bin/bash -l

#SBATCH -D ~/projects/regen_symm/
#SBATCH --mail-type=ALL
#SBATCH -o ~/projects/regen_symm/slurm-logs/cufflinks-o%j.txt
#SBATCH -e ~/projects/regen_symm/slurm-logs/cufflinks-e%j.txt
#SBATCH --mail-user=#enter_email_here 
#SBATCH -J cufflinks
#SBATCH -t 50:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

cd ~/projects/regen_symm/
mkdir assembled_transcripts

# sort sam files
while read f; do
   sort -k 3,3 -k 4,4n hs2_sam/$f.sam > hs2_sam/$f.sam.sorted
done < ~/projects/regen_symm/samples.txt

# run cufflinks on sorted data
while read f; do
 cufflinks -b ~/data/Aurelia.Genome_v1.2_11-27-18.fasta \
   -o assembled_transcripts/$f hs2_sam/$f.sam.sorted
done < ~/projects/regen_symm/samples.txt




  ###########################################################
  # RUN CUFFMERGE TO ASSEMBLE AN EXPERIMENTAL TRANSCRIPTOME #
  ###########################################################


#!/bin/bash -l

#SBATCH -D ~/projects/regen_symm/
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#enter_email_here 
#SBATCH -o ~/projects/regen_symm/slurm-logs/cuffmerge-o%j.txt
#SBATCH -e ~/projects/regen_symm/slurm-logs/cuffmerge-e%j.txt
#SBATCH -J cuffmerge
#SBATCH -t 50:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

#Make the text file containing the list of assembly paths that inputs to cuffmerge
cd ~/projects/regen_symm
while read f; do
 echo "~/projects/regen_symm/assembled_transcripts/$f/transcripts.gtf" >> assembled_transcripts/assemblies_path.txt
done < ~/projects/regen_symm/samples.txt

cd assembled_transcripts

cuffmerge -g ~/data/Aurelia.Genome_Annotation_v1.2_11-28-18.gff3 \
  --ref-sequence ~/data/Aurelia.Genome_v1.2_11-27-18.fasta \
  assemblies_path.txt
  
  
  
  
  ################################################################
  # RUN CUFFQUANT TO QUANTIFY TRANSCRIPT ABUNDANCES PER DATA SET #
  ################################################################


#!/bin/bash -l

#SBATCH -D ~/projects/regen_symm/
#SBATCH --mail-type=ALL
#SBATCH --mail-user=#enter_email_here 
#SBATCH -o ~/projects/regen_symm/slurm-logs/cuffquant-o%j.txt
#SBATCH -e ~/projects/regen_symm/slurm-logs/cuffquant-e%j.txt
#SBATCH -J cuffquant
#SBATCH -t 15:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

cd ~/projects/regen_symm

# Re-sort sam files
mkdir hs2_sam/samtools_sort

#cuffquant requires a re-sorting of the sam files, saved in hs2_sam/samtools_sort/
while read f; do
  samtools sort $f -o hs2_sam/samtools_sort/"$f"_sorted.sam
done < ~/projects/regen_symm/samples.txt

# Save a note about new files (usage/generation)
echo -e "The files in regen_symm/hs2_sam/samtools_sort/ have been sorted from the hisat2 generated
sam files using samtools sort for use with cuffquant (ID_sorted.sam), which gave an 
unsorted error when using the files sorted using cufflinks sort for input to 
cufflinks (ID.sorted.sam)\n\n\
Error sample:
  Error: this SAM file doesn't appear to be correctly sorted!
        current hit is at Seg75:952, last one was at Seg1298:98677
Cufflinks requires that if your file has SQ records in
the SAM header that they appear in the same order as the chromosomes names 
in the alignments.
If there are no SQ records in the header, or if the header is missing,
the alignments must be sorted lexicographically by chromsome
name and by position.
" | cat > hs2_sam/samtools_sort/content_notes.txt

# Run cuffquant on newly sorted sam files
mkdir cuffq_abundances

while read f; do
 cuffquant -u -o cuffq_abundances/"$f" assembled_transcripts/merged_asm/merged.gtf hs2_sam/samtools_sort/"$f"_sorted.sam
done < ~/projects/regen_symm/samples.txt


while read i; do
 cp ~/projects/regen_symm/cuffq_abundances/"$i"/abundances.cxb ~/projects/regen_symm/cuffq_abundances/"$i".cxb
done < ~/projects/regen_symm/samples.txt

# compress to export
tar -zcvf cuffq_abundances.tar.gz ~/projects/regen_symm/cuffq_abundances/*.cxb


  ##########################################################################
  # RUN CUFFNORM TO NORMALIZE TRANSCRIPT EXPRESSION VALUES ACROSS DATA SET #
  ##########################################################################
  
#!/bin/bash -l

#SBATCH -D ~/projects/regen_symm/
#SBATCH --mail-type=ALL
#SBATCH -o ~/projects/regen_symm/slurm-logs/cuffnorm-o%j.txt
#SBATCH -e ~/projects/regen_symm/slurm-logs/cuffnorm-e%j.txt
#SBATCH --mail-user=#enter_email_here 
#SBATCH -J cuffnorm
#SBATCH -t 50:00:00

PATH=$PATH:~/tools/cufflinks-2.2.1.Linux_x86_64

cd ~/projects/regen_symm/cuffq_abundances

# cuffnorm inputs the complete transcriptome (merged.gtf), comma separated labels specific
# to the conditions(-L), and in the same order comma separated replicates of each condition:
cuffnorm --no-update-check -o ./norm \ 
../assembled_transcripts/merged_asm/merged.gtf \
-L WT_T0,SYM_24,REG_24,SYM_36,REG_36,WT_48,SYM_48,REG_48 \  
21020_1_WT_T0_1.cxb,21021_2_WT_T0_2.cxb 21022_3_SYM_24_1.cxb,21023_4_SYM_24_2.cxb,21024_5_SYM_24_3.cxb \
21025_6_REG_24_1.cxb,21026_7_REG_24_2.cxb,21027_8_REG_24_3.cxb 21029_9_SYM_36_1.cxb,21030_10_SYM_36_2.cxb,21031_11_SYM_36_3.cxb \
21032_12_REG_36_1.cxb,21033_13_REG_36_2.cxb,21034_14_REG_36_3.cxb 21035_15_WT_48_1.cxb,21036_16_WT_48_2.cxb \
21038_17_SYM_48_1.cxb,21039_18_SYM_48_2.cxb,21040_19_SYM_48_3.cxb 21041_20_REG_48_1.cxb,21042_21_REG_48_2.cxb,21043_22_REG_48_3.cxb

# compress to export
tar -zcvf norm_abundances.tar.gz ~/projects/regen_symm/cuffq_abundances/norm/*




cd ~/projects/misc/cuffq_abundances

cuffnorm -o ./norm_exc_3mL_27hr \
–output-format cuffdiff
../assembled_transcripts/merged_asm/merged.gtf \
-L 3mL_51hr,NoFood_51hr,NoFood_27hr,Reg,WT \
19198.cxb,19199.cxb 19192.cxb,19193.cxb 19188.cxb,19189.cxb \
17897.cxb,17898.cxb 17895.cxb,17896.cxb


