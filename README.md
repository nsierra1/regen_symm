
Time course differential gene expression between resymmetrizing and regenerating Aurelia polyps

PIPELINE FOR PROCESSING PAIRED END RNA-seq DATA USING HISAT2(MAPPING) 
AND CUFFLINKS(ASSEMBLY AND QUANTIFICATION)
Noemie Sierra, July 2019	

The following script was used to process the data generated for 

To run the following analysis you will need :
	# Hisat2 Aligner 
		# download here: https://ccb.jhu.edu/software/hisat2/manual.shtml
	# Cufflinks Transcriptome Assembly/DE 
		#download here:http://cole-trapnell-lab.github.io/cufflinks/install/
	
	
	# HOME (~) FOLDER CONTAINS: 
	  # PROJECT FOLDER
		# fastq_files - containing PE .fastq.gz data - label_R1 & label_R2
		# slurm-logs - for stdout/stderr output
	  # REFERENCE DATA FOLDER
		# GENOMES and other reference info; all indices and 
	 	# other reference file manipulations are saved here
		
	# Each of the following was run as a separate script on the Peloton Cluster (CSE, UC Davis, 2019)
	# using the following slurm queue submission:
		# sbatch -p (med/high) -n (cores) script
