  RNASeq_pipeline
===============
This script creates a series of commands starting from alignment till gene-counts for each of the samples in Project directory 
as shown below.

              +-- Project_dir
              ¦   +-- Sample_A
              ¦   ¦   +-- A_L001_R1_001.fastq
              ¦   ¦   +-- A_L001_R2_001.fastq
              ¦   +-- Sample_B
              ¦   ¦   +-- B_L001_R1_001.fastq
              ¦   ¦   +-- B_L001_R2_001.fastq
              ¦   +-- Sample_C
              ¦   ¦   +-- C_L001_R1_001.fq
              ¦   ¦   +-- C_L001_R2_001.fq
              ¦   +-- Sample_D
              ¦       +-- D_L001_R1_001.fastq
              ¦       +-- D_L001_R2_001.fastq.gz

The jobs for each of the samples will be launched in parallel onto cluster. The script can be customised to include additional 
executable staements. 
