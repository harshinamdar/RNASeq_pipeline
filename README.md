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
executable statements. The gene\_count table for each sample will be available in folder DESeq\_Project_dir. 

This script currently caters to single Project\_folder. Future update will see processing of multiple Project\_dir.

For Quick tutorial on DESeq [Visit Here](http://harshinamdar.wordpress.com/2014/11/11/quick-tutorial-on-deseq2/)
