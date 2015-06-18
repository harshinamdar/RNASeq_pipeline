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
executable statements. The `.gene_counts` for each sample will be available in folder `DESeq_Project_dir`. 

This script currently caters to single Project_dir. Future update will see processing of multiple Project_dir.

For Quick tutorial on DESeq [Visit Here](http://harshinamdar.wordpress.com/2014/11/11/quick-tutorial-on-deseq2/) or execute the script `run_deseq2.R` once all the `.gene_counts` files are generated for each sample. 
      
The script `run_deseq2.R` requires two arguments 

      Rscript run_deseq2.R /path/to/DESeq_Project_dir metaData.csv
     
      #example metaData.csv
      filename,sample,type
      Sample_Mated1.gene_counts,mated1,mated
      Sample_Mated2.gene_counts,mated2,mated
      Sample_Mated3.gene_counts,mated3,mated
      Sample_Virgin1.gene_counts,virgin1,virgin
      Sample_Virgin2.gene_counts,virgin2,virgin
      Sample_Virgin3.gene_counts,virgin3,virgin
      
      

    
