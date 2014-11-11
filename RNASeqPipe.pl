#!/usr/bin/perl
use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use File::Basename;
my $index;
my $annotation_file;
my $genome_fasta;
my $scripts_out;
my $project_folder;
my $stranded;
GetOptions('index=s'            => \$index,				#bowtie index base
           'annotation_file=s'  => \$annotation_file,			#gtf file
    	   'genome_fasta=s' 	=> \$genome_fasta,			#ref. genome in multi-fasta format
    	   'scripts_out=s'   	=> \$scripts_out,			#Script_out file containing executable commands for all Sample folders
	   'stranded=s'		=> \$stranded,				#strand info required for htseq 
    	   'project_folder=s' 	=> \$project_folder );			#Main Project_Dir containing Sample folders with fastq files; R1 and R2
(!$index or !$annotation_file or !$genome_fasta or !$scripts_out or !$project_folder ) and die ("Usage : job_script.pl --index <bowtie_index> --annotation_file <gtf_file> --genome_fasta <reference_genome_fasta_file> --stranded <yes/no/reverse> --project_folder <project_directory_containing_sample_folders_with_fastq_files> --scripts_out <output_file>\n");


my $ref_genome_index_base = $index;
my $gtf = $annotation_file;
my $parent = $project_folder;
my $strand_info = $stranded;
my $fullpath = $project_folder;
my $base_dir = basename($fullpath);
#print "file = $base_dir\n";
my $DESeq_anl_fol = `mkdir $parent/../DESeq_$base_dir`; 
my $par_dir;
my $sub_dir;
my $R1;
my $R2;
my $full_path;
my $outfile;
my $ref_genome_fasta = $genome_fasta;
my $count_folder;
##################################################################################################################
################## Read recursively through Project_Dir --> Sample_folder ---> fastq_files #######################
##################-------------------------------------------------------------------------#######################

opendir($par_dir,$parent);
open($outfile,">$scripts_out");

while (my $sub_folders = readdir($par_dir)) {
        next if ($sub_folders =~ /^..?$/);  # skip '.' and '..'
        my $path = $parent . '/' . $sub_folders;
        next unless (-d $path); # skip anything that isn't a directory
        opendir ($sub_dir, $path);
        	while (my $file = readdir($sub_dir)) {
		next unless $file=~/\.f*q*/i; ## reads fastq file as < .fastq || .fq || .fastq.gz >
           	$full_path = $path . '/' . $file;
          		if ($full_path=~/\_R1/) {
                	$R1=$full_path;
        }	
          		if ($full_path=~/\_R2/) {
                	$R2=$full_path;
    	}

    }
         	close($sub_dir);
		$count_folder++;  ## counting no. of sample_folders
       		### Now print list of commands to be executed for each sample_folder into output file ==> scripts_out ####
		###------------------------------------------------------------------------------------------------------####
	 	print $outfile "cd \$PBS_O_WORKDIR\n\n"; 
         	print $outfile "module add tophat/2.0.8\n\n"; 
         	print $outfile "module add bowtie2/2.1.0\n\n"; 
	        print $outfile "module add samtools/0.1.18\n\n";
		print $outfile "module add python/2.7.3\n\n";
         	print $outfile "module add htseq/0.5.4p3\n\n";
	        print $outfile "tophat -p 8 -G $gtf -o $path/tophat_out_$sub_folders $ref_genome_index_base $R1 $R2\n\n";

	 	### Can submit other jobs immeditely after tophat alignment : for ex; submitting a simeltaneous cufflinks job as shown below :
		#print $outfile "echo cufflinks -p 8 -N -G $gtf  -M \$DM3_MASK_GTF -b $ref_genome_fasta -o $path/tophat_out_$sub_folders/cufflinks_$sub_folders -u --compatible-hits-norm $path/tophat_out_$sub_folders/accepted_hits.bam| qsub -N cuff_$sub_folders -V -cwd \n\n";

	 	print $outfile "samtools sort -n $path/tophat_out_$sub_folders/accepted_hits.bam $path/tophat_out_$sub_folders/sorted_bam_$sub_folders\n\n";
		## -n : sort w.r.t. to name as required in htseq;  -o : sort w.r.t. to coordinates
		print $outfile "samtools view -h $path/tophat_out_$sub_folders/sorted_bam_$sub_folders.bam >$path/tophat_out_$sub_folders/sorted_bam_$sub_folders.sam\n\n";
          	print $outfile "htseq-count -s $strand_info --order=name $path/tophat_out_$sub_folders/sorted_bam_$sub_folders.sam $gtf >$path/tophat_out_$sub_folders/$sub_folders.gene_count\n\n";
	 	print $outfile "sed -i -e '1igene $sub_folders\' $path/tophat_out_$sub_folders/$sub_folders.gene_count\n\n";
	 	print $outfile "awk -v OFS=\"\\t\" '\$1=\$1' $path/tophat_out_$sub_folders/$sub_folders.gene_count > $path/tophat_out_$sub_folders/$sub_folders.gene_counts\n\n";
	 	print $outfile "mv $path/tophat_out_$sub_folders/$sub_folders.gene_counts $parent/../DESeq_$base_dir\n\n";
}  
		close($outfile);

my $count_lines;	
open (my $infile,"<$scripts_out");
while (<$infile>) {
	chomp $_;
	$count_lines++; ##counting no. of lines in scripts_out
}
my $split_count = $count_lines/$count_folder; ## 
close($par_dir);

#######################################################################################################################################################
################ Split the scripts_out file to create separate submission file for each sample_folder and then submit all of them simeltaneously ######
################---------------------------------------------------------------------------------------------------------------------------------######

    my$split_command=`split -l $split_count  -d $scripts_out NEW_job`; ## split scripts_out to create separate submission file for each sample/replicate
    my@split_array = glob "NEW_job*";
   #print "$split_array[0]\n";
   	 my $pattern = $base_dir;
	 foreach my $split_file(@split_array) {
	 open (my $splitfile,"<$split_file");
	 while (my$splitline = <$splitfile>) {
		if ($splitline =~/$pattern\/\S+\/$/){
			print "...........................$1\n"; 
			}
			last;
	 }
	
   `qsub -l select=1:ncpus=1:mem=12GB -l walltime=70:00:00  $split_file`;
   }

