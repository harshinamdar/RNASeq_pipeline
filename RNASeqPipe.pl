#!/usr/bin/perl
use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

my $index;
my $annotation_file;
my $genome_fasta;
my $scripts_out;
my $project_dir;
my $is_stranded;
my $help;
my $man;
GetOptions('index=s'            => \$index,				
           'annotation_file=s'  => \$annotation_file,			
  	   'genome_fasta=s' 	=> \$genome_fasta,			
    	   'scripts_out=s'   	=> \$scripts_out,			
	   'is_stranded=s'	=> \$is_stranded,			 
	   'help!' 		=> \$help,
	   'man!' 		=> \$man,
 	  'project_dir=s'     => \$project_dir );
	  # 'project_dir=s' 	=> \$project_dir ) or pod2usage(-verbose => 0);
          # 'project_dir=s'      => \$project_dir ) or pod2usage("Try '$0 --help' for more information.") unless $ARGV;
	   pod2usage(-verbose => 1)  if ($help);
	   pod2usage(-verbose => 2)  if ($man);
unless ( defined( $project_dir ) && defined( $index ) )
    {
        Pod::Usage::pod2usage( -exitstatus => 1);
    }
my $ref_genome_index_base = $index;
my $gtf = $annotation_file;
my $parent = $project_dir;
my $strand_info = $is_stranded;
my $fullpath = $project_dir;
my $base_dir = basename($fullpath);
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
		
		my $basedir = basename($path);
                my $jobname  = "job".'_'.$basedir;
                open($outfile,">$jobname");

		### Now print list of commands to be executed for each sample_folder into output file ==> scripts_out ####
		###------------------------------------------------------------------------------------------------------####
	 	print $outfile "cd \$PBS_O_WORKDIR\n\n"; 
         	print $outfile "module add tophat/2.0.11\n\n"; 
         	print $outfile "module add bowtie2/2.1.0\n\n"; 
	        print $outfile "module add samtools/0.1.18\n\n";
		print $outfile "module add python/2.7.3\n\n";
         	print $outfile "module add htseq/0.5.4p3\n\n";
		print $outfile "module add cufflinks/2.2.1\n\n";
                print $outfile "module add rseqc/2.3.7\n\n";
	        print $outfile "tophat --library-type fr-firststrand -r 100 -p 8 -G $gtf -o $path/tophat_out_$sub_folders $ref_genome_index_base $R1 $R2\n\n";

	 	### Can submit other jobs immeditely after tophat alignment : for ex; submitting a simeltaneous cufflinks job as shown below :
		#print $outfile "echo cufflinks -p 8 -N -G $gtf  -M \$DM3_MASK_GTF -b $ref_genome_fasta -o $path/tophat_out_$sub_folders/cufflinks_$sub_folders -u --compatible-hits-norm $path/tophat_out_$sub_folders/accepted_hits.bam| qsub -N cuff_$sub_folders -V -cwd \n\n";

	 	print $outfile "samtools sort -n $path/tophat_out_$sub_folders/accepted_hits.bam $path/tophat_out_$sub_folders/sorted_bam_$sub_folders\n\n";
		## -n : sort w.r.t. to name as required in htseq;  -o : sort w.r.t. to coordinates
		print $outfile "samtools view -h $path/tophat_out_$sub_folders/sorted_bam_$sub_folders.bam >$path/tophat_out_$sub_folders/sorted_bam_$sub_folders.sam\n\n";
          	print $outfile "htseq-count -s $strand_info --order=name $path/tophat_out_$sub_folders/sorted_bam_$sub_folders.sam $gtf >$path/tophat_out_$sub_folders/$sub_folders.gene_count\n\n";
	 	print $outfile "sed -i -e '1igene $sub_folders\' $path/tophat_out_$sub_folders/$sub_folders.gene_count\n\n";
	 	print $outfile "awk -v OFS=\"\\t\" '\$1=\$1' $path/tophat_out_$sub_folders/$sub_folders.gene_count > $path/tophat_out_$sub_folders/$sub_folders.gene_counts\n\n";
	 	print $outfile "mv $path/tophat_out_$sub_folders/$sub_folders.gene_counts $parent/../DESeq_$base_dir\n\n";
		print $outfile " read_NVC.py -i $path/tophat_out_$sub_folders/accepted_hits.bam -o $path/tophat_out_$sub_folders/RseqC_$sub_folders -x \n\n";
		close($outfile);
	       `qsub -l select=1:ncpus=1:mem=12GB -l walltime=70:00:00 $jobname`;
}
__END__

=head1 NAME

RNASeqPipe.pl 

=head1 SYNOPSIS

perl RNASeqPipe.pl [options]

  Options:

   --index
   --annotation_file     
   --genome_fasta     
   --project_dir
   --scripts_out
   --is_stranded

   Example:
	./RNASeqPipe.pl --index /home/RNASeq/genome/Sequence/Bowtie2Index/genome --annotation_file /home/RNASeq/genome/annotation/genes.gtf --genome_fasta /home/RNASeq/genome/Sequence/WholeGenomeFasta/genome.fa --project_dir Project_TEST --scripts_out scripts.out --is_stranded yes

=head1 DESCRIPTION

##
Begining from tophat till htseq-count, this program creates executable commands for all samples in a Project_Dir.
It also launches these jobs onto cluster in parallel.
##

=head1 OPTIONS

=over 4

=item B<--index>

bowtie2 index prefix of the genome

=item B<--annotation_file>

.gtf file; Provide Reference annotation preferably from Ensembl

=item B<--genome_fasta>

Provide reference genome multi-fasta file

=item B<--project_dir>

Input Project directory which contains sample_folders with reads

=item B<--scripts_out>

Provide output file name; executable commands for all Sample folders will be written in this single file

=item B<--is_stranded>

yes or no or reverse as required by htseq-count

=back

=head1 COPYRIGHT

This program is free software and can redistributed and/or modified.

=cut
