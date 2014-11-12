#!/usr/bin/perl
use warnings;
use strict;

# This script will print the transcript with highest fpkm value; parses file "isoforms.fpkm_tracking" produced by running cufflinks
=begin
> head isoforms.fpkm_tracking

tracking_id     class_code      nearest_ref_id  gene_id gene_short_name tss_id  locus   length  coverage        FPKM    FPKM_conf_lo    FPKM_conf_hi    FPKM_status
ENSMUST00000088658      -       -       ENSMUSG00000025912      Mybl1   TSS26267        1:9667414-9700209       4980    11.7529 3.33412 2.90674 3.76149 OK
ENSMUST00000115468      -       -       ENSMUSG00000025912      Mybl1   TSS54786        1:9669536-9700024       2493    2.40412 0.682011        0.391552        0.972471        OK
ENSMUST00000160022      -       -       ENSMUSG00000025912      Mybl1   TSS66543        1:9676131-9678610       543     0.568679        0.161326        0       0.699656        OK
ENSMUST00000160451      -       -       ENSMUSG00000025912      Mybl1   TSS49637        1:9678271-9683452       695     3.47795 0.986642        0.15403 1.81925 OK

=cut
#usage: perl extract_abundant_transcript.pl isoforms.fpkm_tracking

my $infile;
open($infile, "<$ARGV[0]") or die "cant find file\n" ; 
my %multi_hash;
my $header = <$infile>; ## capture header to avoid being read in while loop;
while (my $line = <$infile>) {
    chomp $line;
    my @arr = split(/\t/, $line);
    my $gene_id  = $arr[3];
    my $trans_id = $arr[0];
    my $fpkm     = $arr[9];
    $multi_hash{$gene_id}{$trans_id} = $fpkm;
}

for my $gen_id (sort keys %multi_hash) {
    my %h = %{ $multi_hash{$gen_id} };
    my $key = (reverse sort { $h{$a} <=> $h{$b} } keys %h)[0]; ## [0] indicates first value in array; which is the highest fpkm value
    print "$gen_id\t$key\t$multi_hash{$gen_id}{$key}\n";
}
