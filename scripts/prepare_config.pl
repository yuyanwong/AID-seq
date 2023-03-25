#!/usr/bin/perl -w
use strict;

my $in = shift;
my $clean = shift;
my $result = shift;
my $sample = $in."_sgRNA.txt";
my $genome = shift; ##genome name
my $genome_data= shift; ##genome data dir
my $script=shift;##scripts dir
open IN, "$sample" || die "Cannot open file: $sample\n";

my %hash = ();
while(<IN>){
	chomp;
	my @F = split(/\t/, $_);
	$hash{$F[0]}{$F[1]} = $F[2];
}
close IN;

my $submit = "submit.sh";
my $config = $in."_config.txt";
open SUBMIT, ">$in.submit.sh";
print SUBMIT "#!/bin/bash\n\n";
print SUBMIT "script=$script\n";
print SUBMIT "result=$result\n\n";
print SUBMIT "perl \$script/prepare_new.pl \$result/$config $in"."-run.sh\n";

my $run = $in."_run.sh";
open RUN, ">$run";
print RUN "#!/bin/bash\n\n";
print RUN "nohup sh $in-run.sh &\n";

open CONFIG, ">$config" || die $!;
print CONFIG "bowtie2_index=$genome_data/$genome\n";
print CONFIG "genome=$genome_data/$genome.fa\n";
print CONFIG "data_dir=$clean\n";
print CONFIG "suffix=fastq.gz\n";
print CONFIG "layout=PE\n";
print CONFIG "need=single\n";
print CONFIG "flank_from_feature=20\n";
print CONFIG "required_cliff_ratio=0.9\n";
print CONFIG "min_dup=5\n";
print CONFIG "result_dir=$result\n";
print CONFIG "threads=16\n";
print CONFIG "sample\ttarget_sequence\tbam_predix\toutput_prefix\n";

foreach my $type1(keys %hash){
	foreach my $sgRNA(keys %{$hash{$type1}}){
		print CONFIG "$in\t$hash{$type1}{$sgRNA}\t$in\t$in-$type1\n";
	}
}
# print CONFIG "$in\tGGGAAAGACCCAGCATCCGTNGG\t$in\t$in-S1\n";

close CONFIG;
close RUN;
close SUBMIT;

