#!/usr/bin/perl -w
use strict;

my $list = shift;
my $datadir = shift;
my $resultdir = shift;
my $genome_version = shift;
my $genome = shift; ##hg38, hg19, hg19_hpv18 GRCm38_hpv16
my $min_dup = shift;
my $flank_from_feature = shift;
my $banned_chroms = shift;
my $nuc = shift;
my $genome_data= shift; ##genome data dir
my $script=shift;##scripts dir


## pre QC
if (!(-d $resultdir)){
	`mkdir $resultdir`;
}
if (!(-d "$datadir/clean")){
	`mkdir $datadir/clean`;
}
if (!(-f "$resultdir/sample.xls")){
	`cp $list $resultdir/sample.xls`;
}

## submit
my $submit = "$resultdir/submit.sh";
open SUBMIT, ">$submit" || die "Cannot open file: $submit\n";
print SUBMIT <<FLAG;
#!/bin/bash
script=$script
result=$resultdir

perl \$script/prepare_new.pl \$result/config.txt run.sh
FLAG
close SUBMIT;

## config file
my $config = "$resultdir/config.txt";
open CONFIG, ">$config" || die "Cannot open file: $config\n";

# if($genome =~ /^GRCm/){
	# $genome_version = $genome;
# }elsif($genome =~ /_/){
	# my @F = split(/_/, $genome);
	# $genome_version = $F[0]."_chr_".$F[1];
# }else{
	# $genome_version = $genome."_chr";
# }
if($banned_chroms eq "default"){
	print CONFIG <<FLAG;
bowtie2_index=$genome_data/$genome
genome=$genome_data/$genome.fa
data_dir=$datadir/clean
suffix=fastq.gz
layout=PE
need=single
flank_from_feature=$flank_from_feature
banned_chroms=default
required_cliff_ratio=0.9
min_dup=$min_dup
result_dir=$resultdir
Nuc=$nuc
threads=8
sample	target_sequence	bam_prefix	output_prefix
FLAG
}elsif($banned_chroms eq "NONE"){
	print CONFIG <<FLAG;
bowtie2_index=$genome_data/$genome
genome=$genome_data/$genome.fa
data_dir=$datadir/clean
suffix=fastq.gz
layout=PE
need=single
flank_from_feature=$flank_from_feature
banned_chroms=$banned_chroms
required_cliff_ratio=0.9
min_dup=$min_dup
result_dir=$resultdir
threads=8
sample	target_sequence	bam_prefix	output_prefix
FLAG
}


## pipeline file
my $pipeline = "$resultdir/pipeline.sh";
open OUT, ">$pipeline" || die $!;
print OUT "#!/bin/bash\n\n";

my %T = ();
open IN, "$list" || die "Cannot open file: $list\n";
open SAMPLE, ">$resultdir/sample.txt" || die $!;
while(<IN>){
	chomp;
	my @F = split(/\t/, $_);
	my $sample = $F[3];
	my $type = $F[4];
	my $sgRNA = $F[5];
	if($type eq "T"){
		$T{$sample}=$sgRNA;
		print CONFIG "$sample\t$sgRNA\t$sample\t$sample\n";
		print OUT "sh $script/pre_qc.sh $datadir $sample $script &\n";
		print SAMPLE "$sample\n";
	}else{
		die "has control, use generate_pipeline_control.pl \n";
	}
}
close IN;
close CONFIG;

print OUT "wait\n";
print OUT "sh $resultdir/submit.sh\n";
print OUT "sh $resultdir/run.sh\n";
print OUT "wait\n";

print OUT "perl $script/get_on_or_off.pl $resultdir $resultdir/sample.txt\n";


