#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);

my $config = shift;

open IN, "$config" || die "can not open file: $config\n";

my $bowtie2_index = "";
my $genome = "";
my $data_dir = "";
my $suffix = "";
my $layout = "";
my $result_dir = "";
my $threads = 0;
my $flag = "FALSE";
my $sample_name = "";
my $target_sequence = "";
my $bam_prefix = "";
my $output_prefix = "";
my $need = "";
my $flank_from_feature = "";
my $banned_chroms = "";
my $required_cliff_ratio = "";
my $min_dup = "";
my $nuc="";

my %hash = ();

while(<IN>){
	if($_ =~ m/^bowtie2_index=(\S+)/){
		$bowtie2_index = $1;
	}
	if($_ =~ m/^genome=(\S+)/){
		$genome = $1;
	}
	if($_ =~ m/^data_dir=(\S+)/){
		$data_dir = $1;
	}
	if($_ =~ m/^suffix=(\S+)/){
		$suffix = $1;
	}
	if($_ =~ m/^layout=(\S+)/){
		$layout = $1;
	}
	if($_ =~ m/^need=(\S+)/){
		$need = $1;
	}
	if($_ =~ m/^flank_from_feature=(\S+)/){
		$flank_from_feature = $1;
	}
	if($_ =~ m/^banned_chroms=(\S+)/){
		$banned_chroms = $1;
	}
	if($_ =~ m/^required_cliff_ratio=(\S+)/){
		$required_cliff_ratio = $1;
	}
	if($_ =~ m/^min_dup=(\S+)/){
		$min_dup = $1;
	}
	if($_ =~ m/^result_dir=(\S+)/){
		$result_dir = $1;
	}
	if($_ =~ m/^threads=(\S+)/){
		$threads = $1;
	}
	if($_ =~ m/^sample/){
		$flag = "TRUE";
		next;
	}
	if($_ =~ m/^Nuc=(\S+)/){
		$nuc = $1;
	}
	if($flag eq "TRUE" and $_ =~ m/(\S+)\t(\S+)\t(\S+)\t(\S+)/){
		$sample_name = $1;
		$target_sequence = $2;
		$bam_prefix = $3;
		$output_prefix = $4;
		if(!$hash{$sample_name}){
			$hash{$sample_name} = $output_prefix;
		}else{
			$hash{$sample_name} = $hash{$sample_name}.",".$output_prefix;
		}
		open OUT, ">$output_prefix"."_config.txt" || die $!;
		print OUT "bowtie2_index=$bowtie2_index\n";
		print OUT "genome=$genome\n";
		print OUT "data_dir=$data_dir\n";
		print OUT "sample=$sample_name\n";
		print OUT "bam_prefix=$bam_prefix\n";
		print OUT "output_prefix=$output_prefix\n";
		print OUT "suffix=$suffix\n";
		print OUT "layout=$layout\n";
		print OUT "need=$need\n";
		print OUT "flank_from_feature=$flank_from_feature\n";
		print OUT "banned_chroms=$banned_chroms\n";
		print OUT "required_cliff_ratio=$required_cliff_ratio\n";
		print OUT "min_dup=$min_dup\n";
		print OUT "result_dir=$result_dir\n";
		print OUT "Nuc=$nuc\n";
		print OUT "threads=$threads\n";
		print OUT "target_sequence=$target_sequence\n";
		close OUT;
	}
}
close IN;

open ALL_RUN, ">run.sh" || die $!;
print ALL_RUN "#!/bin/bash\n\n";
foreach my $sam(keys %hash){
	open RUN, ">$sam"."_run.sh" || die $!;

	print RUN "#!/bin/bash\n\n";
	print RUN "script=$Bin\n";
        print RUN "result_dir=$result_dir\n";
        print RUN "export PATH=\$HOME/miniconda3/bin:\$PATH\n";
	print RUN "source activate aidseq\n\n";
        my @F = split(/,/,$hash{$sam});
	foreach my $i (@F){
		if($need eq "single"){
			print RUN "sh \$script/aidseq_only_R1_new.sh -c \$result_dir/"."$i"."_config.txt >>\$result_dir/"."$i.log\n";
		}elsif($need eq "pair"){
			print RUN "sh \$script/aidseq.sh -c $i"."_config.txt >>$i.log\n";
		}
	}
	close RUN;
	print ALL_RUN "sh $sam"."_run.sh &\n";
}
print ALL_RUN "wait\n";
close ALL_RUN;
