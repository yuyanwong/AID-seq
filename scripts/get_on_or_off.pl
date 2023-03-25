#!/usr/bin/perl -w
use strict;

my $result = shift;
my $sample_file = shift;

open IN1, "$sample_file";

open ON, ">$result/on.txt";
open OFF, ">$result/off.txt";
#print ON "Sample\tChromosome\tStart\tEnd\tTotal Reads\tNon-cliff Reads\tCliff Reads\tForward Reads\tReverse Reads\tStrand\tSequence\n";
# print OFF "Sample\tChromosome\tStart\tEnd\tTotal Reads\tNon-cliff Reads\tCliff Reads\tForward Reads\tReverse Reads\tStrand\tSequence\tSubstitutions\tInsertions\tDeletions\n";
#print OFF "Sample\tChromosome\tStart\tEnd\tTotal Reads\tNon-cliff Reads\tCliff Reads\tForward Reads\tReverse Reads\tStrand\tSequence\tSubstitutions\n";
my $NF='';
while(<IN1>){
	chomp;
	my $sample = $_;
        my @prefix = split("-",$_);
        print "$prefix[0]\n";
	open IN, "$result/$sample/$sample"."_identified_matched.txt";

	while(<IN>){
		next if /^#/;
		my @F = split(/\t/,$_);
		my $len=@F;
		$NF=qx(awk '{if(NR==1)print NF}' $result/$sample/$sample*_identified_matched.txt );
		if ($len>26){
		 if($F[13] eq ""){
			 print OFF "$sample\t$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[17]\t$F[20]\t$F[21]\t$F[21]\n";
		 }
                 else{
			if($F[13] == 0){
				print ON "$sample\t$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[12]\n";
			}
			if($F[13] != 0){
				print OFF "$sample\t$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[8]\t$F[9]\t$F[12]\t$F[13]\n";
			}
		}
                    # else{
		#	if($F[11] == 0){
		#		print ON "$sample\t$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[10]\n";
		#	}
		#	if($F[11] != 0){
		#		print OFF "$sample\t$F[0]\t$F[1]\t$F[2]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\t$F[10]\t$F[11]\n";
		#	}}
			
		 }
	}
	close IN;
}
if($NF >26 ){
	qx(sed -i '1i Sample\\tChromosome\\tStart\\tEnd\\tTotal Reads\\tNon-cliff Reads\\tCliff Reads\\tForward Reads\\tReverse Reads\\tStrand\\tSequence' $result/on.txt);
	qx(sed -i '1i Sample\\tChromosome\\tStart\\tEnd\\tTotal Reads\\tNon-cliff Reads\\tCliff Reads\\tForward Reads\\tReverse Reads\\tStrand\\tSequence\\tSubstitutions' $result/off.txt);
}else{
	qx(sed -i '1i Sample\\tChromosome\\tStart\\tEnd\\tTotal Reads\\tForward Reads\\tReverse Reads\\tStrand\\tSequence' $result/on.txt);
	qx(sed -i '1i Sample\\tChromosome\\tStart\\tEnd\\tTotal Reads\\tForward Reads\\tReverse Reads\\tStrand\\tSequence\\tSubstitutions' $result/off.txt);
}
close IN1;
close ON;
close OFF;
