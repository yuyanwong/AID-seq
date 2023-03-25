#!/usr/bin/perl -w
use strict;

my $in = shift;
my $type = shift;
my $type_dir="";

if($type eq "all"){open OUT, ">$in/$in.all.stat.xls" || die $!;$type_dir="";}else{open OUT, ">$in/$type/$in.$type.stat.xls" || die $!;$type_dir=$type;}
print OUT "sample\tsgRNA\ton\toff\ton_reads\toff_reads\toff_reads\toff_reads\toff_reads\toff_reads\toff_reads\toff_reads\toff_reads\toff_reads\toff_reads\n";

open IN, $in."_sgRNA.txt" || die "Cannot open file: $in\n";

my %on = ();
my %off = ();
my @samples = ();

while(my $line1 = <IN>){
	chomp($line1);
	my @F = split(/\t/,$line1);
	my $sample = $F[0];
	my $sgRNA = $F[1];
	push @samples, $sample;
	my $file = "$in/$type_dir/$in"."-$sample"."_identified_matched.txt";
	open FILE, "$file" || die "Cannot open file: $file\n";
	my $on = 0;
	my $off = 0;
	my $reads = 0;
	my %off_reads = ();
	while(my $line = <FILE>){
		next if($line =~ /^#/);
		chomp($line);
		my @F = split(/\t/,$line);
		if($F[13] eq "" or $F[13] != 0){
			$off++;
			$off_reads{$off} = $F[6];
			# if(!$off{$sample}){
				# $off{$sample} = 1;
			# }else{
				# $off{$sample} += 1;
			# }
		}elsif($F[13] == 0){
			$on++;
			if($reads <= $F[6]){
				$reads = $F[6];
			}
			# if(!$on{$sample}){
				# $on{$sample} = 1;
			# }else{
				# $on{$sample} += 1;
			# }
		}
	}
	close FILE;
	
	print OUT "$sample\t$sgRNA\t$on\t$off\t$reads";
	foreach my $i(sort{$a<=>$b} keys %off_reads){
		if($i>10){
			last;
		}
		print OUT "\t$off_reads{$i}";
	}
	print OUT "\n";
}
close IN;
close OUT;
