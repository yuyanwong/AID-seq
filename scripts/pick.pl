#!/usr/bin/perl

my $sample = shift; # G1800
my $sgRNA_file = $sample."_sgRNA.txt";

my %hash_window = ();
my %hash_num = ();
my %hash_all_sgRNA = ();

open SGRNA, "$sgRNA_file";
while(<SGRNA>){
	my @F = split(/\t/,$_);
	my $sgRNA = "$sample-$F[0]";
	$hash_all_sgRNA{$sgRNA} = 1;
	my $file = "${sample}/${sgRNA}_identified_matched.txt";
	open FILE, "$file";
	while(my $line = <FILE>){
		next if($line =~ /^#/);
		
		my @G = split(/\t/,$line);
		my $window = $G[10];
		
		if(not exists($hash_num{$window})){
			$hash_num{$window} = $G[13];
			$hash_window{$window} = $sgRNA;
		}else{
			if($G[13] < $hash_num{$window}){
				$hash_num{$window} = $G[13];
				$hash_window{$window} = $sgRNA;
			}elsif($G[13] == $hash_num{$window}){
				$hash_window{$window} .= ";".$sgRNA;
			}
		}
	}
	close FILE;
}
close SGRNA;

my %hash_sgRNA = ();

foreach my $window (keys %hash_window){
	my @sgRNAs = split(/;/,$hash_window{$window});
	foreach my $sgRNA(@sgRNAs){
		if(!$hash_sgRNA{$sgRNA}){
			$hash_sgRNA{$sgRNA} = $window;
		}else{
			$hash_sgRNA{$sgRNA} = $hash_sgRNA{$sgRNA}.";".$window;
		}
	}
}

if(!(-d "${sample}/pick/") ){
        `mkdir -p ${sample}/pick/`;
}

foreach my $sgRNA(keys %hash_all_sgRNA){
	my $file = "${sample}/${sgRNA}_identified_matched.txt";
	my $out = "${sample}/pick/${sgRNA}_identified_matched.txt";
	open FILE, "$file";
	open OUT, ">$out";
	while(my $line = <FILE>){
		if($line =~ /^#/){
			print OUT "$line";
			next;
		}else{
			if($hash_sgRNA{$sgRNA}){
				my @arr = split(/;/,$hash_sgRNA{$sgRNA});
				my @F = split(/\t/,$line);
				my $window = $F[10];
				if(grep{$_ eq $window} @arr){
					print OUT "$line";
				}
			}
		}
	}
	close FILE;
	close OUT
}
