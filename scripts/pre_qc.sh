#!/bin/bash

datadir=$1
sample=$2
script=$3

#mkdir $datadir/clean

# cutadapt -a AGATCGGAAGAG\
	# -A AGATCGGAAGAG \
	# -g CTCTTCCGATCT \
	# -G CTCTTCCGATCT \
	# -o $datadir/clean/${sample}_R1.fastq.gz \
	# -p $datadir/clean/${sample}_R2.fastq.gz \
	# $datadir/${sample}_R1.fastq.gz \
	# $datadir/${sample}_R2.fastq.gz

cutadapt \
	-a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-O 5 \
	-j 4 \
 	--info-file=$datadir/${sample}_R1.adapter.txt \
	$datadir/${sample}_R1.fastq.gz \
	1>/dev/null 2>/dev/null \
	&&  gzip -c $datadir/${sample}_R1.adapter.txt \
	> $datadir/${sample}_R1.adapter.txt.gz \
	&& rm $datadir/${sample}_R1.adapter.txt &

cutadapt \
	-a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
	-O 5 \
	-j 4 \
	--info-file=$datadir/${sample}_R2.adapter.txt \
	$datadir/${sample}_R2.fastq.gz \
	1>/dev/null 2>/dev/null \
	&&  gzip -c $datadir/${sample}_R2.adapter.txt \
	> $datadir/${sample}_R2.adapter.txt.gz \
	&& rm $datadir/${sample}_R2.adapter.txt &
wait

$script/../software/fqtools sfpe -q 19 -Q 0.15 $datadir/${sample}_R1.fastq.gz $datadir/${sample}_R2.fastq.gz $datadir/${sample}_R1.adapter.txt.gz $datadir/${sample}_R2.adapter.txt.gz $datadir/clean/${sample}_R1.fastq.gz  $datadir/clean/${sample}_R2.fastq.gz 
