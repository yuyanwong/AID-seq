#!/bin/bash

check_opts(){
	if [ -z "$config" ]; then
		echo "use -c to specify the config file"
		exit 1
	fi
}

if [ "$#" -lt 1 ]; then
cat <<HELP
use option -h to get more information. 
HELP
exit 0
fi

helpdoc(){
	cat <<HELP
NAME
	aidseq -- analysis pipeline to find the offtarget using SITE-seq methods
SYNOPSIS
	sh aidseq.sh [OPTION]... [NEMBER]...
DESCRIPTION
	aidseq.sh is used to analysis the data generated from AID-seq methods.
OPTIONS
	-h help text
	-c config.txt
HELP
}

while getopts "c:h" opt
do
	case $opt in
	c)
	config=$OPTARG
	;;
	h)
	helpdoc
	exit 0
	;;
	esac
done

while read line;do
	eval "$line"
done < $config
echo "bowtie2_index: "$bowtie2_index
echo "genome: "$genome
echo "data_dir: "$data_dir
echo "sample: "$sample
echo "layout: "$layout
echo "bam_prefix: "$bam_prefix
echo "output_prefix: "$output_prefix
echo "suffix: "$suffix
echo "need: "$need
echo "flank_from_feature: "$flank_from_feature
echo "required_cliff_ratio: "$required_cliff_ratio
echo "min_dup: "$min_dup
echo "result_dir: "$result_dir
echo "threads: "$threads
echo "target_sequence: "$target_sequence

script_dir=$(dirname $(readlink -f "$0"))
echo "script_dir: "$script_dir

if [ ! -d $result_dir/$bam_prefix ]; then
	mkdir -p $result_dir/$bam_prefix
	cd $result_dir/$bam_prefix
else
	cd $result_dir/$bam_prefix
fi

if [ ! -s ${bam_prefix}.bam ]; then
	if [ ${layout} == 'SE' ]; then
		bowtie2 -x $bowtie2_index -p $threads -U $data_dir/${sample}.${suffix} -S ${bam_prefix}.sam
	elif [ ${layout} == 'PE' ]; then
		bowtie2 -x $bowtie2_index -p $threads -1 $data_dir/${sample}_R1.${suffix} -2 $data_dir/${sample}_R2.${suffix} -S ${bam_prefix}.sam
	fi
	samtools sort -@ $threads ${bam_prefix}.sam -o ${bam_prefix}.bam
	samtools index ${bam_prefix}.bam
fi

if [ -s ${output_prefix}_identified_matched.txt ]; then
	python $script_dir/AID-Seq_core_feature_calling_functions_new.py --layout $layout --target_sequence $target_sequence \
		--matched_file ${output_prefix}_identified_matched.txt --bam_file ${bam_prefix}.bam \
		--genome $genome --out_name ${output_prefix} --flank_from_feature $flank_from_feature --required_cliff_ratio $required_cliff_ratio --min_dup $min_dup --mismatch_threshold 6
elif [ -s ${bam_prefix}_candidate.fa ]; then
	python $script_dir/AID-Seq_core_feature_calling_functions_new.py --layout $layout --target_sequence $target_sequence \
		--candidate_fasta_file ${bam_prefix}_candidate.fa --bam_file ${bam_prefix}.bam \
		--genome $genome --out_name ${output_prefix} --flank_from_feature $flank_from_feature --required_cliff_ratio $required_cliff_ratio  --min_dup $min_dup --mismatch_threshold 6
elif [ -s ${bam_prefix}.peak ]; then
	python $script_dir/AID-Seq_core_feature_calling_functions_new.py --layout $layout --target_sequence $target_sequence \
		--initial_peaks_file ${bam_prefix}.peak --bam_file ${bam_prefix}.bam --genome $genome \
		--out_name ${output_prefix} --flank_from_feature $flank_from_feature --required_cliff_ratio $required_cliff_ratio --min_dup $min_dup --mismatch_threshold 6
else
	python $script_dir/AID-Seq_core_feature_calling_functions_new.py --layout $layout --target_sequence $target_sequence \
		--bam_file ${bam_prefix}.bam --genome $genome --out_name ${output_prefix} --flank_from_feature $flank_from_feature --required_cliff_ratio $required_cliff_ratio --min_dup $min_dup --mismatch_threshold 6
fi
