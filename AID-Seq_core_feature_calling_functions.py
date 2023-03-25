from __future__ import print_function
from __future__ import division
import pysam
import re
import regex
from collections import Counter
import argparse
from pyfaidx import Fasta
import os
import svgwrite

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', '-': '#B3B3B3'}

def find_initial_read_pileups(bam_file,
                              banned_chroms=["alt", "random", "Un"],
                              depth_thresh=5):
    """
    Identifies any pileups in aligned sequencing reads that exceed
    a user specified threshold. Returns the genomic locus (coordinates)
    of contiguous sequence above the pileup threshold.

    Args:
        bam_file (str):
            The full path to a bam file that has been sorted and
            indexed.
        depth_tresh (int):
            The minimum number of reads required to call a pileup.
            Default: 5
        banned_chroms (list):
            A list of strings that cannot exist within the chromosome
            name for a given pileup location.
            Default: ["alt", "random", "Un", "chrM"]
    Returns:
        peak_locations (list):
            A list of genomic coordinates in the form chr1:1234-5678
            that exceed the read pileup threshold.
    """
    peak_locations = []
    bam = pysam.AlignmentFile(bam_file, 'rb')
    pileup = bam.pileup()
    in_peak = False
    # Various parameters for the peak identification.
    peak_start = -1
    peak_end = -1
    last_chrom = -1
    peak_start = -1
    peak_end = -1
    last_pos = -1
    for col in pileup:
        # Ignore undesired chromosomes.
        bad_chroms = [b for b in banned_chroms if b in col.reference_name]
        if any(bad_chroms):
            continue
        if in_peak:
            # If any of these are true, terminate the peak.
            if(int(col.nsegments) < int(depth_thresh)) or \
                    (int(col.reference_pos) != (last_pos + 1)) or \
                    (col.reference_name != last_chrom):
                in_peak = False
                peak_end = last_pos
                # Ignore 1 length peaks.
                # if peak_end == peak_start:
                    # continue
                coord = "{c}:{s}-{e}".format(c=last_chrom,
                                             s=peak_start,
                                             e=peak_end)
                peak_locations.append(coord)
        else:
            if(int(col.nsegments) >= int(depth_thresh)):
                in_peak = True
                peak_start = col.reference_pos
        last_pos = col.reference_pos
        last_chrom = col.reference_name
    sample_name = os.path.splitext(bam_file)[0]
    peak_file = '{0}.peak'.format(sample_name)
    with open(peak_file, "w+") as f:
        for peak in peak_locations:
            f.write(peak + "\n")
    
    return peak_locations


"""def call_site_seq_features(initial_peaks, bam_file,
                           pyfaidx_fasta_obj, out_name,
                           min_dup=5, min_read_length=60,
                           flank_from_feature=20,
                           required_cliff_ratio=0.9,
                           banned_chroms=["alt", "random", "Un", "chrM"]):"""
"""def call_site_seq_features(initial_peaks_file, bam_file,
                           genome, out_name, adjacent,
                           min_dup=5, min_read_length=60,
                           flank_from_feature=20,
                           required_cliff_ratio=0.9,
                           banned_chroms=["alt", "random", "Un", "chrM"]):"""
def call_site_seq_features(layout,
                           initial_peaks_file, bam_file,
                           genome,
                           min_dup=5, min_read_length=60,
                           flank_from_feature=20,
                           required_cliff_ratio=0.9,
                           banned_chroms=["alt", "random", "Un", "chrM"]):
    """
    This function takes peaks and determines whether or
    not they are site-seq features. The steps for feature
    detection are as follows:

    1) Remove peaks that lie within undesired chromosomes.
    2) For each peak, identify the count for locations where sequencing
       reads terminate (cliffs, either 5' or 3') at the same position.
       Reads that are not long enough (determined by min_read_length)
       are removed from the counts.
    3) Ensure that the number of cliff reads for that peak is high enough
       (determined by min_read_length). If the read count for the cliff is
       too low the potential feature is ignored.
    4) Ensure that the ratio of non-cliff reads to cliff reads is not
       too high. If the ratio is too high (determined by required_cliff_ratio)
       then the potential feature is ignored.
    5) All potential features that pass the above filters are reported.

    Args:
        initial_peaks (list):
            A list of peaks already found. These are strings formatted like
            'chr#:#-#'
        bam_file (str):
            The path and name of a sorted and indexed bam file.
        pyfaidx_fasta_obj (pyFaidx Fasta Object):
            A pyFaidx object created from a reference genome.
        #out_name (str):
        #    The name/path of the output fasta file. The fasta file contains
        #    information about the peak in the header as well as the sequence
        #    around
        min_dup (int):
            The number of reads that terminate at the exact same location
            (cliff_reads) needed for a potential feature to be considered.
            Default: 5 reads
        min_read_length (int):
            The minimum read length required. Reads shorter than this will
            be removed.
            Default: 60.
        flank_from_feature (int):
            Sets the sequence around the SITE-Seq feature to return in the
            fasta file output by the function.
            Default: 20.
        required_cliff_ratio (float):
            The minimum ratio of cliff reads to cliff spanning reads
            required in order for a feature to be considered.
            Default: 0.9.
        banned_chroms (list):
            A list of chromosome names that should not be considered when
            attempting to identify features. This may need to be
            changed depending on the reference used.

    Returns:
        out_name (str):
            The path to the fasta file that contains the features.
    """

    if not initial_peaks_file is None:
        with open(initial_peaks_file, 'r') as f:
            initial_peaks = f.readlines()
    else:
        initial_peaks = find_initial_read_pileups(bam_file,
                                                  banned_chroms=["alt", "random", "Un"],
                                                  depth_thresh=min_dup)

    bam = pysam.AlignmentFile(bam_file, 'rb')
    features = []

    pyfaidx_fasta_obj = Fasta(genome)

    for i, peak in enumerate(initial_peaks):
        r1_starts = []
        chrom, start, end = parse_coordinate(peak)
        # print("{0}:{1}-{2}".format(chrom, start, end))
        if any(p for p in banned_chroms if p in chrom):
            # print("banned_chroms")
            continue
        try:
            if start - 5 < 0:
                reads = bam.fetch(chrom, 0, end + 5)
            elif end + 5 > len(pyfaidx_fasta_obj[chrom]):
                reads = bam.fetch(chrom, start - 5, len(pyfaidx_fasta_obj[chrom]))
            else:
                reads = bam.fetch(chrom, start - 5, end + 5)
        except:
            # print("no reads")
            continue
        for read in reads:
            PE_layout = "PE"
            if layout is PE_layout and read.flag & 128:
                continue
            read_start = read.reference_start
            read_end = read.reference_end
            read_on_minus = read.is_reverse
            if read_start is None or read_end is None:
                continue
            aligned_length = abs(int(read_start) - int(read_end))
            if aligned_length < min_read_length:
                continue
            read_append = read_end if read_on_minus else read_start
            r1_starts.append(read_append)
        terminations = Counter(r1_starts)
        if len(terminations) == 0:
            # print("length")
            continue
        # r1_start = terminations.most_common(1)[0][0]
        # r1_count = terminations.most_common(1)[0][1]
        for r1_start, r1_count in terminations.items():
            if r1_count < min_dup:
                continue
            non_cliff_reads_at_cliff = 0
            total_reads = 0
            reverse_reads = 0
            forward_reads = 0
            if not r1_start:
                continue
            if r1_start - 1 < 0:
                reads = bam.fetch(chrom, 0, r1_start + 1)
            elif r1_start + 1 > len(pyfaidx_fasta_obj[chrom]):
                reads = bam.fetch(chrom, r1_start - 1, len(pyfaidx_fasta_obj[chrom]))
            else:
                reads = bam.fetch(chrom, r1_start - 1, r1_start + 1)
            for read in reads:
                total_reads += 1
                if read.reference_start != r1_start and \
                        read.reference_end != r1_start:
                    non_cliff_reads_at_cliff += 1
                else:
                    if read.is_reverse:
                        reverse_reads += 1
                    else:
                        forward_reads += 1
            try:
                cliff_ratio = non_cliff_reads_at_cliff / r1_count
            except ZeroDivisionError:
                cliff_ratio = 0
            # if cliff_ratio >= required_cliff_ratio:
                # print("r1_count:{0}".format(r1_count))
                # print("cliff_ratio:{0}".format(cliff_ratio))
                # continue
            fasta_dict = {}
            fasta_dict["total_reads"] = total_reads
            fasta_dict["reverse_reads"] = reverse_reads
            fasta_dict["forward_reads"] = forward_reads
            fasta_dict["chrom"] = chrom
            fasta_dict["position"] = "{s}".format(s=r1_start)
            if r1_start - int(flank_from_feature) < 0:
                fasta_dict["feature_region"] = "{c}:{s}-{e}".format(
                    c=chrom, s=0,
                    e=r1_start + int(flank_from_feature))
            elif r1_start + int(flank_from_feature) > len(pyfaidx_fasta_obj[chrom]):
                fasta_dict["feature_region"] = "{c}:{s}-{e}".format(
                    c=chrom, s=r1_start - int(flank_from_feature),
                    e=len(pyfaidx_fasta_obj[chrom]))
            else:
                fasta_dict["feature_region"] = "{c}:{s}-{e}".format(
                    c=chrom, s=r1_start - int(flank_from_feature),
                    e=r1_start + int(flank_from_feature))
            fasta_dict["feature_cut"] = "{c}:{s}".format(c=chrom,
                                                         s=r1_start)
            fasta_dict["r1_start_count"] = r1_count
            fasta_dict["reads_near_r1"] = non_cliff_reads_at_cliff
            fasta_dict["sequence"] = retrieve_sequence(fasta_dict["feature_region"], pyfaidx_fasta_obj)
            features.append(fasta_dict)
    sample_name = os.path.splitext(bam_file)[0]
    fasta_file = '{0}_candidate.fa'.format(sample_name)
    fasta_from_features(fasta_file, features, pyfaidx_fasta_obj)
    return features


### HELPER FUNCTIONS ###


def fasta_from_features(out_name, features, pyfaidx_fasta_obj):
    """
    This function requires a reference genome in a fasta file.
    See https://github.com/mdshw5/pyfaidx for more information.
    """
    region_dict = {}  ## avoid print repeat region  ## add in 20200226 by xuwei
    with open(out_name, "w+") as f:
        for feature in features:
            if feature["feature_region"] in region_dict.keys():
                continue
            region_dict[feature["feature_region"]] = feature["feature_region"]
            sequence = retrieve_sequence(feature["feature_region"],
                                         pyfaidx_fasta_obj)
            header = ">{region}|{cut}|{total_reads}|{r1_start_non_cliff}|{start_count}|{forward_reads}|{reverse_reads}".format(
                region=feature["feature_region"], cut=feature["feature_cut"],
                total_reads=feature["total_reads"],
                r1_start_non_cliff=feature["reads_near_r1"],
                start_count=feature["r1_start_count"],
                reverse_reads=feature["reverse_reads"],
                forward_reads=feature["forward_reads"])
            f.write(header + "\n")
            f.write(sequence + "\n")
    return out_name


def retrieve_sequence(coordinate, pyfaidx_fasta_obj):
    """
    Obtain sequence from a pyfaidx fasta object given a
    coordinate.
    """
    chrom, start, stop = parse_coordinate(coordinate)
    if start < 0:
        return pyfaidx_fasta_obj[chrom][0:stop].seq
    elif stop > len(pyfaidx_fasta_obj[chrom]):
        return pyfaidx_fasta_obj[chrom][start:len(pyfaidx_fasta_obj[chrom])].seq
    else:
        return pyfaidx_fasta_obj[chrom][start:stop].seq


def parse_coordinate(coordinate):
    """
    Parse a coordinate formatted like :
    chrZ:#-# OR chrZ:# where Z and # are placeholders.
    """
    coord_range_pattern = r'(.*):\d+-\d+$'
    coord_site_pattern = r'(.*):\d+$'
    if re.match(coord_range_pattern, coordinate):
        chrom = coordinate.split(':')[0]
        start = int(coordinate.split(':')[-1].split('-')[0])
        stop = int(coordinate.split(':')[-1].split('-')[1])
        return chrom, start, stop
    elif re.match(coord_site_pattern, coordinate):
        chrom = coordinate.split(':')[0]
        start = int(coordinate.split(':')[-1])
        return chrom, start
    else:
        print(coordinate)


""" Find actual sequences of potential off-target sites
"""
def output_alignments(out_name, target_sequence, features, mismatch_threshold):
    # dictionary to store the matched reads
    matched_dict = {}
    # dictionary to add read_count for each pair chromosome:start_position among matched reads
    reads_dict = {}

    for feature in features:
        sequence = feature["sequence"]
        region = feature["feature_region"]
        chrom, start, stop = parse_coordinate(region)
        cut = int(feature["position"])

        offtarget_sequence_no_bulge, mismatches, offtarget_sequence_length, chosen_alignment_strand_m, start_no_bulge, end_no_bulge, \
        realigned_target, \
        bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end = \
            alignSequences(target_sequence, sequence, max_score=mismatch_threshold)

        # get genomic coordinates of sequences
        mm_start, mm_end, b_start, b_end = '', '', '', ''
        if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '+':
            mm_start = start + int(start_no_bulge)
            mm_end = start + int(end_no_bulge)
        if offtarget_sequence_no_bulge and chosen_alignment_strand_m == '-':
            mm_start = stop - int(end_no_bulge)
            mm_end = stop - int(start_no_bulge)

        if bulged_offtarget_sequence and chosen_alignment_strand_b == '+':
            b_start = start + int(bulged_start)
            b_end = start + int(bulged_end)
        if bulged_offtarget_sequence and chosen_alignment_strand_b == '-':
            b_start = stop - int(bulged_end)
            b_end = stop - int(bulged_start)

        # define overall start, end and strand. For bed annotation purposes
        if offtarget_sequence_no_bulge:
            target_start_absolute = mm_start
            target_end_absolute = mm_end
            target_strand_absolute = chosen_alignment_strand_m
        elif not offtarget_sequence_no_bulge and bulged_offtarget_sequence:
            target_start_absolute = b_start
            target_end_absolute = b_end
            target_strand_absolute = chosen_alignment_strand_b
        else:
            target_start_absolute = cut
            target_end_absolute = cut
            target_strand_absolute = '*'
        
        name = chrom + ':' + str(target_start_absolute) + '-' + str(target_end_absolute)
        total_reads = int(feature["total_reads"])
        non_cliff_reads_at_cliff = int(feature["reads_near_r1"])
        cliff_read_count = int(feature["r1_start_count"])
        forward_reads = int(feature["forward_reads"])
        reverse_reads = int(feature["reverse_reads"])
        window_name = feature["feature_region"]
        
        if offtarget_sequence_no_bulge or bulged_offtarget_sequence:
            if target_strand_absolute == "+" and abs(target_end_absolute-cut)<=10:
                tag = chrom + ':' + str(target_start_absolute)
                if tag not in reads_dict.keys():
                    reads_dict[tag] = cliff_read_count
                    matched_dict[tag] = [chrom, target_start_absolute, target_end_absolute, name, 
                                         total_reads, non_cliff_reads_at_cliff, cliff_read_count,
                                         forward_reads, reverse_reads, target_strand_absolute, window_name,
                                         sequence, offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end, bulged_offtarget_sequence,
                                         length, score, substitutions,insertions,deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         target_sequence, realigned_target]
                else:
                    current_read_count = reads_dict[tag]
                    if cliff_read_count > current_read_count:
                        reads_dict[tag] = cliff_read_count
                        matched_dict[tag] = [chrom, target_start_absolute, target_end_absolute, name, 
                                             total_reads, non_cliff_reads_at_cliff, cliff_read_count,
                                             forward_reads, reverse_reads, target_strand_absolute, window_name,
                                             sequence, offtarget_sequence_no_bulge, mismatches,
                                             chosen_alignment_strand_m, mm_start, mm_end, bulged_offtarget_sequence,
                                             length, score, substitutions,insertions,deletions,
                                             chosen_alignment_strand_b, b_start, b_end,
                                             target_sequence, realigned_target]
            elif target_strand_absolute == "-" and abs(target_start_absolute-cut)<=10:
                tag = chrom + ':' + str(target_start_absolute)
                if tag not in reads_dict.keys():
                    reads_dict[tag] = cliff_read_count
                    matched_dict[tag] = [chrom, target_start_absolute, target_end_absolute, name, 
                                         total_reads, non_cliff_reads_at_cliff, cliff_read_count,
                                         forward_reads, reverse_reads, target_strand_absolute, window_name,
                                         sequence, offtarget_sequence_no_bulge, mismatches,
                                         chosen_alignment_strand_m, mm_start, mm_end, bulged_offtarget_sequence,
                                         length, score, substitutions,insertions,deletions,
                                         chosen_alignment_strand_b, b_start, b_end,
                                         target_sequence, realigned_target]
                else:
                    current_read_count = reads_dict[tag]
                    if cliff_read_count > current_read_count:
                        reads_dict[tag] = cliff_read_count
                        matched_dict[tag] = [chrom, target_start_absolute, target_end_absolute, name, 
                                             total_reads, non_cliff_reads_at_cliff, cliff_read_count,
                                             forward_reads, reverse_reads, target_strand_absolute, window_name,
                                             sequence, offtarget_sequence_no_bulge, mismatches,
                                             chosen_alignment_strand_m, mm_start, mm_end, bulged_offtarget_sequence,
                                             length, score, substitutions,insertions,deletions,
                                             chosen_alignment_strand_b, b_start, b_end,
                                             target_sequence, realigned_target]

    # Write matched table
    tags_sorted = matched_dict.keys()
    tags_sorted.sort()
    outfile_matched = '{0}_identified_matched.txt'.format(out_name)
    
    o1 = open(outfile_matched, 'w')
    print('#Chromosome', 'Start', 'End', 'Name', 'TotalReadCount', 'NonCliffReadCount', 
          'CliffReadCount', 'ForwardReadCount', 'ReverseReadCount', 'Strand', 
          'WindowName', 'WindowSequence', 'Site_SubstitutionsOnly.Sequence', 
          'Site_SubstitutionsOnly.NumSubstitutions', 'Site_SubstitutionsOnly.Strand',
          'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',
          'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',
          'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',
          'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',
          'TargetSequence', 'RealignedTargetSequence', sep='\t', file=o1)
    o1.close()
    
    offtarget_matched_dict = {}
    ontarget_matched_dict = {}
    
    for key in tags_sorted:
        row = matched_dict[key]
        reads = row[6]
        mismatch = row[13]
        if not mismatch and mismatch == 0:
            ontarget_matched_dict[key] = reads
            # o1 = open(outfile_matched, 'a')
            # print(*row, sep='\t', file=o1)
            # o1.close()
        else:
            offtarget_matched_dict[key] = reads

    with open(outfile_matched, 'a') as o1:
        for key in sorted(ontarget_matched_dict, key=ontarget_matched_dict.__getitem__, reverse=True):
            row = matched_dict[key]
            print(*row, sep='\t', file=o1)
        for key in sorted(offtarget_matched_dict, key=offtarget_matched_dict.__getitem__, reverse=True):
            row = matched_dict[key]
            print(*row, sep='\t', file=o1)

    return outfile_matched

    
"""
Given a targetsite and window, use a fuzzy regex to align the targetsite to
the window. Returns the best match(es).
"""
def alignSequences(targetsite_sequence, window_sequence, max_score):

    window_sequence = window_sequence.upper()
    query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)
    # Try both strands
    alignments_mm, alignments_bulge = list(), list()
    alignments_mm.append(('+', 'standard', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
    alignments_mm.append(('-', 'standard', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
    alignments_bulge.append(('+', 'gapped', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
    alignments_bulge.append(('-', 'gapped', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))
    
    lowest_distance_score, lowest_mismatch = 100, max_score + 1
    chosen_alignment_b, chosen_alignment_m, chosen_alignment_strand_b, chosen_alignment_strand_m = None, None, '', ''
    
    # Use regex to find the best match allowing only for mismatches
    for aln_m in alignments_mm:
        strand_m, alignment_type_m, match_m = aln_m
        if match_m != None:
            mismatches, insertions, deletions = match_m.fuzzy_counts
            if mismatches < lowest_mismatch:
                chosen_alignment_m = match_m
                chosen_alignment_strand_m = strand_m
                lowest_mismatch = mismatches
                
    # Use regex to find the best match allowing for gaps, so that its edit distance is strictly lower than the
    # total number of mismatches of the sequence founded (if any) allowing only for mismatches.
    for aln_b in alignments_bulge:
        strand_b, alignment_type_b, match_b = aln_b
        if match_b != None:
            substitutions, insertions, deletions = match_b.fuzzy_counts
            if insertions or deletions:
                distance_score = substitutions + (insertions + deletions) * 3
                edistance = substitutions + insertions + deletions
                if distance_score < lowest_distance_score and edistance < lowest_mismatch:
                    chosen_alignment_b = match_b
                    chosen_alignment_strand_b = strand_b
                    lowest_distance_score = distance_score
                    
    if chosen_alignment_m:
        offtarget_sequence_no_bulge = chosen_alignment_m.group()
        mismatches = chosen_alignment_m.fuzzy_counts[0]
        start_no_bulge = chosen_alignment_m.start()
        end_no_bulge = chosen_alignment_m.end()
    else:
        offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''
        
    bulged_offtarget_sequence, score, length, substitutions, insertions, deletions, bulged_start, bulged_end, realigned_target = \
        '', '', '', '', '', '', '', '', 'none'
    if chosen_alignment_b:
        realigned_target, bulged_offtarget_sequence = realignedSequences(targetsite_sequence, chosen_alignment_b, max_score)
        if bulged_offtarget_sequence:
            length = len(chosen_alignment_b.group())
            substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
            score = substitutions + (insertions + deletions) * 3
            bulged_start = chosen_alignment_b.start()
            bulged_end = chosen_alignment_b.end()
        else:
            chosen_alignment_strand_b = ''
    
    return [offtarget_sequence_no_bulge, mismatches, len(offtarget_sequence_no_bulge), chosen_alignment_strand_m, start_no_bulge, end_no_bulge,
            realigned_target,
            bulged_offtarget_sequence, length, score, substitutions, insertions, deletions, chosen_alignment_strand_b, bulged_start, bulged_end]
    
"""
Recreat A!!! sequence in the window_sequence that matches the conditions given for the fuzzy regex.
Currently only working for off-targets with at most one bulge !!!
"""
def realignedSequences(targetsite_sequence, chosen_alignment, errors):
    match_sequence = chosen_alignment.group()
    substitutions, insertions, deletions = chosen_alignment.fuzzy_counts
    
    # get the .fuzzy_counts associated to the matching sequence after adjusting fot indels, where 0 <= INS, DEL <= 1
    realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))
    
    # recreate the target sequence, with a '-' in the case of an DNA-bulge
    if insertions:
        targetsite_realignments = [targetsite_sequence[:i] + '-' + targetsite_sequence[i:] for i in range(1, len(targetsite_sequence))]
    else:
        targetsite_realignments = [targetsite_sequence]
        
    realigned_target_sequence, realigned_offtarget_sequence = None, '' # in case the matching sequence is not founded
    
    for seq in targetsite_realignments:
        # recreate off-target sequence (with a '-' in the case of an RNA-bulge) and pattern matching the realigned target sequence
        if deletions:
            match_realignments = [match_sequence[:i + 1] + '-' + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
            match_pattern = [match_sequence[:i + 1] + seq[i + 1] + match_sequence[i + 1:] for i in range(len(match_sequence) - 1)]
        else:
            match_realignments = match_pattern = [match_sequence]
        
        x = extendedPattern(seq, errors)
        for y_pattern, y_alignment in zip(match_pattern, match_realignments):
            m = regex.search(x, y_pattern, regex.BESTMATCH)
            if m and m.fuzzy_counts == realigned_fuzzy:
                realigned_target_sequence, realigned_offtarget_sequence = seq, y_alignment
    return realigned_target_sequence, realigned_offtarget_sequence

"""
Allow for '-' in our search, but do not allow insertions or deletions. 
"""
def extendedPattern(seq, errors):
    IUPAC_notation_regex_extended = {'N': '[ATCGN]','-': '[ATCGN]','Y': '[CTY]','R': '[AGR]','W': '[ATW]','S': '[CGS]','A': 'A','T': 'T','C': 'C','G': 'G'}
    realign_pattern = ''
    for c in seq:
        realign_pattern += IUPAC_notation_regex_extended[c]
    return '(?b:' + realign_pattern + ')' + '{{s<={0}}}'.format(errors)



""" Reverse complement DNA sequence
"""
def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])

def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
    seq = seq.upper()
    """
    Given a sequence with ambiguous base characters, returns a regex that matches for
    the explicit (unambiguous) base characters
    """
    IUPAC_notation_regex = {'N': '[ATCGN]',
                            'Y': '[CTY]',
                            'R': '[AGR]',
                            'W': '[ATW]',
                            'S': '[CGS]',
                            'A': 'A',
                            'T': 'T',
                            'C': 'C',
                            'G': 'G'}

    pattern = ''

    for c in seq:
        pattern += IUPAC_notation_regex[c] # N match ATCGN, any nucleotide in N of NGG is not as mismatch, mismatches in other place is <=7

    if lookahead:
        pattern = '(?b:' + pattern + ')'

    pattern_standard = pattern + '{{s<={0}}}'.format(errors)
    pattern_gap = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
    # print(pattern_standard)
    # print(pattern_gap)
    return pattern_standard, pattern_gap

def visualizeOfftargets(infile, outfile, title):

    #output_folder = os.path.dirname(outfile)
    #if not os.path.exists(output_folder):
    #    os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets, target_seq, total_seq = parseSitesFile(infile)
    
    most_reads = 5
    for j, seq in enumerate(offtargets):
        reads_count = len(str(seq['reads']))
        if reads_count >= most_reads:
            most_reads = reads_count

    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full', size=(box_size*len(target_seq)+55+most_reads*9, 66+(total_seq+1)*box_size))

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    # Draw ticks
    if target_seq.find('N') >= 0:
        p = target_seq.index('N')
        if p > len(target_seq) / 2:  # PAM on the right end
            tick_locations = [1, len(target_seq)] + range(p, len(target_seq))  # limits and PAM
            tick_locations += [x + p - 20 + 1 for x in range(p)[::10][1:]]  # intermediate values
            tick_locations = list(set(tick_locations))
            tick_locations.sort()
            tick_legend = [p, 10, 1] + ['P', 'A', 'M']
        else:
            tick_locations = range(2, 6) + [14, len(target_seq)]  # complementing PAM and limits
            tick_legend = ['P', 'A', 'M', '1', '10'] + [str(len(target_seq) - 4)]

        for x, y in zip(tick_locations, tick_legend):
            dwg.add(dwg.text(y, insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))
    else:
        tick_locations = [1, len(target_seq)]  # limits
        tick_locations += range(len(target_seq) + 1)[::10][1:]
        tick_locations.sort()
        for x in tick_locations:
            dwg.add(dwg.text(str(x), insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw reference sequence row
    for i, c in enumerate(target_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(target_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 1  # leave some extra space after the reference row
    line_number = 0  # keep track of plotted sequences
    for j, seq in enumerate(offtargets):
        realigned_target_seq = offtargets[j]['realigned_target_seq']
        no_bulge_offtarget_sequence = offtargets[j]['seq']
        bulge_offtarget_sequence = offtargets[j]['bulged_seq']

        if no_bulge_offtarget_sequence != '':
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(no_bulge_offtarget_sequence, target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
        if bulge_offtarget_sequence != '':
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(bulge_offtarget_sequence, realigned_target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(realigned_target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1

        if no_bulge_offtarget_sequence == '' or bulge_offtarget_sequence == '':
            reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(reads_text)
        else:
            reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(target_seq) + 1) + 20, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(reads_text)
            reads_text02 = dwg.text(u"\u007D", insert=(box_size * (len(target_seq) + 1) + 7, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:23px; font-family:Courier")
            dwg.add(reads_text02)
        
        if seq['chr'].find("PERV") >= 0:
            reads_text = dwg.text(str(seq['chr']), insert=(box_size * (len(target_seq) + 1) + 20 + 50, y_offset + box_size * (line_number + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(reads_text)
    dwg.save()

def parseSitesFile(infile):
    offtargets = []
    total_seq = 0
    ontarget = []
    alltargets = []
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            line = line.rstrip('\n')
            line_items = line.split('\t')
            chr = line_items[0]
            offtarget_reads = line_items[6]
            no_bulge_offtarget_sequence = line_items[12]
            substitutions = line_items[13]
            bulge_offtarget_sequence = line_items[17]
            target_seq = line_items[26]
            realigned_target_seq = line_items[27]

            if no_bulge_offtarget_sequence != '':
                total_seq += 1
            if bulge_offtarget_sequence != '':
                total_seq += 1
            if no_bulge_offtarget_sequence != '' or bulge_offtarget_sequence != '':
                if no_bulge_offtarget_sequence:
                    #total_seq += 1
                    if int(substitutions) == 0:
                        ontarget.append({'chr': chr.strip(),
                                         'seq': no_bulge_offtarget_sequence.strip(),
                                         'bulged_seq': bulge_offtarget_sequence.strip(),
                                         'reads': int(offtarget_reads.strip()),
                                         'target_seq': target_seq.strip(),
                                         'realigned_target_seq': realigned_target_seq.strip()
                                         })
                    else:
                        offtargets.append({'chr': chr.strip(),
                                           'seq': no_bulge_offtarget_sequence.strip(),
                                           'bulged_seq': bulge_offtarget_sequence.strip(),
                                           'reads': int(offtarget_reads.strip()),
                                           'target_seq': target_seq.strip(),
                                           'realigned_target_seq': realigned_target_seq.strip()
                                           })
                elif bulge_offtarget_sequence:
                    #total_seq += 1
                    offtargets.append({'chr': chr.strip(),
                                       'seq': no_bulge_offtarget_sequence.strip(),
                                       'bulged_seq': bulge_offtarget_sequence.strip(),
                                       'reads': int(offtarget_reads.strip()),
                                       'target_seq': target_seq.strip(),
                                       'realigned_target_seq': realigned_target_seq.strip()
                                       })
    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    alltargets = ontarget + offtargets
    return alltargets, target_seq, total_seq

def call_site_seq_features_form_file(candidate_fasta_file, genome):
    features = []

    pyfaidx_fasta_obj = Fasta(genome)

    with open(candidate_fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                fasta_dict = {}
                line = line.rstrip('\n')
                line_items = line.split('|')
                fasta_dict["total_reads"] = line_items[2]
                fasta_dict["reads_near_r1"] = line_items[3]
                fasta_dict["r1_start_count"] = line_items[4]
                fasta_dict["reverse_reads"] = line_items[5]
                fasta_dict["forward_reads"] = line_items[6]
                cut = line_items[1]
                fasta_dict["feature_cut"] = cut
                cut_items = fasta_dict["feature_cut"] = cut.split(':')
                fasta_dict["chrom"] = cut_items[0]
                fasta_dict["position"] = cut_items[1]
                region = line_items[0]
                region_items = region.split('>')
                fasta_dict["feature_region"] = region_items[1]
                fasta_dict["sequence"] = retrieve_sequence(fasta_dict["feature_region"], pyfaidx_fasta_obj).replace("\r","")
                features.append(fasta_dict)
    return features

def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--layout', '-l', help='The layout of the bam file.', default="SE")
    parser.add_argument('--target_sequence', '-t', help='sgRNA and PAM sequence.', required=True)
    parser.add_argument('--mismatch_threshold', '-m', help='maxmum mismatch threshold.', default=7, type=int)
    parser.add_argument('--initial_peaks_file', '-p', help='A list of peaks already found. These are strings formatted like\'chr#:#-#\'', default=None)
    parser.add_argument('--candidate_fasta_file', '-cf', help='candidate fasta file.', default=None)
    parser.add_argument('--matched_file', '-mf', help='identified matched file.', default=None)
    parser.add_argument('--bam_file', '-b', help='The path and name of a sorted and indexed bam file.', required=True)
    parser.add_argument('--genome', '-g', help='A reference genome.', required=True)
    parser.add_argument('--out_name', '-o', help='The name/path of the output fasta file. The fasta file contains information about the peak in the header as well as the sequence around', required=True)
    parser.add_argument('--min_dup', '-md', help='The number of reads that terminate at the exact same location (cliff_reads) needed for a potential feature to be considered.', default=5, type=int)
    parser.add_argument('--min_read_length', '-ml', help='The minimum read length required. Reads shorter than this will be removed.', default=60, type=int)
    parser.add_argument('--flank_from_feature', '-f', help='Sets the sequence around the SITE-Seq feature to return in the fasta file output by the function.', default=20, type=int)
    parser.add_argument('--required_cliff_ratio', '-r', help='The minimum ratio of cliff reads to cliff spanning reads required in order for a feature to be considered.', default=0.9, type=float)
    parser.add_argument('--banned_chroms', '-c', help='A list of chromosome names that should not be considered when attempting to identify features. This may need to be changed depending on the reference used.', default=["alt", "random", "Un", "chrM"])
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    if not args.matched_file is None:
        outfile_matched = args.matched_file
    else:
        if not args.candidate_fasta_file is None:
            features = call_site_seq_features_form_file(args.candidate_fasta_file, args.genome)
        else:
            features = call_site_seq_features(args.layout, 
                                              args.initial_peaks_file, 
                                              args.bam_file, 
                                              args.genome, 
                                              args.min_dup,
                                              args.min_read_length, 
                                              int(args.flank_from_feature), 
                                              args.required_cliff_ratio,
                                              args.banned_chroms)
        outfile_matched = output_alignments(args.out_name, args.target_sequence, features, args.mismatch_threshold)
    
    visualizeOfftargets(outfile_matched, args.out_name, args.out_name)
    
    outfile_matched_example = '{0}.example'.format(outfile_matched)
    os.system('head -16 {0} >{1}'.format(outfile_matched, outfile_matched_example))
    example_svg = '{0}_example'.format(args.out_name)
    visualizeOfftargets(outfile_matched_example, example_svg, args.out_name)


if __name__ == '__main__':
    main()
