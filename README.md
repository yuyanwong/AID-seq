# AID-seq
The test data and script for AID-seq
## Dependencies
- To establish the environment when running this pipeline, you can use this conda environment ```software/aidseq.yaml```

## command:
```
perl scripts/generate_pipeline.pl samplesheet.txt \ 
<raw_data_dir> <out_dir> \ 
<genome_version> <genome> \ 
<min_dup> <flank_from_feature> \ 
<banned_chroms> <nuc> \
<genome_dir> <scripts_dir> 
```
```
nohup sh pipeline.sh &
```
#### If multiple sgRNAs are to be analyzed for the same sample
```
perl scripts/prepare_config.pl sample <clean_data_dir> <out_dir> <genome> <genome_dir> <scripts_dir>
perl scripts/prepare_new.pl sample_config.txt
sh sample_run.sh &
perl scripts/pick.pl sample 
perl scripts/stat.pl sample pick
perl scripts/stat.pl sample all
```

## arguments:
-  ```samplesheet.txt ```

    You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline.
    
    | sample | HEK293 | Homo Sapiens | sample | T    | XXXXXXNGG
    |---|---|---|---|---|---|
    
   
- ```sample_sgRNA.txt```

    You will need to create a ```sample```_sgRNA.txt with information about the sample's sgRNAs before running the pooling pipeline.
    
    |sample |id|XXXXXXNGG|
    |---|---|---|
    
- ```<raw_data_dir>``` : The directory of raw data.
- ```<clean_data_dir>```: The directory of clean data.
- ```<out_dir>```: The output directory where the results will be saved.
- ```<genome_version>```: genome version (GRChg38).
- ```<genome>```: genome name（hg38）.
- ```<min_dup>```: The number of reads that terminate at the exact same location (cliff_reads) needed for a potential feature to be considered (5).
- ```<flank_from_feature>```: Sets the sequence around the AID-Seq feature to return in the fasta file output by the function (20).
- ```<banned_chroms>```: A list of chromosome names that should not be considered when attempting to identify features. This may need to be changed depending on the reference used (default).
- ```<nuc>```: Nucleases (cas9).
- ```<genome_dir>```: The local path to genome.
- ```<scripts_dir>```: The path to this pipeline's scripts.
