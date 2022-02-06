# Demultiplexing Tool

## Overview
This tool demultiplexes pooled sample fastq files. It outputs records to new fastq files 
based on their index designation, and whether the indexes for each read were dual-matched. 
Reads with indexes found to be swapped are output to separate files, and records with unknown 
indexes (due to sequencing error) are kept separate as well. A report summarizing index designations 
is provided for the user on completion. 

## Requirements
- Python 3.6+

## Input File Requirements
- **Fastq files**: Read 1 and read 2 biological fastq files
- **Index files**: Read 1 and read 2 index fastq files
- **Index sequence list**: File containing list of unique index sequences

## Arguments
The following arguments are required for running the program:
- ```-ir1```, ```--ind1_seq_file```: specifies absolute path to input FASTA file containing read 1 index 
sequences
- ```-ir2```, ```--ind2_seq_file```: specifies absolute path to input FASTA file containing read 2 index  
sequences
- ```-br1```, ```--read1_seq_file```: specifies absolute path to input FASTA file containing read 1 biological
sequences
- ```-br2```, ```--read2_seq_file```: specifies absolute path to input FASTA file containing read 2 biological
sequences
- ```-ind_file```, ```--index_file```: specifies absolute path to file containing list of unique index sequences

The following arguments are optional:
- ```-q```, ```--quality_cutoff```: quality score cutoff used to filter index and biological sequences reads out
- ```-m```, ```--method```: method used to calculate index sequence quality; use *avg* to calculate and average 
score across the index quality score string, or use *ind* to calculate each individual base quality score separately 
and compare against the threshold (record is thrown out if any base quality score is below threshold in this mode)
- ```-outdir```, ```--output_dir```: output directory to store new fastq files and report