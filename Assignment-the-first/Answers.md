# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. ![](https://github.com/2020-bgmp/demultiplexing-bwinnacott/blob/master/plots/read1.png?raw=true)
       ![](https://github.com/2020-bgmp/demultiplexing-bwinnacott/blob/master/plots/index1.png?raw=true)
       ![](https://github.com/2020-bgmp/demultiplexing-bwinnacott/blob/master/plots/index2.png?raw=true)
       ![](https://github.com/2020-bgmp/demultiplexing-bwinnacott/blob/master/plots/read2.png?raw=true)
       
    2.
    ```
    An appropriate average quality score cutoff for the index sequences is at least Q30 (1 in 1,000 probability of incorrect base call). We want to make sure the 
    index sequences are of high quality due to their importance in assigning reads to their correct libraries. If misassignment occurs, downstream analysis 
    can be compromised. It may even be more prudent to apply a Q30 cutoff for each individual base call in each index sequence. However, this might remove
    an excess of good data that can support further analysis. For biological reads, a slightly lower quality score cutoff might be appropriate, such as Q25, to 
    retain more sequence information, without compromising quality. Having higher coverage with shorter read sequencing experiments provides overlap which 
    can resolve consensus sequences and error correct. 
    ```
    
    3. Command: 
    
       ```
       ls -1 1294_S1_L008_R[23].001.fastq.gz | while read line; do echo $line; zcat $line | sed -n "2~4p" | grep 'N' | wc -l; done
       ```
       
       Output: 
       
       ```
       1294_S1_L008_R2_001.fastq.gz
       3976613
       1294_S1_L008_R3_001.fastq.gz
       3328051
       ```
    
## Part 2
1. Define the problem
    ```
    Multiplexing is used to generate sequencing data for multiple libraries simultaneously on a single flow cell. This is accomplished using adapters 
    containing unique index sequences (barcodes) that map reads back to their respective libraries. To accomplish this, a technique called demultiplexing 
    is applied, which bins reads based on their index information. However, a phenomenon known as index hopping can arise from sequencing runs that 
    utilize multiplexing. This occurs most frequently on patterned flow cells where the library fragment from one well anneals to a different fragment 
    (w/different index) in another nearby well. Additionally, having excess free floating adapters in solution can increase the chances of index hopping. 
    If not accounted for, this can result in complications in downstream analysis. For this assignment, our goal is to demultiplex paired end sequencing 
    data that may have experienced index misassignment. 
    ```
2. Describe output
   ```
   Output from demultiplexing consists of fastq files that contain categorized reads from the original raw fastq files. Two reads for each unique
   index (forward and reverse) will be generated, as well as two additional output files (forward and reverse) for each of the remaining index categories:
   index swapped and low quality/unknown indexes. For this assignment, a total of 52 output fastq files will be generated (24 unique indexes x 2 = 48 + 2
   (swapped) = 50 + 2 (low quality/unknown) = 52).
   ```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
   ```
   See ../TEST-input_FASTQ and ../TEST-output_FASTQ
   ```
4. Pseudocode
   ```
   See pseudocode.txt
   ```
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
   ```
   See pseudocode.txt
   ```
