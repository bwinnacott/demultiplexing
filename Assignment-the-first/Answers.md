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
    is applied, which bins reads  based on their index information. However, a phenomenon known as index hopping can arise from sequencing runs that 
    utilize multiplexing. This occurs most frequently on patterned flow cells where the library fragment from one well anneals to a different fragment 
    (w/different index) in another nearby well. Additionally, having excess free floating adapters in solution can increase the chances of index hopping. 
    If not accounted for, this can result in complications in downstream analysis. For this assignment, our goal is to demultiplex paired end sequencing 
    data that may have experienced index misassignment. 
    ```
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
