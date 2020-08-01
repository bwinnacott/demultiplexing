#!/usr/bin/env python
import numpy as np
import gzip
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(prog='mean_qscore_position',description='Generates a distribution of mean qscores by base position across records.')
parser.add_argument('-f','--in_file',type=str,action='store',help='Specifies input fastq file.')
parser.add_argument('-o','--out_plot',type=str,action='store',help='Specifies output filename for plot.')
parser.add_argument('--index',default=False,action='store_true',help='Specifies whether or not input file (-f) is an index file.')
args = parser.parse_args()

def convert_phred(letter):
    """Converts a single character (letter) into a phred score"""
    score = ord(letter) - 33
    
    return score

def mean_qscore_distribution():
    if args.index:
        position_array = np.zeros(8,dtype=int)
    else:
        position_array = np.zeros(101,dtype=int)
    LN = 1
    with gzip.open(args.in_file,'rt') as f:
        for line in f:
            if LN % 4 == 0:
                line = line.strip()
                for ind,char in enumerate(line):
                    position_array[ind] += convert_phred(char)
                LN += 1
            else:
                LN += 1
                continue

    total_records = LN/4
    position_array = position_array/total_records

    return position_array

def plot_distribution(array):
    plt.bar(range(len(array)),array)
    plt.xlabel('Base Position')
    plt.ylabel('Mean Quality Score')
    sequence = args.out_plot.split('.')[0]
    plt.title('Mean Quality Score by Base Position_{}'.format(sequence))
    plt.savefig(args.out_plot)

scores = mean_qscore_distribution()
plot_distribution(scores)