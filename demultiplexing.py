#!/usr/bin/env python
import argparse
import gzip
import itertools
from collections import defaultdict
import re
import os

def generate_index_dictionaries(index_file):
    '''This function takes in the index file as an argument and generates three dictionaries for use throughout this program.
    The first dictionary stores the index sequences as keys and a list of two 'None' objects as values. The second and third dictionaries
    will be built from the first dictionary. Dictionary 2 will consist of swapped index pair permutations as keys (using itertools) and values 
    will be set to 0. The third dictionary will be the dual matched indexes as keys and values also set to 0. Returns all dictionaries.'''
    # initiate default dictionary for referencing index info for writing out; use default dictionary so that appending to list can 
    # occur without initializing list for each new key
    index_reference = defaultdict(list)
    i = 0
    # open index file
    with open(index_file,'r') as f:
        for line in f:
            # don't include the header line
            if i == 0:
                i += 1
                continue
            # split the line into column elements and access index ID and sequence; add to dictionary
            else:
                line = line.strip().split('\t')
                index_reference[line[4]].extend([None,None])
                i += 1
    
    # uses dictionary comprehension to create two new dictionaries with tuples of index permutations as keys and counter as values;
    # one for swapped indexes and the other for dual matched; uses information from 'index_reference' keys
    swapped_permutations = {tup:0 for tup in itertools.permutations(index_reference.keys(),2)}
    dual_matched_dict = defaultdict(int)
    for key,_ in index_reference.items():
        dual_matched_dict[(key,key)] = 0

    return index_reference,swapped_permutations,dual_matched_dict

def rev_comp(sequence):
    '''This function serves to return the reverse complement of the "input_sequence" provided as argument.'''
    # reverse the input sequence
    rev_sequence = sequence[::-1]
    # use str.translate to change bases in 'rev_sequence' to complementary bases
    rev_complement = rev_sequence.translate(str.maketrans('ACGT','TGCA'))

    return rev_complement

def convert_phred(letter):
    """Converts a single character (letter) into a phred score"""
    # converts the ASCII symbol to phred equivalent (-33)
    score = ord(letter) - 33
    
    return score

def open_output_files(index_dict,directory):
    '''This function takes in a dictionary (keys are 24 unique index sequences, values are index numbers for each)
    and generates a file handle for the forward and reverse output files for each matched index. Two additional files 
    are opened for each of the last categories: index hopping and unknown indexes. These are also added to the dictionary 
    (categories as keys and file handles as values) for future reference. Total of 52 files are created.'''
    # checks if specified output directory exists, if not, create it
    if not os.path.isdir(directory):
        os.makedirs(directory)
    # iterate over the input dictionary (containing index sequences as keys and lists of size 2 containing 'None' as 
    # values; open file handles for outputting both foward and reverse records for each matched index
    for key,_ in index_dict.items():
        index_dict[key][0] = open('{}/{}_forward.fq'.format(directory,key),'a')
        index_dict[key][1] = open('{}/{}_reverse.fq'.format(directory,key),'a')
    # open four additional file handles for writing out to both forward and reverse for unknown and unmatched indexes
    unmatched_forward = open('{}/unmatched_forward.fq'.format(directory),'a')
    unmatched_reverse = open('{}/unmatched_reverse.fq'.format(directory),'a')
    unknown_forward = open('{}/unknown_lowqual_forward.fq'.format(directory),'a')
    unknown_reverse = open('{}/unknown_lowqual_reverse.fq'.format(directory),'a')
    # add these file handles to the dictionary
    index_dict['unknown'] = [unknown_forward,unknown_reverse]
    index_dict['unmatched'] = [unmatched_forward,unmatched_reverse]

def close_output_files(index_dict):
    '''Takes in the same dictionary used to generate file handles in 'create_output_files' and closes all file handles. 
    Closes the additional four file handles for writing out to unmatched and unknown files.'''
    # iterate over the input dictionary (containing index sequences as keys and file handles for both forward and reverse 
    # as values); close each file handle
    for val in index_dict.values():
        val[0].close()
        val[1].close()

def check_qscore(qscore_string,qual_cutoff,method):
    '''Function to either calculate the average qscore for a given quality score string 'qscore_string', or check whether 
    each base qscore in 'qscore_string' is above a certain threshold. User provides option; default is average. Returns boolean 
    depending on if quality string meets criteria based on the user defined quality score cutoff. Utilizes the 'convert_phred' 
    function to translate from ASCII to phred score.'''
    # runs this section if finding average quality score for index is desired
    if method == 'avg':
        # initiate counter for storing total score
        total_score = 0
        # iterate through each character in string and add converted value to score counter
        for char in qscore_string:
            total_score += convert_phred(char)

        # obtain the average score
        avg_score = total_score/len(qscore_string)
        # check if average score is at least cutoff 'qual_cutoff'
        if avg_score >= qual_cutoff:
            return True
        else:
            return False
    # runs this section if user wants to check whether each individual base quality score is above cutoff 'qual_cutoff'
    if method == 'ind':
        for char in qscore_string:
            score = convert_phred(char)
            if score >= qual_cutoff:
                continue
            else:
                return False
        return True

def add_index_header(index1_seq,index2_seq,read1_record,read2_record):
    '''Function that takes in the current record for each read as well as the sequences for each index (input for 'index2_seq' 
    should already be rev_comp!), and adds the index sequence pair (separated by '-') to each read's header. Modifies the header 
    in place, so returns nothing.'''
    # combine the indexes (ex: 'ATGC-ATGC')
    combined_index = index1_seq + '-' + index2_seq
    # add the combined index to end of each read's header line
    read1_record[0] = read1_record[0] + '_{}'.format(combined_index)
    read2_record[0] = read2_record[0] + '_{}'.format(combined_index)

def write_out_record(read1_record,read2_record,index_dict,swapped=False,unknown=False):
    '''Function used to write out two records (one for forward and one for reverse; contained in 'read1_record' and 'read2_record') to 
    their appropriate output files. 'index_dict' is used to reference the specific file handles associated with the current record's index 
    (optional, only for matched index records). 'swapped' and 'unknown' are arguments specifying whether the current record's indexes have
    been identified as mismatched or unknown/low quality. Does not return anything, just updates output files with new record info.'''
    # if 'swapped' is set to True, write out records to mismatched files
    if swapped:
        files = index_dict['unmatched']
        for line in read1_record:
            files[0].write(line + '\n')
        for line in read2_record:
            files[1].write(line + '\n')
    # if 'unknown' is set to True, write out records to unknown files
    elif unknown:
        files = index_dict['unknown']
        for line in read1_record:
            files[0].write(line + '\n')
        for line in read2_record:
            files[1].write(line + '\n')
    # if neither 'swapped' nor 'unknown' are specified as arguments, records are matched based on index and wrote out to index specific files
    else:
        # extract the specific index associated with the current records
        index = re.search(r'_[A-Z]+',read1_record[0])
        # obtain the appropriate files handles based on 'index' from 'index_dict'
        files = index_dict[index.group()[1:]]
        for line in read1_record:
            files[0].write(line + '\n')
        for line in read2_record:
            files[1].write(line + '\n')

def parse_files(index_reference,swapped_permutations,dual_matched_dict,cutoff,method,i1_file,i2_file,br1_file,br2_file):
    '''Goes through the four specified input FASTQ files (2 index and 2 biological sequences) and assigns each set of records to 
    a new fq file based on the record's index category designation. There are three categories: 1) matched indexes, 2) swapped indexes,
    and 3) unknown/low quality indexes.'''
    # open four input files (gzipped files)
    with gzip.open(br1_file,'rt') as read1, gzip.open(i1_file,'rt') as index1, gzip.open(i2_file,'rt') as index2, gzip.open(br2_file,'rt') as read2:
        # initiate variables used to send records to correct output files (line count, counters for each index category, record information 
        # holders, and category indicator booleans)
        LN = 0
        total_records,total_matched,total_unmatched,total_unknown = 0,0,0,0
        read_forward = []
        read_reverse = []
        swapped,unknown = False,False
        # iterate over each line in all four files
        for read1_line in read1:
            read1_line = read1_line.strip()
            read2_line = read2.readline().strip()
            index1_line = index1.readline().strip()
            index2_line = index2.readline().strip()
            # enter if when sequence line is reached for each record
            if LN % 4 == 1:
                # get reverse complement of index2
                index2_rev_comp = rev_comp(index2_line)
                # check if both indexes are any of unique indexes used for library prep
                if index1_line in index_reference and index2_rev_comp in index_reference:
                    # check to see if the indexes are not dual matched, if so, adjust 'unmatched' variables; increment swapped dict counter
                    if index1_line != index2_rev_comp:
                        swapped = True
                        total_unmatched += 1
                        swapped_permutations[(index1_line,index2_rev_comp)] += 1
                    # if matched, increment matched dict counter
                    else:
                        dual_matched_dict[(index1_line,index2_rev_comp)] += 1
                # if one or more of the indexes are not in the index file, adjust 'unknown' variables
                else:
                    unknown = True
                    total_unknown += 1
                # append sequence lines (only for biological reads!) to the record information list variables
                read_forward.append(read1_line)
                read_reverse.append(read2_line)
                # call 'add_index_header' to add the index sequences from both forward and reverse reads to header of biological reads
                add_index_header(index1_line,index2_rev_comp,read_forward,read_reverse)
                # increment line count
                LN += 1
            # enter if when qscore line is reached for each record
            elif LN % 4 == 3:
                # append this line (only for biological reads) to the lists containing the record information
                read_forward.append(read1_line)
                read_reverse.append(read2_line)
                # if 'swapped' variable is True (record identified as having swapped indexes), enter if statement
                if swapped:
                    # check the quality score based on methods available in function 'main' (see 'help')
                    if check_qscore(index1_line,cutoff,method) and check_qscore(index2_line,cutoff,method):
                        # if the quality is acceptable, write out the current biological records to the unmatched index file
                        write_out_record(read_forward,read_reverse,index_reference,swapped=True)
                    # if quality score is below cutoff, change variables accordingly (no longer considered swapped, considered unknown)
                    else:
                        swapped = False
                        unknown = True
                        total_unmatched -= 1
                        total_unknown += 1
                        # decrement previously recorded occurrence in swapped permutation dictionary (no longer considered swapped; instead
                        # considered unknown)
                        index_perm = re.search(r'[A-Z]+-[A-Z]+',read_forward[0])
                        swapped_permutations[tuple(index_perm.group().split('-'))] -= 1
                # if 'unknown' variable is True (record identified as having unknown or low quality index sequences), write out current 
                # records to unknown files
                if unknown:
                    write_out_record(read_forward,read_reverse,index_reference,unknown=True)
                # if the current sequences are found to be dual matched, enter if
                if not swapped and not unknown:
                    # check the quality score based on methods available in function 'main' (see 'help')
                    if check_qscore(index1_line,cutoff,method) and check_qscore(index2_line,cutoff,method):
                        # if quality is acceptable, write out the current biological records to matched records based on specific index
                        total_matched += 1
                        write_out_record(read_forward,read_reverse,index_reference)
                    # if low quality, index is considered unknown (no longer matched), increment unknown counter, and write out records to 
                    # unknown files
                    else:
                        total_unknown += 1
                        write_out_record(read_forward,read_reverse,index_reference,unknown=True)
                        # decrement previously recorded occurrence in matched permutation dictionary (no longer considered matched; instead
                        # considered unknown)
                        index_perm = re.search(r'[A-Z]+-[A-Z]+',read_forward[0])
                        dual_matched_dict[tuple(index_perm.group().split('-'))] -= 1
                # reset variables for next record, increment line counter
                swapped,unknown = False,False
                read_forward = []
                read_reverse = []
                LN += 1
                total_records += 1
            # if other two lines, append to lists and increment line count
            else:
                read_forward.append(read1_line)
                read_reverse.append(read2_line)
                LN += 1
    
    return total_records,total_matched,total_unmatched,total_unknown

def generate_report(swapped_permutations,dual_matched_dict,total_records,num_matched,num_unmatched,num_unknown):
    '''Generates a comprehensive report of total number of matched, mismatched, and unknown index events. Additionally, a table detailing 
    specific dual matched and index swapping occurrences (and percentages) is provided.'''
    # open report file to write results out to (includes totals for each index category as well as percentage of total records for each)
    with open('index_report.txt','w') as out_f:
        out_f.write('Total number of records processed: ' + str(total_records) + '\n')
        out_f.write('Total number of matched indexes: ' + str(num_matched) + '\n')
        out_f.write('Percentage of matched indexes: ' + str(round((num_matched/total_records)*100,2)) + '\n')
        out_f.write('Total number of swapped indexes: ' + str(num_unmatched) + '\n')
        out_f.write('Percentage of swapped indexes: ' + str(round((num_unmatched/total_records)*100,2)) + '\n')
        out_f.write('Total number of unknown indexes: ' + str(num_unknown) + '\n')
        out_f.write('Percentage of unknown indexes: ' + str(round((num_unknown/total_records)*100,2)) + '\n')
        out_f.write('\n' + 'Summary of Dual Matched Indexes' + '\n')
        out_f.write('\n' + '|' + ' Index Pair ' + ' '*14 + '|' + ' Number of Occurrences ' + '|' + ' Percent Occurence ' + '|' + '\n')
        # create a table to display the dual matched index data separate from swapped
        for key,val in dual_matched_dict.items():
            out_f.write('| ' + str(key) + ' | ' + str(val) + ' '*(22-len(str(val))) + '| ' + str(round((val/total_records)*100,3)) + ' '*(18-len(str(round((val/total_records)*100,3)))) + '|' + '\n')
        out_f.write('\n' + 'Summary of Swapped Indexes' + '\n')
        out_f.write('\n' + '|' + ' Index Pair ' + ' '*14 + '|' + ' Number of Occurrences ' + '|' + ' Percent Occurence ' + '|' + '\n')
        # create an additional table to report swapped index data
        for key,val in swapped_permutations.items():
            out_f.write('| ' + str(key) + ' | ' + str(val) + ' '*(22-len(str(val))) + '| ' + str(round((val/total_records)*100,3)) + ' '*(18-len(str(round((val/total_records)*100,3)))) + '|' + '\n')

def main():
    # initiate command line arguments/options
    parser = argparse.ArgumentParser(description='Demultiplexes FASTQ sequencing files (separate index and biological read files) \
                                    and outputs results report.')
    parser.add_argument('-ir1','--ind1_seq_file',action='store',type=str,help='Specifies file containing sequence reads for index 1.')
    parser.add_argument('-ir2','--ind2_seq_file',action='store',type=str,help='Specifies file containing sequence reads for index 2.')
    parser.add_argument('-br1','--read1_seq_file',action='store',type=str,help='Specifies file containing sequence reads for biological read 1.')
    parser.add_argument('-br2','--read2_seq_file',action='store',type=str,help='Specifies file containing sequence reads for biological read 2.')
    parser.add_argument('-ind_file','--index_file',action='store',type=str,help='Specifies file containing list of unique indexes (tab separated).')
    parser.add_argument('-q','--quality_cutoff',action='store',type=int,default=30,help='Specifies quality score cutoff used to filter index \
                        reads (and biological record information) out.')
    parser.add_argument('-m','--method',action='store',choices=['avg','ind'],help='If "avg", tells program to calculate an average score across \
                        the index quality score string, for which to compare specified cutoff (-q) value to. If "ind", calculates each \
                        individual base quality score separately and compares against the cutoff (-q). If ANY base is below the user-defined \
                        cutoff, then record is thrown out.')
    parser.add_argument('-outdir','--output_dir',action='store',type=str,help='Specifies output directory for which to store files. If directory \
                        doesnt exist, its created.')
    args = parser.parse_args()
    # runs the program
    index_reference,swapped_permutations,dual_matched_dict = generate_index_dictionaries(args.index_file)
    open_output_files(index_reference,args.output_dir)
    total,matched,unmatched,unknown = parse_files(index_reference,swapped_permutations,dual_matched_dict,args.quality_cutoff,args.method,args.ind1_seq_file,
                                                args.ind2_seq_file,args.read1_seq_file,args.read2_seq_file)
    close_output_files(index_reference)
    generate_report(swapped_permutations,dual_matched_dict,total,matched,unmatched,unknown)

if __name__ == "__main__":
    main()
