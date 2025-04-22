######## Compare sequence with HG002 assembly ########
import pandas as pd
import numpy as np
import math
import pysam
import sys
from tqdm import tqdm


n = sys.argv[1]
input = sys.argv[2]

HG002_bam = pysam.AlignmentFile("/projects/ps-gymreklab/helia/HipSTR_LR/TRGT/sequence_compare/HG002_assembly_hg38_sorted.bam", "rb")

def extract_alleles_from_assembly(chrom, repeat_pos, repeat_end, padding):
    iter_ = HG002_bam.fetch(chrom, repeat_pos, repeat_end) #fetching reads overlapping repeat
    
    alleles = []
    info = []
    maternal_seen = 0
    paternal_seen = 0
    for read in iter_:
        if maternal_seen == 1 and paternal_seen == 1:
            break
        if read.reference_end < repeat_end or read.reference_start > repeat_pos: # not spanning
            continue
        if read.query_name == chrom + "_MATERNAL":
            maternal_seen = 1
        elif read.query_name == chrom + "_PATERNAL":
            paternal_seen = 1
        else:
            continue
        cigar = read.cigartuples
        cigar_list = []
        for c in cigar: # convert list of tuples to list of lists
            cigar_list.append([c[0],c[1]])            

        cigar_cnt = 0
        current_pos = read.reference_start
        ltrim = 0 # how many characters should be cut from left
        while True:
            #print(cigar_cnt, cigar_list[cigar_cnt])
            if cigar_cnt >= len(cigar_list):
                break
            if current_pos == repeat_pos - 1:
                #print(cigar_list[cigar_cnt+1:cigar_cnt+10])
                break
            current_cigar = cigar_list[cigar_cnt]
            
            if int(current_cigar[0]) in [0,7,8]:
                if current_pos + current_cigar[1] < repeat_pos - padding:
                    current_pos += current_cigar[1]
                    cigar_cnt += 1
                    ltrim += current_cigar[1]
                    continue
                current_pos += 1
                ltrim += 1
            elif current_cigar[0] == 2:
                if (current_pos >= repeat_pos - padding) and (current_pos < repeat_pos):
                    info.append(("BEFORE_DEL", current_cigar[1], repeat_pos - current_pos))  
                if current_pos + current_cigar[1] < repeat_pos - padding:
                    current_pos += current_cigar[1]
                    cigar_cnt += 1
                    continue
                current_pos += 1
            elif current_cigar[0] in [1,4]:
                if (current_pos >= repeat_pos - padding) and (current_pos < repeat_pos):
                    info.append(("BEFORE_INS", current_cigar[1], repeat_pos - current_pos))  
                ltrim += current_cigar[1]
                cigar_cnt += 1
                continue
            elif current_cigar[0] == 5:
                cigar_cnt += 1
                continue
            current_cigar[1] -= 1
            if current_cigar[1] == 0:
                cigar_cnt +=1

        cigar_cnt = len(cigar_list) - 1
        current_pos = read.reference_end
        rtrim = 0 # how many characters should be cut from right
        while True:
            #print("hh", cigar_cnt)
            if cigar_cnt <= -1:
                break
            if current_pos == repeat_end:
                break
            current_cigar = cigar_list[cigar_cnt]
            if current_cigar[0] in [0,7,8]:
                if current_pos - current_cigar[1] > repeat_end + padding:
                    current_pos -= current_cigar[1]
                    cigar_cnt -=1
                    rtrim += current_cigar[1]
                    continue
                current_pos -= 1
                rtrim += 1
            elif current_cigar[0] == 2:
                if (current_pos > repeat_end) and (current_pos <= repeat_end + padding):
                    info.append(("AFTER_DEL", current_cigar[1], current_pos - repeat_end))    
                if current_pos - current_cigar[1] > repeat_end + padding:
                    current_pos -= current_cigar[1]
                    cigar_cnt -= 1
                    continue
                current_pos -= 1
            elif current_cigar[0] in [1,4]:
                if (current_pos > repeat_end) and (current_pos <= repeat_end + padding):
                    info.append(("AFTER_INS", current_cigar[1], current_pos - repeat_end))                    
                rtrim += current_cigar[1]
                cigar_cnt -= 1
                continue
            elif current_cigar[0] == 5:
                cigar_cnt -= 1
                continue
            current_cigar[1] -= 1
            if current_cigar[1] == 0:
                cigar_cnt -=1

        seq = read.query_sequence
        if seq == None:
            continue
        alleles.append(seq[ltrim: len(seq) - rtrim])
    if chrom != 'chrX' and chrom != 'chrY' and (maternal_seen == 0 or paternal_seen == 0):
        return [], ["no alignment from both haplotypes"]
        
    return alleles, info


different_seq_df = pd.read_csv(input, header=None, delim_whitespace=True)
different_seq_df = different_seq_df[different_seq_df[0] == "chr" + n]

for index, row in tqdm(different_seq_df.iterrows()):
    alleles,info = extract_alleles_from_assembly(row[0], int(row[1]), int(row[2]), 25)
    alleles = sorted(alleles)

    print("\t".join([row[0], str(row[1]), str(row[2]), str(row[5]), str(alleles), str(info)]))


    
