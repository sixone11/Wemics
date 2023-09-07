#!/bin/env python3
import os,sys
import numpy as np
import pandas as pd
import argparse
from utils import *
from function import *

# CpG context
CG = ['Z','z']
CHG = ['X','x']
CHH = ['H','h']
XXX = ['U','u']
C_context_include = []
C_context = ['CG']

if "CG" in C_context: C_context_include.extend(CG)
if "CHG" in C_context: C_context_include.extend(CHG)
if "CHH" in C_context: C_context_include.extend(CHH)
if 'XXX' in C_context: C_context_include.extend(XXX)
C_context_exclude = list(set(['Z','z','H','h','X','x','U','u']).difference(set(C_context_include)))

# debug variable
infile = sys.argv[1]
contig = sys.argv[2]
overlap = False
outfile = sys.argv[3]

reads = pd.read_csv(infile,sep="\t",header=None)
reads.columns = ["chr","read1_start","read1_strand","read1_cigar","read1","read2_start","read2_strand","read2_cigar","read2"]
data = reads[reads["chr"] == contig]
data.reset_index(drop=True, inplace=True)

def correct_read_cigar(cigar,reads_Cseq):
    cigar_list = re.split('(\d+)',cigar)
    for C_context in C_context_exclude: 
        reads_Cseq = reads_Cseq.replace(C_context,".")
#     # read doesn't contain C or G
#     if not reads_Cseq.replace(".","") : 
#         return False 
    if len(cigar_list) == 3 and cigar_list[2] == "M" :
            seq = reads_Cseq 
    else:
        pos = 0 
        seq = ""
        for process in range(2,len(cigar_list),2):
            # 0 MEANS: MATCH M
            cigar_list[process-1] = int(cigar_list[process-1])
            if cigar_list[process] == "M" :
                # pdb.set_trace()
                seq += reads_Cseq[pos:(pos + cigar_list[process-1])]
                pos += cigar_list[process-1]
            # 1 MEANS: INSERT I
            elif cigar_list[process] == "I" : 
                pos += cigar_list[process-1]
            # 2 MEANS: DEL D
            elif cigar_list[process] == "D" :
                seq += "." * cigar_list[process-1]
            # 3 MEANS： similarlly like D, skipped N bases
            elif cigar_list[process] == "N" :
                seq += "." * cigar_list[process-1]
            else:
                print("There was a character not MIDN in CIGAR!") 
                return False
    return seq

## 00 Initialization
# disp("process contig " + args.contig + " start")
contig_methylation = {}
contig_length = 250000000
for element in ['coverage','meth','CHALM_meth','CAMDA_meth']:
    # contig_methylation[element] = np.zeros(contigs_length[args.contig], dtype=int)   
    contig_methylation[element] = np.zeros(contig_length, dtype=int) 


start_pos_list = list()
end_pos_list = list()
seq_list = list()
strand_list = list()
# C_incontext_dis_up_list = list()
# C_incontext_dis_down_list = list()
C_incontext_pos_list = list()


for index in range(0,data.shape[0]):
    if index % 10000 == 0:
            disp("hit 10000 times")
    read1_seq = correct_read_cigar(data["read1_cigar"][index],data["read1"][index])
    read2_seq = correct_read_cigar(data["read2_cigar"][index],data["read2"][index])
    read1_start = data["read1_start"][index]
    read2_start = data["read2_start"][index]
    seq_read1_read2 = {'read1':read1_seq,'read2':read2_seq}
    start_read1_read2 = {'read1':read1_start,'read2':read2_start}
    if read1_start <= read2_start :
        read_front_str = "read1"
        read_behind_str = "read2"
    else:
        read_front_str = "read2"
        read_behind_str = "read1"
    read_len = len(seq_read1_read2[read_front_str])
    overlap = read_len + start_read1_read2[read_front_str] - start_read1_read2[read_behind_str]
    # read position1: overlap 
    if overlap > 0:
        # using the calls from the first read which is presumably the one with a lowest error rate
        if read_front_str == "read1": 
            seq = read1_seq + read2_seq[overlap:]
        else:
            gap1_end = read_len - overlap 
            # gap1 = seq_read1_read2[read_front_str][0:gap1_end]
            seq = read2_seq[0:gap1_end] + read1_seq
    # read position2: no overlap between read and mate 
    else:
        gap_dis = -overlap
        gap2 = "." * gap_dis
        seq = seq_read1_read2[read_front_str] + gap2 + seq_read1_read2[read_behind_str]
    # bam file and python and bed file are different based

    start_pos = start_read1_read2[read_front_str] - 1 
    # end_pos = start_pos + len(seq) - 1 
    # 0-base for further python script use
    C_incontext_pos = [ start_pos + index for index in range(0,len(seq)) if seq[index].isalpha() ]
    if len(C_incontext_pos) != 0 :
        C_incontext_pos_list.append(C_incontext_pos)
        for C_context in C_context_exclude:
            seq = seq.replace(C_context,".")
        # relative pos
        C_incontext_pos = [ index for index in range(0,len(seq)) if seq[index].isalpha() ]
        # 0-based
        start_pos_list.append(start_pos)
        # 1-based
        end_pos_list.append(start_pos + len(seq))
        # calculate distance
        # C_incontext_pos_up = C_incontext_pos[1:]
        # C_incontext_pos_down = C_incontext_pos[:-1]
        # C_incontext_dis = [  C_incontext_pos_up[index] - C_incontext_pos_down[index] for index in range(0,(len(C_incontext_pos_up)))]
        # C_incontext_dis.insert(0,np.nan)
        # C_incontext_dis.append(np.nan)
        # C_incontext_dis_up = C_incontext_dis[:-1]
        # C_incontext_dis_up_list.append(C_incontext_dis_up)
        # C_incontext_dis_down = C_incontext_dis[1:]
        # C_incontext_dis_down_list.append(C_incontext_dis_down)
        seq = seq.replace(".","")
        for C_context in C_context_include:
            if C_context.isupper():
                seq = seq.replace(C_context,"1")
            else:
                seq = seq.replace(C_context,"0")
        seq_list.append([int(char) for char in seq])
        strand_list.append(data["read1_strand"][index])

# out = pd.DataFrame({"chr":contig,"start":start_pos_list,"end":end_pos_list,"strand":strand_list,"seq":seq_list,"C_pos":C_incontext_pos_list,"dis_up":C_incontext_dis_up_list,"dis_down":C_incontext_dis_down_list})
out = pd.DataFrame({"chr":contig,"start":start_pos_list,"end":end_pos_list,"strand":strand_list,"seq":seq_list,"C_pos":C_incontext_pos_list})

# reads must include C incontext pos
out = out[out["seq"] != ""]
out.reset_index(drop=True, inplace=True)
out.to_csv(outfile,index=False,sep="\t",encoding='gbk',na_rep='NaN',header=None)
