#!/usr/bin/env python3
# -*- coding: utf-8 -*-
############################################################################
import numpy as np
############################################################################
def count_read_meth(read_seq,read_start,CHALM_threshold,contig_methylation):
    start_pos = read_start - 1
    methylated_pos= [ start_pos + index  for index in range(0,len(read_seq)) if read_seq[index].isupper() ] # and refmark[args.contig][start_pos + index] in C_context_num_include 
    C_incontext_pos = [ start_pos + index  for index in range(0,len(read_seq)) if read_seq[index].isalpha() ]
    methylated_num = len(methylated_pos)
    totalC_num = len(C_incontext_pos)
    if totalC_num != 0:
        # for coverage record
        np.put(contig_methylation['coverage'], C_incontext_pos,contig_methylation['coverage'][C_incontext_pos] + 1)
        # for methy record
        if methylated_num != 0 :
            # 01 trad mean methy record
            np.put(contig_methylation['meth'],methylated_pos,contig_methylation['meth'][methylated_pos] + 1)
            # 02 CHALM methy record
            if methylated_num >= CHALM_threshold:
                CHALM_methylated_pos = C_incontext_pos
            else :
                CHALM_methylated_pos = methylated_pos
            np.put(contig_methylation['CHALM_meth'],CHALM_methylated_pos,contig_methylation['CHALM_meth'][CHALM_methylated_pos] + 1)
            # 03 CAMDA methy record
            if methylated_num != totalC_num :
                CAMDA_unmethylated_pos = [ start_pos + index for index in range(0,len(read_seq)) if read_seq[index].islower() ]
                np.put(contig_methylation['CAMDA_meth'],CAMDA_unmethylated_pos,contig_methylation['CAMDA_meth'][CAMDA_unmethylated_pos] + 1)
    return contig_methylation
############################################################################
def count_part_read2_meth(read_seq,read_start,scan_start,scan_end,CHALM_threshold,contig_methylation):
    start_pos = read_start - 1
    methylated_pos= [ start_pos + index  for index in range(scan_start,scan_end) if read_seq[index].isupper() ] # and refmark[args.contig][start_pos + index] in C_context_num_include 
    C_incontext_pos = [ start_pos + index  for index in range(scan_start,scan_end) if read_seq[index].isalpha() ]
    methylated_num = len(methylated_pos)
    totalC_num = len(C_incontext_pos)
    # for coverage record
    if totalC_num!=0:
        np.put(contig_methylation['coverage'], C_incontext_pos,contig_methylation['coverage'][C_incontext_pos] + 1)
        # for methy record
        if methylated_num != 0 :
            # 01 trad mean methy record
            np.put(contig_methylation['meth'],methylated_pos,contig_methylation['meth'][methylated_pos] + 1)
            # 02 CHALM methy record
            if methylated_num >= CHALM_threshold:
                CHALM_methylated_pos = C_incontext_pos
            else :
                CHALM_methylated_pos = methylated_pos
            np.put(contig_methylation['CHALM_meth'],CHALM_methylated_pos,contig_methylation['CHALM_meth'][CHALM_methylated_pos] + 1)
            # 03 CAMDA methy record
            if methylated_num != totalC_num :
                CAMDA_unmethylated_pos = [ start_pos + index for index in range(scan_start,scan_end) if read_seq[index].islower() ]
                np.put(contig_methylation['CAMDA_meth'],CAMDA_unmethylated_pos,contig_methylation['CAMDA_meth'][CAMDA_unmethylated_pos] + 1)
    return contig_methylation
############################################################################
def count_whole_read2_meth(read_seq,read_start,scan_start,scan_end,CHALM_threshold,contig_methylation):
    start_pos = read_start - 1
    methylated_pos= [ start_pos + index  for index in range(0,len(read_seq)) if read_seq[index].isupper() ] # and refmark[args.contig][start_pos + index] in C_context_num_include 
    C_incontext_pos = [ start_pos + index  for index in range(0,len(read_seq)) if read_seq[index].isalpha() ]
    part_methylated_pos = [ index for index in methylated_pos if index >= start_pos + scan_start and index < start_pos + scan_end ]
    part_C_incontext_pos = [ index for index in C_incontext_pos if index >= start_pos + scan_start and index < start_pos + scan_end ]
    methylated_num = len(methylated_pos)
    totalC_num = len(C_incontext_pos)
    part_methylated_num = len(part_methylated_pos)
    part_totalC_num = len(part_C_incontext_pos)
    # for coverage record
    if part_totalC_num != 0:
        np.put(contig_methylation['coverage'], part_C_incontext_pos,contig_methylation['coverage'][part_C_incontext_pos] + 1)
        # for methy record
        if methylated_num != 0 :
            # 02 CHALM methy record
            if methylated_num >= CHALM_threshold:
                CHALM_methylated_pos = part_C_incontext_pos
            else:
                CHALM_methylated_pos = part_methylated_pos
        if part_methylated_num != 0 :
            # 01 trad mean methy record
            np.put(contig_methylation['meth'],part_methylated_pos,contig_methylation['meth'][part_methylated_pos] + 1)
            # 02 CHALM methy record
            np.put(contig_methylation['CHALM_meth'],CHALM_methylated_pos,contig_methylation['CHALM_meth'][CHALM_methylated_pos] + 1)       
        # 03 CAMDA methy record
        if methylated_num != totalC_num and part_methylated_num != part_totalC_num:
            CAMDA_unmethylated_pos = [ start_pos + index for index in range(scan_start,scan_end) if read_seq[index].islower() ]
            np.put(contig_methylation['CAMDA_meth'],CAMDA_unmethylated_pos,contig_methylation['CAMDA_meth'][CAMDA_unmethylated_pos] + 1)
    return contig_methylation
############################################################################