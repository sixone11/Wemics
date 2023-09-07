#!/bin/env python3
import os,sys
import numpy as np
import pandas as pd
import math
from itertools import groupby,chain
# import pysam
import argparse
from utils import *
from function import *
#####################################
def cal_add_number(row):
    if(row["ele"]!=0):
        return [row["count"]*i for i in row["count"]*[len(seq)]]
    else:
        return row["count"]*[len(seq)]
#####################################   
contig_methylation = {}
contig_length = 250000000
#####################################
# debug variable
infile = sys.argv[1]
contig = sys.argv[2]
outfile = sys.argv[3]
###########################################################################
data = pd.read_csv(infile,sep="\t",header=None)
data.columns = ["chr","start","end","strand","seq","C_pos"]
data_contig = data[data["chr"] == contig]
data_contig.reset_index(drop=True,inplace=True)
###########################################################################
# Initialization
for element in ['coverage','CHALM_coverage','CHALM_meth']:  
    contig_methylation[element] = np.zeros(contig_length, dtype=int)
for element in ['strand']:
    contig_methylation[element] = np.chararray(contig_length)
###########################################################################
for index1 in range(0,data_contig.shape[0]):
    if index1 % 10000 == 0:
        disp("hit 10000 times")
    seq = np.asarray(eval(data_contig["seq"][index1]))
    strand = data_contig["strand"][index1]
    pos_list = data_contig["C_pos"][index1]
    pos_list = eval(pos_list)
    C_incontext_pos = pos_list
    if 0 in set(seq) and len(set(seq)) == 1:
        np.put(contig_methylation['CHALM_coverage'],C_incontext_pos,contig_methylation['CHALM_coverage'][C_incontext_pos] + 1)
    else:
        methylated_index = [ i for i in range(0,len(seq)) if seq[i] == 1 ]
        consecutive_ele_count  = pd.DataFrame(np.array([[key,sum(1 for ele in group)] for key,group in groupby(seq)]),columns=["ele","count"])
        add_number = consecutive_ele_count.apply(cal_add_number,axis=1)
        add_number = np.array(list(chain(*add_number)))
        np.put(contig_methylation['CHALM_coverage'],C_incontext_pos,contig_methylation['CHALM_coverage'][C_incontext_pos] + add_number)
        np.put(contig_methylation['CHALM_meth'],C_incontext_pos,contig_methylation['CHALM_meth'][C_incontext_pos] + add_number)
    np.put(contig_methylation['coverage'],C_incontext_pos,contig_methylation['coverage'][C_incontext_pos] + 1)
    np.put(contig_methylation['strand'],C_incontext_pos,strand)
C_incontext_pos = (np.where(contig_methylation['coverage']!= 0)[0]).tolist()
# the length of C_incontext_pos could be 0!
out = pd.DataFrame({"chr":contig,
              "start":C_incontext_pos,
              "end":np.array(C_incontext_pos) + 1,
              "depth":contig_methylation['coverage'][C_incontext_pos],
              "CHALM_weighted":contig_methylation['CHALM_meth'][C_incontext_pos],
              "CHALM_coverage":contig_methylation['CHALM_coverage'][C_incontext_pos],
              "strand":contig_methylation['strand'][C_incontext_pos]})
out['CHALM_weighted_ratio'] = out['CHALM_weighted']/out['CHALM_coverage']
out.drop("depth",axis=1,inplace=True)
out["strand"] = out["strand"].apply(lambda x:x.decode("utf-8"))
out.to_csv(outfile,sep='\t',encoding='gbk',header=None,index=False)
###########################################################################
