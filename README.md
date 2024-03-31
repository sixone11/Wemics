# Wemics
Code for the reference "Wemics: A Single-Base Resolution Methylation Quantification Method for Enhanced Prediction of Epigenetic Regulation".

# Introduction
Proper quantification of RRBS data is the most important first step in studying DNA methylation. The methylation level of a genomic region has long been quantified by C/C+T(i.e., the average proportion of methylated cytosines in a genomic region).

A good DNA methylation quantification method should reflect a strong negative correlation between gene promoter methylation level and its gene expression, thereby increasing our understanding of epigenetic mechanism by which DNA methylation regulates gene expression.

In this repository, we introduced a single-base resolution methylation quantification method by weighting methylation of consecutive CpG sites (Wemics) in genomic regions. Wemics reveals cell heterogeneity in DNA methylation by distinguishing methylated reads, unmethylated reads and paritially methylated reads in bulk sequenced cells. It provides reliable quantification of DNA methylation levels by borrowing information from consecutive CpG sites as there is a strong dependency of DNA methylation among consecutive CpG sites.

# Code
This repository contains scripts for extracting methylation levels at CpG sites under Wemics quantification.

Code contains:
- `pipeline.bash`
- `process_reads.py`
- `Wemics.py`
- `utils.py`
- `function.py`

You only need to follow the command lines in pipeline.bash to get the output result file (i.e., bed file) under Wemics quantification.

The input is a bam file processed by bismark and sorted by read name(i.e.,` samtools sort -n` ). We also provide an example bam file to give a better understanding of Wemics. Through the `pipeline.bash`, you could get a bed file that contains 7 columns: chr, start, end, Wemics methylated counts, Wemics depth, and methylation level under Wemics quantification.

# Dependencies
python3
- itertools
- math
- numpy
- pandas
- argparse

# Citation

If you find our research useful, please consider citing: Liu Y, Yi J, Wu P, Zhang J, Li X, Li J, Zhou L, Liu Y, Xu H, Chen E, Zhang H, Liang M, Liu P, Pan X, Lu Y. Wemics: A Single-Base Resolution Methylation Quantification Method for Enhanced Prediction of Epigenetic Regulation. Adv Sci (Weinh). 2024 Mar 28:e2308884. doi: 10.1002/advs.202308884

# Contact

yiliu11@zju.edu.cn
