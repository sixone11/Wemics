# Wemics
The code for the manuscript "Wemics: a single-base resolution methylation quantification method by weighting methylation of consecutive CpG sites".

# Introduction
Proper quantification of RRBS data is the most important first step in studying DNA methylation. The methylation level of a genomic region has long been quantified by C/C+T(i.e., the average proportion of methylated cytosines in a genomic region).

A good DNA methylation quantification method should reflect a strong negative correlation between gene promoter methylation level and its gene expression, thereby increasing our understanding of epigenetic mechanism by which DNA methylation regulates gene expression.

In this repository, we introduced a single-base resolution methylation quantification method by weighting methylation of consecutive CpG sites (Wemics) in genomic regions. Wemics reveals cell heterogeneity in DNA methylation by distinguishing methylated reads, unmethylated reads and paritially methylated reads in bulk sequenced cells. It provides reliable quantification of DNA methylation levels by borrowing information from consecutive CpG sites as there is a strong dependency of DNA methylation among consecutive CpG sites.

# Code
This repository contains the scripts used to extract the methylation level from CpG sites under Wemics quantification.

Code contains:
- `pipeline.bash`
- `process_reads.py`
- `Wemics.py`
- `utils.py`
- `function.py`

You just need follow the command lines in `pipeline.bash` to get the `.bed` file under Wemics quantification.

The input is the bam file which was processed by bismark and sorted by read name (i.e.,` samtools sort -n` ). We also provide an example bam file so that you can better understand the Wemics.  Through the `pipeline.bash`, you could get a bed file, which contains 7 columns: chr, start, end, the Wemics methylated counts, the Wemics depth, and the methylation level under the quantification of Wemics.

# Dependencies
python3
- itertools
- math
- numpy
- pandas
- argparse

# Citation

If you find our research useful, please consider citing: Wemics: a single-base resolution methylation quantification method by weighting methylation of consecutive CpG sites(2023).

# Contact

yiliu11@zju.edu.cn
