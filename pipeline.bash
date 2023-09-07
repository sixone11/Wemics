#!/bin/bash
# modifying by yourself
wkd=/data/sixone/lllab/RRBS/downstream_analysis
bam=${wkd}/Wemics_example/example.bam
sample=example 

# Work directory 
cd ${wkd}/Wemics_example

# Process the bam file from bismark(bam file was sorted by name)
# Separate the information of read1 and read2 from bam file
samtools view -@10 ${sample}.bam|\
awk -v read1=${sample}_read1.out  -v read2=${sample}_read2.out '{if(NR%2==1) print $0 >>read1;else print $0 >>read2}' 

for read in read1 read2 
do 
cat ${sample}_${read}.out|awk 'BEGIN{FS=OFS="\t"} {gsub("XM:Z:","",$14);if($9>0) $9="+";else $9="-"; print $3,$4,$9,$6,$14}'> ${sample}_${read}.bed &
done 

wait 

paste -d "\t" ${sample}_read1.bed ${sample}_read2.bed|grep ^chr[1-9]|sort -k1,1V -k2,2n|awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$7,$8,$9,$10}'> ${sample}_reads12.bed

# Using chr21 as an example
contig=chr21 
cat ${sample}_reads12.bed |grep ${contig} > ${sample}_${contig}.bed 

# Process reads
python3 ${wkd}/Wemics_example/bin/process_reads.py \
	${sample}_${contig}.bed  \
	${contig} \
	${sample}_${contig}_processed.bed

# Sort reads
cat ${sample}_${contig}_processed.bed|sort -k1,1V -k2,2n -k3,3n > ${sample}_readmeth.bed 


# Wemics qutification
python3 ${wkd}/Wemics_example/bin/Wemics.py \
	    ${sample}_readmeth.bed \
	    ${contig} \
	    ${sample}_${contig}_Wemics.bed 
		

# Merge all contig file into one
# cat $(for contig in chr{1..22};do echo "${sample}_${contig}_Wemics.bed ";done) \
# |awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$5,$4,$7}' > ${sample}_Wemics.bed

# Because we only have chr21 in our bam file, we could just use command below
cat ${sample}_${contig}_Wemics.bed|awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$5,$4,$7}' > ${sample}_Wemics.bed


# Remove temp file
rm ${sample}_${contig}_Wemics.bed \
${sample}_read*.out ${sample}_read*.bed

